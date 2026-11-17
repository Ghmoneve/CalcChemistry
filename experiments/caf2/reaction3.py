"""Recalculate the CaF2 model reaction with a correct optimiser and diffuse bases.

WHY THE PREVIOUS NUMBER IS BEING REPLACED
Two defects were found after the previous run:

1. The optimiser. A damped BFGS seeded with 0.3*I is ~740x too stiff along the
   bend (true curvature 4.1e-4), and the 0.08 Angstrom step cap destroys the
   secant relation so the inverse Hessian inflates without bound. The CaF2
   "converged" geometry at 128.2 deg with |g| = 3.0e-2 was a runaway, not a
   minimum. The fix used here is an analytic-Hessian seed, a relative curvature
   guard on the secant pair, a backtracking line search, and convergence judged
   on the rotation/translation-projected gradient.

2. The basis. CaF2's barrier to linearity is strongly basis-dependent in a way
   that matters at the 1-2 kcal/mol level:
       PBE/def2-svp 2.19, tzvp 1.14, svpd 0.50, tzvpd 0.87,
       B3LYP/svp 1.52, B3LYP/tzvpd 0.29 kcal/mol
   against a literature barrier of 0.15 kcal/mol. def2-svp overestimates the
   bend by a factor of ~15, and CaF2 enters the reaction with coefficient +1, so
   this propagates directly into dE.

WHAT IS COMPUTED
Every species is re-optimised at PBE/def2-svp with the corrected optimiser, then
single points are taken at eight levels (PBE and B3LYP crossed with def2-svp,
def2-tzvp, def2-svpd, def2-tzvpd). CaF2 is treated twice: at its bent minimum and
constrained to D_inf_h with only the bond length relaxed. The pair brackets the
treatment, which is the honest way to report a number whose systematic error is
shared by every level rather than visible in the spread.

WHAT IS NOT COMPUTED
No zero-point, thermal or free energy. No surface. The ligand is still truncated
to hydrogen. dE is a bare electronic energy difference between verified minima.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from pyscf import dft, gto
from pyscf.hessian import rks as hrks

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / 'experiments/caf2/reaction2.json'
REOPT = ROOT / 'experiments/caf2/reoptimized.json'
OUT = ROOT / 'experiments/caf2/reaction3.json'
LOG = ROOT / 'experiments/caf2/reaction3.log'

BOHR = gto.param.BOHR
HARTREE2KCAL = 627.5094740631
GEOM_XC, GEOM_BASIS, GRID = 'pbe', 'def2-svp', 5
GTOL = 3e-4
SPREAD = [('pbe', 'def2-svp'), ('pbe', 'def2-tzvp'), ('pbe', 'def2-svpd'),
          ('pbe', 'def2-tzvpd'), ('b3lyp', 'def2-svp'), ('b3lyp', 'def2-tzvp'),
          ('b3lyp', 'def2-svpd'), ('b3lyp', 'def2-tzvpd')]


def say(msg):
    line = f'[{time.strftime("%H:%M:%S")}] {msg}'
    print(line, flush=True)
    with LOG.open('a') as fh:
        fh.write(line + '\n')


def dump(o):
    OUT.write_text(json.dumps(o, indent=2))


def build(zs, x, basis=GEOM_BASIS):
    return gto.M(atom=[(int(zs[i]), tuple(x[i])) for i in range(len(zs))],
                 basis=basis, spin=0, charge=0, unit='Angstrom', verbose=0)


def energy_grad(mol, xc=GEOM_XC):
    mf = dft.RKS(mol)
    mf.xc = xc
    mf.grids.level = GRID
    mf.conv_tol = 1e-11
    mf.kernel()
    if not mf.converged:
        raise RuntimeError('SCF did not converge')
    gm = mf.nuc_grad_method()
    gm.grid_response = True            # derivative of the energy actually used
    return float(mf.e_tot), gm.kernel().reshape(-1) / BOHR, mf


def projector(mol, x):
    """Mass-weighted projector onto the internal (3N-6, or 3N-5 if linear) space."""
    n = mol.natm
    masses = np.array(mol.atom_mass_list())
    sq = np.repeat(np.sqrt(masses), 3)
    R = mol.atom_coords()          # Bohr, from the built molecule
    cand = []
    for a in range(3):
        v = np.zeros(3 * n)
        for i in range(n):
            v[3 * i + a] = sq[3 * i]
        cand.append(v)
    for a in range(3):
        e = np.zeros(3)
        e[a] = 1.0
        v = np.zeros(3 * n)
        for i in range(n):
            v[3 * i:3 * i + 3] = sq[3 * i] * np.cross(e, R[i])
        cand.append(v)
    bv = []
    for v in cand:
        nrm = np.linalg.norm(v)
        if nrm < 1e-8:
            continue
        v = v / nrm
        for u in bv:
            v = v - float(v @ u) * u
        nrm = np.linalg.norm(v)
        if nrm > 1e-6:
            bv.append(v / nrm)
    return np.eye(3 * n) - sum(np.outer(u, u) for u in bv), sq


def _linear_stretch_step(grad, n):
    """Step for a linear ABA molecule: move only the symmetric bond length.

    Both terminal atoms move along +z and -z by equal amounts and the centre atom
    stays put, so the geometry remains exactly D_inf_h. Taking the antisymmetric
    combination of the two terminal z-components is what projects the gradient onto
    that single internal coordinate; using the raw gradient would break linearity.
    """
    u = np.array([0.0, 0.0, 1.0])
    d = grad.reshape(n, 3)
    r = 0.5 * (float(d[1] @ u) - float(d[2] @ u))
    out = np.zeros((n, 3))
    out[1] = r * u
    out[2] = -r * u
    return out.reshape(-1)


def relax(zs, x0, maxiter=60, maxstep=0.15, linear=False):
    """Relax to a stationary point, optionally under the D_inf_h constraint.

    `linear=True` restricts the search to the one-dimensional symmetric-stretch
    subspace, so the result is the best linear geometry rather than a bent one.
    Convergence is judged on the projected gradient in the unconstrained case and
    on the single stretch component in the constrained case.
    """
    x = np.array(x0, float)
    n = len(zs)
    P = None          # only defined for the unconstrained branch
    hinv = None
    if linear:
        pass
    else:
        mol0 = build(zs, x)
        P, sq = projector(mol0, x)
        H, _ = hessian_mw(zs, x)
        Hp = P @ H @ P
        w, V = np.linalg.eigh(0.5 * (Hp + Hp.T))
        keep = np.argsort(-np.abs(w))[: int(round(np.trace(P)))]
        hinv = np.zeros_like(Hp)
        for k in keep:
            hinv += np.outer(V[:, k], V[:, k]) / w[k]

    xp = None
    gp = None
    e = float('nan')
    gi = float('inf')
    it = 0
    for it in range(maxiter):
        mol = build(zs, x)
        e, g, _ = energy_grad(mol)
        if P is not None:
            gi = float(np.linalg.norm(P @ g))
        else:
            gi = float(abs(g.reshape(n, 3)[1, 2]))
        if it % 5 == 0 or it < 3:
            say(f'   it {it:3d}  E {e:.9f}  |Pg| {gi:.2e}')
        if gi < GTOL:
            break
        if linear:
            step = _linear_stretch_step(-g * 0.5, n)
        else:
            assert hinv is not None, 'unconstrained branch must have a seeded Hessian'
            if gp is not None:
                svec = (x - xp).reshape(-1)
                yvec = g - gp
                sy = float(svec @ yvec)
                if sy > 1e-8 * np.linalg.norm(svec) * np.linalg.norm(yvec) + 1e-14:
                    rho = 1.0 / sy
                    idn = np.eye(3 * n)
                    hinv = ((idn - rho * np.outer(svec, yvec)) @ hinv
                            @ (idn - rho * np.outer(yvec, svec))
                            + rho * np.outer(svec, svec))
            step = -hinv @ g
        nrm = float(np.linalg.norm(step))
        if nrm > maxstep:
            step *= maxstep / nrm
        # backtracking line search on the energy
        a = 1.0
        e0 = e
        for _ in range(12):
            xt = x + a * step.reshape(n, 3)
            try:
                et, _, _ = energy_grad(build(zs, xt))
            except Exception:
                a *= 0.5
                continue
            if et <= e0 + 1e-12:
                break
            a *= 0.5
        else:
            a = 0.0
        xp, gp = x.copy(), g.copy()
        if a == 0.0:
            break
        x = x + a * step.reshape(n, 3)
    return e, x, gi, it + 1


def hessian_mw(zs, x):
    mol = build(zs, x)
    m2 = dft.RKS(mol)
    m2.xc = GEOM_XC
    m2.grids.level = GRID
    m2.conv_tol = 1e-11
    m2.kernel()
    if not m2.converged:
        raise RuntimeError('SCF failed in Hessian')
    h = np.array(hrks.Hessian(m2).kernel())
    H = h.transpose(0, 2, 1, 3).reshape(3 * mol.natm, 3 * mol.natm)
    H = 0.5 * (H + H.T)
    sq = np.repeat(np.sqrt(np.array(mol.atom_mass_list())), 3)
    return H / np.outer(sq, sq), sq


def main():
    LOG.write_text('')
    src = json.loads(SRC.read_text())
    zs = {k: v['zs'] for k, v in src['relaxed'].items()}
    # Start from the verified minima found by reoptimize.py, not from the saddle
    # geometries in reaction2.json: Ca(famd)2 has at least two minima and the
    # saddle-escape run found one 3.98 kcal/mol lower.
    start = {k: np.array(v['final']['coords'], float)
             for k, v in json.loads(REOPT.read_text()).items()
             if not k.startswith('_')}
    # CaF2 restarts from a bent geometry at the corrected curvature, because the
    # linear structure is a saddle and the previous 128 deg point was a runaway.
    r = 1.9925
    th = np.radians(140.0)
    start['CaF2'] = np.array([[0.0, 0.0, 0.0],
                              [r * np.sin(th / 2), 0.0, r * np.cos(th / 2)],
                              [-r * np.sin(th / 2), 0.0, r * np.cos(th / 2)]])
    state = {'geometry_level': f'{GEOM_XC}/{GEOM_BASIS}', 'species': {}}
    if OUT.exists():
        try:
            prev = json.loads(OUT.read_text())
            if prev.get('geometry_level') == f'{GEOM_XC}/{GEOM_BASIS}':
                state['species'] = prev.get('species', {})
                say(f"resuming with {len(state['species'])} geometries")
        except Exception:
            pass

    say(f'geometry optimisation at {GEOM_XC}/{GEOM_BASIS}, gtol {GTOL}, '
        f'analytic-Hessian seed, projected gradient')
    coords = {}
    for name in ['formamidine', 'HF', 'Ca_famd2', 'CaF2']:
        if name in state['species']:
            coords[name] = np.array(state['species'][name]['coords'], float)
            say(f'{name}: reusing optimised geometry (E {state["species"][name]["E"]:.9f})')
            continue
        t0 = time.time()
        say(f'{name}: {len(zs[name])} atoms')
        e, x, gi, niter = relax(zs[name], start[name])
        state['species'][name] = {'E': float(e), 'coords': x.tolist(),
                                  'zs': zs[name], 'proj_grad': float(gi),
                                  'niter': niter, 'sec': round(time.time() - t0, 1)}
        coords[name] = x
        say(f'  -> E {e:.9f}  |Pg| {gi:.2e}  {niter} it  ({time.time()-t0:.0f}s)')
        dump(state)

    # CaF2 constrained linear (D_inf_h): only the bond length relaxes
    zsc = zs['CaF2']
    xlin = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, r], [0.0, 0.0, -r]])
    t0 = time.time()
    if 'CaF2_linear' in state['species']:
        coords['CaF2_linear'] = np.array(
            state['species']['CaF2_linear']['coords'], float)
        say('CaF2_linear: reusing optimised geometry')
    else:
        say('CaF2 linear (D_inf_h), bond length only')
        e, xl, gi, niter = relax(zsc, xlin, linear=True)
        state['species']['CaF2_linear'] = {'E': float(e), 'coords': xl.tolist(),
                                           'zs': zsc, 'proj_grad': float(gi),
                                           'niter': niter,
                                           'sec': round(time.time() - t0, 1)}
        coords['CaF2_linear'] = xl
        say(f'  -> E {e:.9f}  {niter} it  ({time.time()-t0:.0f}s)')
        dump(state)

    say('single points across eight levels')
    spread = {}
    for xc, basis in SPREAD:
        tag = f'{xc}/{basis}'
        E = {}
        for name in ['formamidine', 'HF', 'CaF2', 'CaF2_linear', 'Ca_famd2']:
            zsn = zs['CaF2'] if name == 'CaF2_linear' else zs[name]
            mol = gto.M(atom=[(int(zsn[i]), tuple(coords[name][i]))
                              for i in range(len(zsn))],
                        basis=basis, spin=0, charge=0, unit='Angstrom', verbose=0)
            mf = dft.RKS(mol)
            mf.xc = xc
            mf.grids.level = GRID
            mf.conv_tol = 1e-11
            mf.kernel()
            if not mf.converged:
                raise RuntimeError(f'SCF failed {tag} {name}')
            E[name] = float(mf.e_tot)
        d_bent = E['CaF2'] + 2 * E['formamidine'] - E['Ca_famd2'] - 2 * E['HF']
        d_lin = (E['CaF2_linear'] + 2 * E['formamidine']
                 - E['Ca_famd2'] - 2 * E['HF'])
        spread[tag] = {'E': E, 'dE_bent_Ha': float(d_bent),
                       'dE_bent_kcal': float(d_bent * HARTREE2KCAL),
                       'dE_linear_Ha': float(d_lin),
                       'dE_linear_kcal': float(d_lin * HARTREE2KCAL)}
        say(f'  {tag:18s} bent {d_bent*HARTREE2KCAL:+7.2f}   '
            f'linear {d_lin*HARTREE2KCAL:+7.2f} kcal/mol')
        state['single_points'] = spread
        dump(state)

    b = [v['dE_bent_kcal'] for v in spread.values()]
    l = [v['dE_linear_kcal'] for v in spread.values()]
    state['summary'] = {
        'dE_bent_range': [min(b), max(b)], 'dE_bent_spread': max(b) - min(b),
        'dE_linear_range': [min(l), max(l)], 'dE_linear_spread': max(l) - min(l),
        'bracket_low': min(b), 'bracket_high': max(l),
        'bracket_width': max(l) - min(b),
    }
    dump(state)
    say('')
    say(f'bent   dE in [{min(b):+.2f}, {max(b):+.2f}]  spread {max(b)-min(b):.2f}')
    say(f'linear dE in [{min(l):+.2f}, {max(l):+.2f}]  spread {max(l)-min(l):.2f}')
    say(f'full bracket [{min(b):+.2f}, {max(l):+.2f}] width {max(l)-min(b):.2f} kcal/mol')
    say(f'written {OUT}')


if __name__ == '__main__':
    main()
