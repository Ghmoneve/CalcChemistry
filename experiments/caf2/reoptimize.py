"""Re-optimise every species to a verified minimum, then recompute the reaction.

WHY THE ORIGINAL RUN IS NOT ENOUGH
reaction2.py relaxed with a damped BFGS and stopped when the gradient norm fell
below 6e-4. That criterion is satisfied by saddle points as well as minima, and a
Hessian run (experiments/caf2/frequencies.py) found imaginary frequencies in three
of the four species:

    formamidine   -245.4 cm^-1
    CaF2          -103.5 cm^-1  (doubly degenerate bend)
    Ca_famd2       -74.5 cm^-1
    HF            clean

A reaction energy built on saddle points is not a reaction energy, so the number in
the paper has to be replaced rather than merely qualified.

WHAT THIS SCRIPT DOES
For each species: relax with a tighter gradient tolerance and a smaller step cap;
compute the analytic Hessian; project out translations and rotations; if any mode is
imaginary, displace the geometry along that mode in mass-weighted coordinates and
relax again. Repeat until no imaginary mode survives or the round budget is spent.
Every round is appended to JSON so a long run is inspectable while it proceeds.

CaF2 is started from a bent geometry because the linear one is a saddle: an explicit
angle scan puts the minimum near 140 degrees, 2.06 kcal/mol below linear at
PBE/def2-svp. That the computed minimum is bent while gas-phase CaF2 is conventionally
linear is a known difficulty with the very flat bending potential of the heavier
alkaline-earth dihalides; it is reported, not hidden.

WHAT THIS DOES NOT DO
No zero-point energy, thermal correction or free energy is computed. The reaction
energy reported here is a bare electronic energy difference between verified minima,
which is what the paper claims and nothing more.
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
OUT = ROOT / 'experiments/caf2/reoptimized.json'
LOG = ROOT / 'experiments/caf2/reoptimized.log'

XC, BASIS, GRID = 'pbe', 'def2-svp', 4
CONV = 5140.02
HARTREE2KCAL = 627.5094740631
MAX_ROUNDS = 3
TARGET_RMS = 0.12          # Angstrom, displacement used to leave a saddle


def say(msg: str) -> None:
    line = f'[{time.strftime("%H:%M:%S")}] {msg}'
    print(line, flush=True)
    with LOG.open('a') as fh:
        fh.write(line + '\n')


def dump(obj) -> None:
    OUT.write_text(json.dumps(obj, indent=2))


def build(zs, x):
    return gto.M(atom=[(zs[i], tuple(x[i])) for i in range(len(zs))],
                 basis=BASIS, spin=0, charge=0, unit='Angstrom', verbose=0)


def scf(mol):
    """SCF with retries; refuses to return an unconverged energy."""
    attempts = []
    for damping, init in ((0.0, None), (0.3, None), (0.0, 'atom')):
        mf = dft.RKS(mol)
        mf.xc = XC
        mf.grids.level = GRID
        mf.conv_tol = 1e-11
        if damping:
            mf.damp = damping
        if init:
            mf.init_guess = init
        try:
            mf.kernel()
        except Exception as exc:            # noqa: BLE001 - recorded below
            attempts.append(f'{damping}/{init}: {type(exc).__name__}')
            continue
        if mf.converged:
            return float(mf.e_tot), mf
        attempts.append(f'{damping}/{init}: not converged')
    raise RuntimeError(f'SCF failed for {mol.natm} atoms: {attempts}')


def relax(zs, x0, max_iter=200, gtol=1e-5, step_cap=0.08):
    x = np.array(x0, dtype=float).copy()
    n = len(zs)
    hinv = np.eye(3 * n) * 0.3
    gp = None
    xp = None
    e = 0.0
    gnorm = float('inf')
    for it in range(max_iter):
        mol = build(zs, x)
        e, mf = scf(mol)
        g = mf.nuc_grad_method().kernel().reshape(-1) / gto.param.BOHR
        gnorm = float(np.linalg.norm(g))
        if gnorm < gtol:
            break
        if gp is not None:
            s = (x - xp).reshape(-1)
            y = g - gp
            sy = float(s @ y)
            if sy > 1e-12:
                rho = 1.0 / sy
                idn = np.eye(3 * n)
                hinv = ((idn - rho * np.outer(s, y)) @ hinv
                        @ (idn - rho * np.outer(y, s)) + rho * np.outer(s, s))
        step = -hinv @ g
        ns = float(np.linalg.norm(step))
        if ns > step_cap:
            step *= step_cap / ns
        xp, gp = x.copy(), g.copy()
        x = x + step.reshape(n, 3)
    return e, x, gnorm, it + 1


def hess_modes(zs, x):
    """Return (frequencies in cm^-1, mass-weighted eigenvectors) for the physical modes."""
    mol = build(zs, x)
    _, mf = scf(mol)
    h = np.array(hrks.Hessian(mf).kernel())
    H = h.transpose(0, 2, 1, 3).reshape(3 * mol.natm, 3 * mol.natm)
    H = 0.5 * (H + H.T)
    masses = np.array(mol.atom_mass_list())
    sq = np.repeat(np.sqrt(masses), 3)
    Hm = H / np.outer(sq, sq)

    R = mol.atom_coords()
    N = mol.natm
    cands = []
    for a in range(3):
        v = np.zeros(3 * N)
        for i in range(N):
            v[3 * i + a] = sq[3 * i]
        cands.append(v)
    for a in range(3):
        e = np.zeros(3)
        e[a] = 1.0
        v = np.zeros(3 * N)
        for i in range(N):
            v[3 * i:3 * i + 3] = sq[3 * i] * np.cross(e, R[i])
        cands.append(v)
    bv = []
    for v in cands:
        nrm = np.linalg.norm(v)
        if nrm < 1e-8:
            continue
        v = v / nrm
        for u in bv:
            v = v - float(v @ u) * u
        nrm = np.linalg.norm(v)
        if nrm > 1e-6:
            bv.append(v / nrm)
    P = np.eye(3 * N) - sum(np.outer(u, u) for u in bv)
    Hp = P @ Hm @ P
    w, vecs = np.linalg.eigh(0.5 * (Hp + Hp.T))
    # the physical modes are the 3N - len(bv) eigenvalues of largest magnitude
    keep = np.argsort(-np.abs(w))[: 3 * N - len(bv)]
    keep = keep[np.argsort(w[keep])]
    freqs = np.sign(w[keep]) * CONV * np.sqrt(np.abs(w[keep]))
    return freqs, vecs[:, keep], sq


def displace(zs, x, vec_mw, sq, rms=TARGET_RMS):
    """Move along a mass-weighted eigenvector with a controlled Cartesian RMS."""
    n = len(zs)
    dx = (vec_mw / sq).reshape(n, 3)
    mag = float(np.sqrt((dx ** 2).sum() / n))
    if mag < 1e-12:
        return x
    dx *= rms / mag
    return np.array(x) + dx


def to_minimum(name, zs, x0, state):
    train = []
    x = np.array(x0, dtype=float)
    e, gnorm = 0.0, float('inf')
    freqs = None
    for rnd in range(1, MAX_ROUNDS + 1):
        t0 = time.time()
        e, x, gnorm, niter = relax(zs, x)
        freqs, vecs, sq = hess_modes(zs, x)
        imag = freqs[freqs < -1.0]
        train.append({'round': rnd, 'E': float(e), 'grad_norm': gnorm,
                      'niter': niter, 'imaginary_count': int(len(imag)),
                      'imaginary_cm': [float(v) for v in imag],
                      'lowest_real_cm': float(freqs[freqs >= -1.0].min())
                      if (freqs >= -1.0).any() else None,
                      'highest_cm': float(freqs.max()),
                      'coords': x.tolist(), 'sec': round(time.time() - t0, 1)})
        say(f'  {name} round {rnd}: E={e:.6f} |g|={gnorm:.2e} '
            f'imag={len(imag)} '
            f'{np.round(imag,1) if len(imag) else ""} '
            f'({time.time()-t0:.0f}s)')
        if len(imag) == 0:
            say(f'  {name}: MINIMUM reached')
            break
        if rnd == MAX_ROUNDS:
            say(f'  {name}: still has imaginary modes after {MAX_ROUNDS} rounds')
            break
        # displace along the most negative mode
        idx = int(np.argmin(freqs))
        x = displace(zs, x, vecs[:, idx], sq)
    state[name] = {'final': train[-1], 'rounds': train,
                   'is_minimum': train[-1]['imaginary_count'] == 0}
    dump(state)
    return np.array(train[-1]['coords']), train[-1]['imaginary_count'] == 0


def main():
    LOG.write_text('')
    src = json.loads(SRC.read_text())
    state = {}
    if OUT.exists():
        try:
            state = json.loads(OUT.read_text())
            state = {k: v for k, v in state.items() if not k.startswith('_')}
            say(f'resuming from {OUT.name}: {list(state)}')
        except Exception:
            state = {}
    coords = {}

    say(f're-optimising at {XC}/{BASIS}, gtol 1e-5, up to {MAX_ROUNDS} saddle escapes')

    # CaF2 starts bent: the linear geometry is a saddle and an angle scan puts the
    # minimum near 140 degrees.
    start = {}
    for name, rec in src['relaxed'].items():
        start[name] = np.array(rec['coords'], dtype=float)
    zs_ca = src['relaxed']['CaF2']['zs']
    r = 1.9925
    th = np.radians(140.0)
    start['CaF2'] = np.array([
        [0.0, 0.0, 0.0],
        [r * np.sin(th / 2), 0.0, r * np.cos(th / 2)],
        [-r * np.sin(th / 2), 0.0, r * np.cos(th / 2)],
    ])

    for name in ['formamidine', 'HF', 'CaF2', 'Ca_famd2']:
        zs = src['relaxed'][name]['zs']
        # resume: reuse a species already verified as a minimum in a previous run
        if name in state and state[name].get('is_minimum'):
            say(f'{name}: already a verified minimum, reusing')
            coords[name] = np.array(state[name]['final']['coords'], dtype=float)
            continue
        say(f'{name}: {len(zs)} atoms')
        x, ok = to_minimum(name, zs, start[name], state)
        coords[name] = x

    say('')
    ok_all = all(state[k]['is_minimum'] for k in state)
    if not ok_all:
        bad = [k for k in state if not state[k]['is_minimum']]
        say(f'WARNING: not all minima: {bad}')
    else:
        say('all species are verified minima')

    # reaction energy at the relaxed level
    E = {k: state[k]['final']['E'] for k in state}
    dE = E['CaF2'] + 2 * E['formamidine'] - E['Ca_famd2'] - 2 * E['HF']
    state['_reaction'] = {
        'equation': 'Ca(famd)2 + 2 HF -> CaF2 + 2 formamidine',
        'level': f'{XC}/{BASIS}',
        'E': E, 'dE_Ha': float(dE), 'dE_kcal': float(dE * HARTREE2KCAL),
        'all_minima': bool(ok_all),
    }
    dump(state)
    say(f'reaction dE = {dE*HARTREE2KCAL:+.2f} kcal/mol at {XC}/{BASIS}')

    # single points across the spread at the NEW geometries
    say('single points at the re-optimised geometries')
    spreads = {}
    for xc, basis in [('pbe', 'def2-svp'), ('pbe', 'def2-tzvp'),
                      ('b3lyp', 'def2-svp'), ('b3lyp', 'def2-tzvp')]:
        tag = f'{xc}/{basis}'
        E2 = {}
        for name, xco in coords.items():
            zs = src['relaxed'][name]['zs']
            mol = gto.M(atom=[(zs[i], tuple(xco[i])) for i in range(len(zs))],
                        basis=basis, spin=0, charge=0, unit='Angstrom', verbose=0)
            mfx = dft.RKS(mol)
            mfx.xc = xc
            mfx.grids.level = GRID
            mfx.conv_tol = 1e-11
            mfx.kernel()
            if not mfx.converged:
                raise RuntimeError(f'SCF failed {tag} {name}')
            E2[name] = float(mfx.e_tot)
        dE2 = (E2['CaF2'] + 2 * E2['formamidine'] - E2['Ca_famd2'] - 2 * E2['HF'])
        spreads[tag] = {'E': E2, 'dE_Ha': float(dE2),
                        'dE_kcal': float(dE2 * HARTREE2KCAL)}
        say(f'  {tag}: dE = {dE2*HARTREE2KCAL:+.2f} kcal/mol')
        state['_single_points'] = spreads
        dump(state)

    vals = [v['dE_kcal'] for v in spreads.values()]
    state['_spread'] = {'min': min(vals), 'max': max(vals),
                        'spread_kcal': max(vals) - min(vals)}
    dump(state)
    say(f'spread across four levels = {max(vals)-min(vals):.2f} kcal/mol')
    say(f'written {OUT}')


if __name__ == '__main__':
    main()
