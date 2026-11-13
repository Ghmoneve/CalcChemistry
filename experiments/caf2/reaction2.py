"""Gas-phase model reaction for CaF2 ALD ligand protonation.

    Ca(formamidinate)2 + 2 HF  ->  CaF2 + 2 formamidine

Electron balance: 66 + 2*10 = 86 = 38 + 2*24.

DESIGN (the first version of this run proved unobservable for 105 minutes):
  1. Relax each molecule ONCE at PBE/def2-svp. Geometry depends weakly on the
     functional, so re-relaxing at every combination multiplies cost for no
     information.
  2. Single-point every molecule at those geometries for each (functional, basis)
     pair, so the reported spread is purely electronic.
  3. Write JSON after every molecule and every single point, and log every BFGS
     iteration. A run that only writes at the end cannot be told apart from a hung
     one, which is what went wrong before.

MODELLING STATUS, stated here because it matters:
  * The real precursor is bis(N,N'-di-i-propylformamidinate)calcium(II); each N
    carries an i-propyl group. This model replaces every i-propyl by H.
    Substituent effects on a protonation energy are real and are NOT captured, so
    the number produced is a model number for the truncated ligand.
  * No geometry optimiser ships with this pyscf, so a damped BFGS on pyscf
    analytic nuclear gradients is used. Each iteration is one SCF plus one
    gradient; the final gradient norm is reported with every energy.
  * One functional and one basis cannot bracket the electronic answer, so a spread
    is reported instead of a single number presented as exact.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from pyscf import gto, dft

HARTREE2KCAL = 627.5094740631
OUT = Path('experiments/caf2/reaction2.json')
LOG = Path('experiments/caf2/reaction2.log')

RELAX_LEVEL = ('pbe', 'def2-svp')
SPREAD = [('pbe', 'def2-svp'), ('pbe', 'def2-tzvp'),
          ('b3lyp', 'def2-svp'), ('b3lyp', 'def2-tzvp')]

GEOMS = {
    'formamidine': ('C 0.00 0.00 0.00; N 0.00 1.35 0.00; N 0.00 -1.35 0.00; '
                    'H 1.05 0.00 0.00; H 0.00 2.30 0.00; '
                    'H 0.90 -1.90 0.00; H -0.90 -1.90 0.00'),
    'HF': 'H 0.00 0.00 0.00; F 0.00 0.00 0.917',
    'CaF2': 'Ca 0.00 0.00 0.00; F 0.00 0.00 2.10; F 0.00 0.00 -2.10',
    'Ca_famd2': ('Ca 0.00 0.00 0.00; '
                 'N 2.45 0.80 0.00; N 2.45 -0.80 0.00; C 3.60 0.00 0.00; '
                 'H 4.50 0.00 0.00; H 2.70 1.80 0.00; H 2.70 -1.80 0.00; '
                 'N -2.45 0.80 0.00; N -2.45 -0.80 0.00; C -3.60 0.00 0.00; '
                 'H -4.50 0.00 0.00; H -2.70 1.80 0.00; H -2.70 -1.80 0.00'),
}
ORDER = ['formamidine', 'HF', 'CaF2', 'Ca_famd2']


def say(msg):
    line = f'[{time.strftime("%H:%M:%S")}] {msg}'
    print(line, flush=True)
    with LOG.open('a') as fh:
        fh.write(line + '\n')


def dump(obj):
    OUT.write_text(json.dumps(obj, indent=2))


def scf(mol, xc):
    """SCF with an explicit convergence guard and two fallbacks.

    A plain `kernel()` occasionally fails to converge on a small molecule, and the
    failure is environment-dependent: thread count changes the summation order in
    the XC grid. Rather than trust one attempt, retry with damping and then from a
    Hartree-Fock density. Every attempt must report `converged`, otherwise no
    energy is returned at all.
    """
    from pyscf import scf as pscf

    attempts = []

    def run(kwargs, dm0=None, tag='default'):
        mf = dft.RKS(mol)
        mf.xc = xc
        mf.max_cycle = 300
        mf.conv_tol = 1e-9
        for k, v in kwargs.items():
            setattr(mf, k, v)
        e = mf.kernel(dm0=dm0) if dm0 is not None else mf.kernel()
        attempts.append((tag, bool(mf.converged), int(mf.cycles)))
        return (e, mf) if mf.converged else None

    out = run({})
    if out is not None:
        return out
    out = run({'damp': 0.3, 'diis_space': 12}, tag='damped')
    if out is not None:
        return out
    mf_h = pscf.RHF(mol)
    mf_h.conv_tol = 1e-9
    mf_h.max_cycle = 300
    mf_h.kernel()
    out = run({}, dm0=mf_h.make_rdm1(), tag='from_RHF')
    if out is not None:
        return out
    raise RuntimeError(
        f'SCF not converged ({xc}) after {len(attempts)} attempts {attempts}; '
        'refusing to use this energy')


def relax(name, spec, xc, basis, max_iter=60, gtol=6e-4):
    """Damped BFGS on nuclear coordinates, logging every iteration.

    UNIT TRAP: `mol.atom_coord` returns Bohr, while `gto.M(atom=...)` interprets
    bare numbers as Angstrom. Rebuilding the molecule from raw `atom_coord`
    therefore shrinks every bond by a factor of `BOHR` and the SCF stops
    converging -- which is what happened here before this was fixed. Coordinates
    are converted once on the way in and once on the way out.
    """
    mol0 = gto.M(atom=spec, basis=basis, spin=0, charge=0, verbose=0)
    n = mol0.natm
    zs = [int(mol0.atom_charge(i)) for i in range(n)]
    # Angstrom, which is what gto.M expects when we rebuild.
    x = np.array([mol0.atom_coord(i) for i in range(n)]) * gto.param.BOHR
    hinv = np.eye(3 * n) * 0.4
    gp = None
    xp = None
    e = 0.0
    gnorm = float('inf')
    niter = 0

    for it in range(max_iter):
        niter = it + 1
        mol = gto.M(atom=[(zs[i], tuple(x[i])) for i in range(n)],
                    basis=basis, spin=0, charge=0, verbose=0)
        e, mf = scf(mol, xc)
        # The nuclear gradient is dE/d(Bohr). Coordinates are held in Angstrom, so
        # convert the gradient to dE/d(Angstrom) before using it for a step;
        # otherwise every step is too large by a factor of BOHR.
        g = mf.nuc_grad_method().kernel().reshape(-1) / gto.param.BOHR
        gnorm = float(np.linalg.norm(g))
        say(f'  {name} iter {it:2d}  E={e:.6f}  |g|={gnorm:.2e}')
        if gnorm < gtol:
            break
        if gp is not None:
            s = (x - xp).reshape(-1)
            y = g - gp
            sy = float(s @ y)
            if sy > 1e-12:
                rho = 1.0 / sy
                ident = np.eye(3 * n)
                hinv = ((ident - rho * np.outer(s, y)) @ hinv
                        @ (ident - rho * np.outer(y, s)) + rho * np.outer(s, s))
        step = -hinv @ g
        ns = float(np.linalg.norm(step))
        if ns > 0.20:
            step *= 0.20 / ns
        xp, gp = x.copy(), g.copy()
        x = x + step.reshape(n, 3)

    return e, x, gnorm, zs, niter

def main():
    LOG.write_text('')
    state = {'relaxed': {}, 'single_points': {}, 'reactions': {}}
    dump(state)

    say(f'STEP 1: relax each molecule once at {RELAX_LEVEL[0]}/{RELAX_LEVEL[1]}')
    for name in ORDER:
        t0 = time.time()
        say(f'relaxing {name} ...')
        e, x, gn, zs, niter = relax(name, GEOMS[name], *RELAX_LEVEL)
        state['relaxed'][name] = {
            'E': float(e), 'grad_norm': gn, 'converged_grad': gn < 6e-4,
            'coords': x.tolist(), 'zs': zs, 'niter': niter,
            'level': list(RELAX_LEVEL), 'sec': round(time.time() - t0, 1)}
        dump(state)
        say(f'  {name}: E={e:.6f} Ha  |g|={gn:.2e}  ({time.time()-t0:.0f}s)')

    say('STEP 2: single points across the method spread')
    for xc, basis in SPREAD:
        tag = f'{xc}/{basis}'
        state['single_points'][tag] = {}
        for name in ORDER:
            t0 = time.time()
            rec = state['relaxed'][name]
            xyz, zs = rec['coords'], rec['zs']
            mol = gto.M(atom=[(zs[i], tuple(xyz[i])) for i in range(len(zs))],
                        basis=basis, spin=0, charge=0, verbose=0)
            e, _ = scf(mol, xc)
            state['single_points'][tag][name] = {
                'E': float(e), 'sec': round(time.time() - t0, 1)}
            dump(state)
            say(f'  {tag} {name}: E={e:.6f} Ha ({time.time()-t0:.0f}s)')
        sp = state['single_points'][tag]
        dE = (sp['CaF2']['E'] + 2 * sp['formamidine']['E']
              - sp['Ca_famd2']['E'] - 2 * sp['HF']['E'])
        state['reactions'][tag] = {'dE_Ha': dE, 'dE_kcal': dE * HARTREE2KCAL}
        dump(state)
        say(f'  {tag} REACTION dE = {dE:+.6f} Ha = {dE*HARTREE2KCAL:+.2f} kcal/mol')

    vals = [r['dE_kcal'] for r in state['reactions'].values()]
    say('')
    say('SUMMARY  dE over the method spread (kcal/mol):')
    for tag, r in state['reactions'].items():
        flag = '' if state['relaxed'] else ''
        say(f'  {tag:22s} {r["dE_kcal"]:+8.2f}{flag}')
    say(f'  spread = {max(vals)-min(vals):.2f} kcal/mol over {len(vals)} methods')
    say('NOTE: geometry relaxed at PBE/def2-svp only; truncated ligand (i-Pr -> H).')


if __name__ == '__main__':
    main()
