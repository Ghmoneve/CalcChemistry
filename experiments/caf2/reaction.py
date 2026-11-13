"""Gas-phase model reaction for CaF2 ALD ligand protonation.

    Ca(formamidinate)2 + 2 HF  ->  CaF2 + 2 formamidine

Electron balance: 66 + 2*10 = 86 = 38 + 2*24.

MODELLING STATUS, stated plainly because it matters:
  * The real precursor is bis(N,N'-di-i-propylformamidinate)calcium(II): the
    amidinate carries two i-propyl groups on each N. This model replaces every
    i-propyl by H. Substituent effects on a protonation energy are real and are
    NOT captured, so the number produced is a model number for the truncated
    ligand, not for the ALD precursor.
  * No geometry optimiser ships with this pyscf installation, so a small BFGS
    minimiser is implemented here on pyscf analytic nuclear gradients.
  * A single functional and basis cannot bracket the electronic answer. The same
    reaction is run for several functionals and bases, and the spread is
    reported instead of one number presented as exact.
"""
from __future__ import annotations

import json
import time

import numpy as np
from pyscf import gto, dft

HARTREE2KCAL = 627.5094740631
FUNCTIONALS = ['pbe', 'b3lyp']
BASES = ['def2-svp', 'def2-tzvp']

GEOMS = {
    # 24 electrons, C1H4N2, planar HC(=NH)(NH2)
    'formamidine': ('C 0.00 0.00 0.00; N 0.00 1.35 0.00; N 0.00 -1.35 0.00; '
                    'H 1.05 0.00 0.00; H 0.00 2.30 0.00; '
                    'H 0.90 -1.90 0.00; H -0.90 -1.90 0.00'),
    # 10 electrons
    'HF': 'H 0.00 0.00 0.00; F 0.00 0.00 0.917',
    # 38 electrons, linear
    'CaF2': 'Ca 0.00 0.00 0.00; F 0.00 0.00 2.10; F 0.00 0.00 -2.10',
    # 66 electrons: Ca with two bidentate formamidinate ligands
    'Ca_famd2': ('Ca 0.00 0.00 0.00; '
                 'N 2.45 0.80 0.00; N 2.45 -0.80 0.00; C 3.60 0.00 0.00; '
                 'H 4.50 0.00 0.00; H 2.70 1.80 0.00; H 2.70 -1.80 0.00; '
                 'N -2.45 0.80 0.00; N -2.45 -0.80 0.00; C -3.60 0.00 0.00; '
                 'H -4.50 0.00 0.00; H -2.70 1.80 0.00; H -2.70 -1.80 0.00'),
}


def make_mol(spec, basis):
    return gto.M(atom=spec, basis=basis, spin=0, charge=0, verbose=0)


def scf_energy(mol, xc):
    mf = dft.RKS(mol)
    mf.xc = xc
    mf.max_cycle = 300
    mf.conv_tol = 1e-9
    e = mf.kernel()
    return e, mf


def relax(spec, xc, basis, max_iter=50, gtol=4e-4):
    """BFGS minimisation of the electronic energy over nuclear coordinates."""
    mol0 = make_mol(spec, basis)
    n = mol0.natm
    zs = [mol0.atom_charge(i) for i in range(n)]
    x = np.array([mol0.atom_coord(i) for i in range(n)])
    hinv = np.eye(3 * n) * 0.5
    grad_prev = None
    x_prev = None
    e = None
    converged_grad = None

    for _ in range(max_iter):
        mol = gto.M(atom=[(zs[i], tuple(x[i])) for i in range(n)],
                    basis=basis, spin=0, charge=0, verbose=0)
        e, mf = scf_energy(mol, xc)
        g = mf.nuc_grad_method().kernel().reshape(-1)
        gn = float(np.linalg.norm(g))
        converged_grad = gn
        if gn < gtol:
            break
        if grad_prev is not None:
            s = (x - x_prev).reshape(-1)
            y = g - grad_prev
            sy = float(s @ y)
            if sy > 1e-14:
                rho = 1.0 / sy
                ident = np.eye(3 * n)
                hinv = ((ident - rho * np.outer(s, y)) @ hinv
                        @ (ident - rho * np.outer(y, s))
                        + rho * np.outer(s, s))
        step = -hinv @ g
        nstep = float(np.linalg.norm(step))
        if nstep > 0.25:
            step = step * (0.25 / nstep)
        x_prev = x.copy()
        grad_prev = g.copy()
        x = x + step.reshape(n, 3)

    mol = gto.M(atom=[(zs[i], tuple(x[i])) for i in range(n)],
                basis=basis, spin=0, charge=0, verbose=0)
    e, mf = scf_energy(mol, xc)
    return e, x, converged_grad, bool(mf.converged)


def main():
    out = {}
    for basis in BASES:
        for xc in FUNCTIONALS:
            t0 = time.time()
            res = {}
            for name in GEOMS:
                try:
                    e, x, gnorm, conv = relax(GEOMS[name], xc, basis)
                    res[name] = {'E': e, 'grad_norm': gnorm,
                                 'scf_converged': conv, 'coords': x.tolist()}
                except Exception as exc:
                    res[name] = {'error': f'{type(exc).__name__}: {str(exc)[:90]}'}
            if all('error' not in v for v in res.values()):
                dE = (res['CaF2']['E'] + 2 * res['formamidine']['E']
                      - res['Ca_famd2']['E'] - 2 * res['HF']['E'])
                res['reaction'] = {'dE_Ha': dE, 'dE_kcal': dE * HARTREE2KCAL}
                worst = max(v['grad_norm'] for k, v in res.items()
                            if k != 'reaction')
                print(f'{xc:6s}/{basis:10s} dE = {dE:+.6f} Ha = '
                      f'{dE * HARTREE2KCAL:+8.2f} kcal/mol  '
                      f'max|grad| = {worst:.1e}  [{time.time()-t0:.0f}s]',
                      flush=True)
            else:
                bad = {k: v['error'] for k, v in res.items() if 'error' in v}
                print(f'{xc:6s}/{basis:10s} FAILED {bad}', flush=True)
            out[f'{xc}/{basis}'] = res
            json.dump(out, open('experiments/caf2/reaction.json', 'w'), indent=2)
    print('\nwritten experiments/caf2/reaction.json')


if __name__ == '__main__':
    main()
