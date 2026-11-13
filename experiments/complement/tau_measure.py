"""Measure tau: the off-diagonal bound the complement certificate needs.

The certificate (CalcChemistry/ComplementBound.lean) needs both
    D_i >= d                       for the diagonal,
    |Q_O(x)| <= tau * ||x||^2      for the off-diagonal part,
to conclude lambda_min(H_C) >= d - tau. Direction 2 established that `d` can be
the reference determinant's own energy. This script measures `tau`.

WHY PYSCF'S FCI AND NOT HAND-ROLLED SLATER-CONDON
A first version wrote the determinant-basis Hamiltonian by hand and reported a
ground state of -108.86 Ha against pyscf's -75.01 Ha, wrong by 33.8 Ha: the sign
conventions for double excitations and the one-sided fill were both incorrect. A
hand-built FCI matrix is not worth re-deriving under time pressure when a verified
implementation exists, so the matrix is now obtained by applying pyscf's own
`contract_2e` to each unit vector, and checked against `fci.FCI` with an abort on
disagreement.

WHAT IS MEASURED
  * d, the diagonal minimum, and the gap to the next determinant;
  * the true spectral norm ||O||_2 of the off-diagonal part (exact here);
  * cheap proxies for it (Gershgorin max-row-abs-sum, max row 2-norm, Frobenius,
    max row nonzeros), so the distance between a cheap estimate and the truth is
    visible;
  * the trend of tau as the complement shrinks, which is what decides whether a
    complement certificate is decisive or vacuous.

SCOPE  A model measurement on water/STO-3G, whose 441 determinants are few enough
to diagonalise exactly. It says what tau IS where the truth is knowable. It
establishes nothing about an untruncated space, and no scaling law is inferred.
"""
from __future__ import annotations

import json
import sys
from math import comb
from pathlib import Path

import numpy as np
from pyscf import ao2mo, fci, gto, scf

ATOM = 'O 0 0 0.117; H 0 0.757 -0.467; H 0 -0.757 -0.467'
BASIS = 'sto-3g'
OUT = Path('experiments/complement/tau_measure.json')
HARTREE2EV = 27.211386245988


def setup():
    mol = gto.M(atom=ATOM, basis=BASIS, charge=0, spin=0, verbose=0)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-11
    mf.kernel()
    if not mf.converged:
        raise RuntimeError('RHF not converged; refusing to continue')
    return mol, mf


def fci_matrix(mol, mf):
    """Total-energy Hamiltonian in the determinant basis, via pyscf's FCI kernels."""
    norb = mf.mo_coeff.shape[1]
    nelec = mol.nelectron
    na, nb = nelec // 2 + nelec % 2, nelec // 2
    h1e = mf.mo_coeff.T @ mf.get_hcore() @ mf.mo_coeff
    eri = ao2mo.restore(1, ao2mo.kernel(mol, mf.mo_coeff), norb)
    h2e = fci.direct_spin1.absorb_h1e(h1e, eri, norb, (na, nb), 0.5)
    nda, ndb = comb(norb, na), comb(norb, nb)
    ndet = nda * ndb

    H = np.empty((ndet, ndet))
    for j in range(ndet):
        unit = np.zeros(ndet)
        unit[j] = 1.0
        H[:, j] = np.asarray(
            fci.direct_spin1.contract_2e(h2e, unit.reshape(nda, ndb),
                                         norb, (na, nb))).ravel()
    H += mol.energy_nuc() * np.eye(ndet)
    return H, (norb, na, nb, nda, ndb)


def main():
    mol, mf = setup()
    H, shape = fci_matrix(mol, mf)
    norb, na, nb, nda, ndb = shape
    ndet = H.shape[0]

    e_ref, _ = fci.FCI(mol, mf.mo_coeff).kernel()
    evals = np.linalg.eigvalsh(H)
    if abs(evals[0] - e_ref) > 1e-7:
        raise SystemExit(f'FCI matrix disagrees with fci.FCI: '
                         f'{evals[0]:.8f} vs {e_ref:.8f} Ha')
    print(f'water/{BASIS}: norb={norb} nelec={mol.nelectron} determinants={ndet}')
    print(f'FCI matrix validated against fci.FCI: {evals[0]:.8f} Ha '
          f'(diff {evals[0]-e_ref:+.1e})')

    D = np.diag(H).copy()
    O = H - np.diag(D)
    Ds = np.sort(D)
    gap = float(Ds[1] - Ds[0])
    corr = float(D.min() - evals[0])
    print(f'\ndiagonal: min={D.min():.6f} max={D.max():.6f} Ha')
    print(f'diagonal gap to the next determinant: {gap:.6f} Ha'
          f' = {gap*HARTREE2EV:.2f} eV')
    print(f'correlation energy (d - E_FCI)      : {corr:.6f} Ha'
          f' = {corr*HARTREE2EV:.2f} eV')

    tau_true = float(np.abs(np.linalg.eigvalsh(O)).max())
    print('\noff-diagonal part:')
    print(f'  ||O||_2 (true spectral norm)   = {tau_true:.6f} Ha'
          f' = {tau_true*HARTREE2EV:.2f} eV')
    print(f'  Frobenius                      = {np.linalg.norm(O, "fro"):.6f} Ha')
    print(f'  max row 2-norm                 = {np.linalg.norm(O, axis=1).max():.6f} Ha')
    print(f'  Gershgorin max row abs sum     = {np.abs(O).sum(axis=1).max():.6f} Ha')
    print(f'  max row nonzeros               = '
          f'{int((np.abs(O) > 1e-12).sum(axis=1).max())} of {ndet}')
    verdict = ('< 1: certificate decisive' if tau_true < gap
               else '> 1: sound but vacuous')
    print(f'\n  tau / diagonal gap             = {tau_true/gap:.2f}   ({verdict})')

    order = np.argsort(D)
    rows = []
    print('\ntrend as the complement shrinks (active = the k lowest-diagonal dets):')
    for k in (1, 5, 21, 50, 100, 200, 350):
        if k >= ndet:
            continue
        rest = order[k:]
        Oc = O[np.ix_(rest, rest)]
        tauc = float(np.abs(np.linalg.eigvalsh(Oc)).max())
        dc = float(D[rest].min())
        rows.append({'k_active': k, 'n_complement': int(rest.size),
                     'd': dc, 'tau': tauc, 'certified': dc - tauc})
        print(f'  k={k:4d} n_C={rest.size:4d}  d={dc:.4f}  tau={tauc:.4f}'
              f'  d-tau={dc-tauc:.4f} Ha')

    OUT.write_text(json.dumps({
        'system': {'atom': ATOM, 'basis': BASIS},
        'n_determinants': ndet,
        'rhf_energy': float(mf.e_tot),
        'fci_energy': float(evals[0]),
        'fci_reference': float(e_ref),
        'diagonal_min': float(D.min()),
        'diagonal_gap_to_next': gap,
        'correlation_energy': corr,
        'tau_spectral': tau_true,
        'tau_frobenius': float(np.linalg.norm(O, 'fro')),
        'tau_max_row_2norm': float(np.linalg.norm(O, axis=1).max()),
        'tau_gershgorin': float(np.abs(O).sum(axis=1).max()),
        'max_row_nnz': int((np.abs(O) > 1e-12).sum(axis=1).max()),
        'ratio_tau_over_gap': tau_true / gap,
        'sector_trend': rows,
    }, indent=2))
    print(f'\nwritten {OUT}')


if __name__ == '__main__':
    sys.exit(main())
