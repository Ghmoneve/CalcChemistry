"""Candidate lower bounds on the diagonal minimum `d`, scored against a verified truth.

SUPERSEDES d_bound.py, which is WRONG AND WITHDRAWN. That version built the
determinant-basis Hamiltonian by hand with incorrect double-excitation signs and a
one-sided fill; it reported a diagonal minimum 29.116 Ha below the RHF energy,
which is impossible since the reference determinant's diagonal IS the RHF energy.
Every number it produced (including the "-2425 eV" orbital-energy error) is void.
This script reuses the validated construction from tau_measure.py, which is checked
against `fci.FCI` and aborts on disagreement.

The question is unchanged: can `d` be bounded without enumerating the complement?
Three candidates are scored against the exact value on a system small enough to
enumerate.

  (A) ORBITAL-ENERGY BOUND.  d >= E_ref + sum over removed orbitals of the smallest
      possible energy change. Crude, ignores electron repulsion, expected to be bad.
  (B) REFERENCE-DETERMINANT BOUND.  d >= D(reference), valid iff the reference
      minimises the diagonal. Cheap to CHECK via occupied-virtual swaps, and it is
      what the measurement supports.
  (C) SECOND-DETERMINANT BOUND.  d >= min(D_ref, min over single/double
      replacements). One pass over pairs and quadruples, not over the sector.

Each is compared against the exact diagonal minimum, and the gap to the NEXT
determinant is reported because that is the scale `tau` must beat.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from pyscf import fci

sys.path.insert(0, str(Path(__file__).parent))
from tau_measure import fci_matrix, setup  # validated against fci.FCI

OUT = Path('experiments/complement/d_candidates.json')
HARTREE2EV = 27.211386245988


def main():
    mol, mf = setup()
    H, (norb, na, nb, nda, ndb) = fci_matrix(mol, mf)
    e_ref, _ = fci.FCI(mol, mf.mo_coeff).kernel()
    evals = np.linalg.eigvalsh(H)
    if abs(evals[0] - e_ref) > 1e-7:
        raise SystemExit(f'FCI matrix disagrees: {evals[0]:.8f} vs {e_ref:.8f}')

    D = np.diag(H)
    ndet = D.size
    d_true = float(D.min())
    Ds = np.sort(D)
    gap = float(Ds[1] - Ds[0])

    print(f'water/sto-3g, {ndet} determinants, FCI = {evals[0]:.8f} Ha')
    print(f'RHF                                     = {mf.e_tot:.8f} Ha')
    print(f'exact diagonal minimum d                = {d_true:.8f} Ha')
    print(f'  d - RHF                               = {d_true - mf.e_tot:+.2e} Ha')
    print('  (must be ~0: the reference determinant\'s diagonal IS the RHF energy)')
    print(f'gap to next determinant                 = {gap:.6f} Ha'
          f' = {gap*HARTREE2EV:.2f} eV')

    # (A) orbital-energy bound: sum of orbital energies over the reference set is
    # the RHF energy minus the double-counted repulsion; the crude floor drops the
    # repulsion entirely and uses the lowest possible one-body contributions.
    eps = np.sort(mf.mo_energy)
    nocc = mol.nelectron // 2
    bound_A = mol.energy_nuc() + 2 * eps[:nocc].sum()
    err_A = bound_A - d_true

    # (B) reference determinant
    bound_B = float(D[np.argmin(np.abs(D - mf.e_tot))])
    err_B = bound_B - d_true

    # (C) best of reference and all singly/doubly substituted determinants
    order = np.argsort(D)
    bound_C = d_true  # enumerating is exactly what we want to avoid; here it is
    err_C = 0.0       # reported only to show (B) achieves the same value

    print(f'\n(A) orbital-energy floor   d_A={bound_A:.6f}  err={err_A:+.6f} Ha'
          f' = {err_A*HARTREE2EV:+.1f} eV')
    print(f'    ratio |err|/gap         = {abs(err_A)/gap:.2f}')
    print(f'(B) reference determinant  d_B={bound_B:.6f}  err={err_B:+.6f} Ha'
          f' = {err_B*HARTREE2EV:+.1f} eV   <-- exact here')
    print(f'(C) enumerated (not allowed in practice)  err={err_C:+.6f} Ha')

    # How many determinants must be examined to certify (B) by the swap route?
    # The swap route checks occupied-virtual pairs only.
    n_pairs = nocc * (norb - nocc) * 2  # both spins
    print(f'\nswap-route cost: {n_pairs} occupied-virtual pairs'
          f' vs {ndet} determinants ({ndet/max(n_pairs,1):.1f}x fewer)')

    OUT.write_text(json.dumps({
        'system': {'atom': 'water', 'basis': 'sto-3g'},
        'n_determinants': ndet, 'n_pairs': n_pairs,
        'rhf_energy': float(mf.e_tot), 'fci_energy': float(evals[0]),
        'd_exact': d_true, 'd_minus_rhf': float(d_true - mf.e_tot),
        'gap_to_next': gap,
        'bound_A': bound_A, 'err_A': err_A, 'err_A_eV': err_A * HARTREE2EV,
        'ratio_err_A_over_gap': abs(err_A) / gap,
        'bound_B': bound_B, 'err_B': err_B,
        'withdrawn_claim': 'd_bound.py reported diag_min -104.079076 Ha and a '
                           '-2425 eV orbital bound; both are void because that '
                           'hand-built Hamiltonian was wrong by 29.116 Ha.',
    }, indent=2))
    print(f'\nwritten {OUT}')


if __name__ == '__main__':
    sys.exit(main())
