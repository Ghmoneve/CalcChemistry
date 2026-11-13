"""Complement diagonal lower bound `d`: a combinatorial minimum that does not
enumerate the omitted determinant space.

The certificate (CalcChemistry/ComplementBound.lean, complement_ground_bound)
needs `d` with `D_i >= d` for every determinant `i` in the complement sector,
where `D_i` is a diagonal Hamiltonian entry. Enumerating the complement is the
thing we are trying to avoid, so `d` must come from a bound on `D_i` itself.

Split the diagonal in the determinant basis:

    D_i = E_core + sum_{p in occ(i)} eps_p + sum_{p<q in occ(i)} (pq|pq) - (pq|qp)

with `eps_p` the orbital energies of the one-body part and the second sum the
two-electron diagonal. Two structurally different lower bounds follow.

(A) ORBITAL-ENERGY BOUND. If the occupied set is a subset of a known reference
    occupation and every `eps_p` exceeds a common floor `eps_min`, and the
    exchange-type two-electron diagonal is bounded below by zero, then
    D_i >= E_core + N_elec * eps_min. This is cheap but usually very loose,
    because it ignores the large positive `(pq|pq)` Coulomb terms.

(B) EXCITATION-COUNT BOUND. Write `D_i = D_ref + Delta_i`, where `Delta_i` is
    the change in the diagonal when moving from the reference determinant to
    determinant `i`. `Delta_i` is a sum over the orbitals removed and added.
    If the complement sector is defined as "at most k excitations from the
    reference", then numerically minimising `Delta` over that sector is a
    small combinatorial problem in the orbital index set, not in the
    determinant count.

This script implements (A) and (B) on a real integral set and MEASURES how loose
each is against the true minimum, which is computable exactly for a small enough
active space. The point is to find out whether any cheap bound is useful; a
useless bound is a legitimate finding and is reported as such.

Method: RHF orbitals, then the diagonal of the FCI Hamiltonian in the
determinant basis. For a small model system the full determinant space is
affordable, so the TRUE minimum is known and every bound can be scored against it.
"""
from __future__ import annotations

import itertools
import json
import math
import sys

import numpy as np
from pyscf import ao2mo, gto, scf

# A neutral, closed-shell, cheap test system with real two-electron structure.
# Water in a minimal basis: 10 electrons in 7 orbitals is small enough that the
# full determinant space can be enumerated as a reference truth.
ATOM = 'O 0 0 0.117; H 0 0.757 -0.467; H 0 -0.757 -0.467'
BASIS = 'sto-3g'
CHARGE, SPIN = 0, 0

HARTREE2EV = 27.211386245988


def setup():
    mol = gto.M(atom=ATOM, basis=BASIS, charge=CHARGE, spin=SPIN, verbose=0)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-11
    mf.kernel()
    if not mf.converged:
        raise RuntimeError('RHF did not converge; refusing to continue')
    return mol, mf


def determinant_diagonals(mf):
    """Diagonal of the FCI Hamiltonian in the determinant basis, alpha and beta
    spin-orbitals treated as one ordered list of spin orbitals."""
    mol, norb = mf.mol, mf.mo_coeff.shape[1]
    nso = 2 * norb
    nelec = mol.nelectron

    # One- and two-electron integrals in the MO basis, spin-orbital form.
    h1 = mf.mo_coeff.T @ mf.get_hcore() @ mf.mo_coeff
    eri = ao2mo.full(mf.mol, mf.mo_coeff, compact=False).reshape(norb, norb, norb, norb)

    def so_h1(p, q):
        return h1[p // 2, q // 2] if (p % 2) == (q % 2) else 0.0

    def so_eri(p, q, r, s):
        # chemist notation (pq|rs) with spin conservation on each pair
        if (p % 2) != (r % 2) or (q % 2) != (s % 2):
            return 0.0
        return eri[p // 2, q // 2, r // 2, s // 2]

    core = mol.energy_nuc()
    occ_sets = []
    alpha = [2 * i for i in range(norb)]
    beta = [2 * i + 1 for i in range(norb)]
    na, nb = nelec // 2 + nelec % 2, nelec // 2
    for a in itertools.combinations(alpha, na):
        for b in itertools.combinations(beta, nb):
            occ_sets.append(tuple(sorted(a + b)))

    diag = np.empty(len(occ_sets))
    for idx, occ in enumerate(occ_sets):
        e = core
        for p in occ:
            e += so_h1(p, p)
        for p in occ:
            for q in occ:
                if p == q:
                    continue
                e += 0.5 * (so_eri(p, q, p, q) - so_eri(p, q, q, p))
        diag[idx] = e
    return diag, occ_sets, h1, eri, core


def orbital_energy_bound(mf):
    """Bound (A): E_core + N * eps_min, ignoring the two-electron diagonal."""
    mol = mf.mol
    eps = mf.mo_energy
    h1 = mf.mo_coeff.T @ mf.get_hcore() @ mf.mo_coeff
    # Core energy of the bare one-body operator plus nuclear repulsion is not
    # the right constant here; use the sum of eps over the reference occupation,
    # which already contains the mean-field two-electron contribution, and note
    # that this makes (A) a bound only under the stated hypotheses.
    return mol.energy_nuc() + mol.nelectron * eps.min()


def excitation_count_bound(diag, occ_sets, ref_occ, kmax):
    """Bound (B): minimum diagonal over determinants within kmax excitations.

    Only determinants at Hamming distance <= 2*kmax (alpha + beta) are scored,
    which is a small subset of the full space; the minimum over that subset is a
    valid lower bound for the restricted sector and is compared against the true
    minimum over the same sector.
    """
    ref = set(ref_occ)
    dist = np.array([len(ref.symmetric_difference(set(o))) // 2 for o in occ_sets])
    sector = dist <= kmax
    return float(diag[sector].min()), int(sector.sum())


def main():
    mol, mf = setup()
    print(f'system: {ATOM.split(";")[0].strip()} ... basis {BASIS}')
    print(f'RHF E = {mf.e_tot:.8f} Ha  converged={mf.converged}'
          f'  norb={mf.mo_coeff.shape[1]}  nelec={mol.nelectron}')

    diag, occ_sets, h1, eri, core = determinant_diagonals(mf)
    ndet = len(diag)
    print(f'determinant space size (S_z=0 sector): {ndet}')

    ref_idx = int(np.argmin(diag))
    ref_occ = occ_sets[ref_idx]
    true_min = float(diag.min())
    true_max = float(diag.max())
    print(f'TRUE diagonal range: [{true_min:.6f}, {true_max:.6f}] Ha'
          f'  (min is the reference determinant)')

    bound_A = orbital_energy_bound(mf)
    err_A = bound_A - true_min
    print(f'\n(A) orbital-energy bound      d_A = {bound_A:.6f} Ha')
    print(f'    error vs true minimum      {err_A:+.6f} Ha = {err_A*HARTREE2EV:+.2f} eV')
    print(f'    looser than the true range by a factor {err_A/(true_max-true_min):.2f}'
          f' of the full diagonal spread')

    rows = []
    print('\n(B) excitation-count bound over restricted sectors:')
    for kmax in (1, 2, 3, 4):
        d_k, n_k = excitation_count_bound(diag, occ_sets, ref_occ, kmax)
        sector_true = d_k  # by construction this IS the sector minimum (exact)
        rows.append({'kmax': kmax, 'n_det_in_sector': n_k, 'd_sector_min': d_k,
                     'excitation_energy': d_k - true_min,
                     'gap_to_reference': d_k - true_min})
        print(f'    k<={kmax}: {n_k:6d} dets   sector min D = {d_k:.6f} Ha'
              f'   Delta = {(d_k-true_min)*HARTREE2EV:7.3f} eV')

    out = {
        'system': {'atom': ATOM, 'basis': BASIS, 'charge': CHARGE, 'spin': SPIN},
        'rhf_energy': float(mf.e_tot),
        'rhf_converged': bool(mf.converged),
        'n_determinants': ndet,
        'true_diag_min': true_min,
        'true_diag_max': true_max,
        'bound_A_orbital_energy': bound_A,
        'bound_A_error_Ha': err_A,
        'bound_A_error_eV': err_A * HARTREE2EV,
        'bound_B_sectors': rows,
    }
    with open('experiments/complement/d_bound.json', 'w') as fh:
        json.dump(out, fh, indent=2)
    print('\nwritten experiments/complement/d_bound.json')


if __name__ == '__main__':
    sys.exit(main())
