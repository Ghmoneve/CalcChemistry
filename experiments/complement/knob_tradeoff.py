"""How does the cost of certifying the Schur hypothesis depend on the two knobs?

THE KNOBS
1. Shift `t`: the bound  lambda_1 >= (d - t) - b^T (Hr - (d-t)I)^{-1} b  is valid for
   any shift with `Hr - (d-t)I > 0`. Larger `t` makes the shift easier to certify
   positive definite but subtracts `t` from the bound.
2. Reference space size `k`: split off the k lowest determinants as a block. A larger
   reference space captures more correlation and pushes the remaining spectrum up,
   which should make the positive-definiteness condition easier.

WHY THIS MEASUREMENT DECIDES THE NEXT STEP
Gershgorin certifies `A + tI` only for `t > 1.0207`, which costs 640 kcal/mol: the
shift route with a diagonal criterion is dead. What is not known is whether a larger
reference space moves the diagonal-dominance threshold down far enough to matter. If
it does, a cheap certificate exists; if it does not, only factorisation works and the
O(n^3) cost is intrinsic to certification rather than an artefact of our choice.

WHAT IS AND IS NOT CLAIMED
This script measures margins. It does not claim any of the resulting bounds are
useful --- that is what it is measuring.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from tau_measure import fci_matrix, setup

KC = 627.5094740631
OUT = Path('experiments/complement/knob_tradeoff.json')


def gershgorin_shift_threshold(M: np.ndarray) -> float:
    """Smallest t >= 0 making M + t*I diagonally dominant (Gershgorin-positive)."""
    dg = np.diag(M)
    rs = np.abs(M).sum(1) - np.abs(dg)
    return float(np.max(rs - dg))


def main():
    mol, mf = setup()
    H, _ = fci_matrix(mol, mf)
    n = H.shape[0]
    D = np.diag(H)
    order = np.argsort(D)
    e0 = float(np.linalg.eigvalsh(H)[0])

    results = {'n': n, 'E_FCI': e0, 'reference_space_scan': [], 'shift_scan': []}

    print(f'n = {n}, E_FCI = {e0:.8f}\n')
    print('=== knob 2: reference space size k (t = 0) ===')
    print(f'{"k":>4}{"lambda_min(A)":>16}{"Gersh t-thresh":>18}'
          f'{"bound at t=0":>18}{"slack kcal":>14}')
    for k in (1, 2, 4, 8, 16, 32, 64):
        ref = order[:k]
        rest = order[k:]
        if len(rest) < 10:
            continue
        Hr = H[np.ix_(rest, rest)]
        B = H[np.ix_(rest, ref)]           # (n-k) x k coupling
        dk = np.linalg.eigvalsh(H[np.ix_(ref, ref)])
        lam = float(dk[0])

        # A = Hr - lam*I ; need A > 0 for the Schur step at this reference energy
        A = Hr - lam * np.eye(len(rest))
        lmin = float(np.linalg.eigvalsh(A)[0])
        thr = gershgorin_shift_threshold(A)

        # bound: lambda_min(H) >= lam - tr(B^T A^{-1} B)   (2x2 Schur: 1/(lam) form)
        # Use the exact block Schur value when A > 0.
        if lmin > 0:
            try:
                M = np.linalg.solve(A, B)
                corr = float(np.trace(B.T @ M))
                bound = lam - corr
            except np.linalg.LinAlgError:
                bound = float('nan')
        else:
            bound = float('nan')

        slack = (bound - e0) * KC if np.isfinite(bound) else float('nan')
        results['reference_space_scan'].append({
            'k': k, 'ref_energy': lam, 'lambda_min_A': lmin,
            'gershgorin_shift_threshold': thr, 'bound': bound, 'slack_kcal': slack,
        })
        print(f'{k:>4}{lmin:>16.6f}{thr:>18.6f}{bound:>18.6f}{slack:>14.2f}')

    # knob 1: shift at k = 1
    print('\n=== knob 1: shift t at k = 1 (single reference) ===')
    iref = int(order[0])
    rest = np.array([i for i in range(n) if i != iref])
    Hr = H[np.ix_(rest, rest)]
    b = H[rest, iref]
    d = float(D[iref])
    A = Hr - d * np.eye(len(rest))
    print(f'{"t":>10}{"lambda_min(A+tI)":>20}{"bound":>18}{"slack kcal":>14}')
    for t in (0.0, 0.05, 0.2, 0.5, 1.0207, 2.0):
        At = A + t * np.eye(len(rest))
        lmin = float(np.linalg.eigvalsh(At)[0])
        if lmin <= 0:
            print(f'{t:>10.4f}{lmin:>20.6f}{"not PD":>18}{"-":>14}')
            continue
        corr = float(b @ np.linalg.solve(At, b))
        bound = (d - t) - corr
        slack = (bound - e0) * KC
        results['shift_scan'].append({'t': t, 'lambda_min': lmin,
                                      'bound': bound, 'slack_kcal': slack})
        print(f'{t:>10.4f}{lmin:>20.6f}{bound:>18.6f}{slack:>14.2f}')

    OUT.write_text(json.dumps(results, indent=2))
    print(f'\nwritten {OUT}')


if __name__ == '__main__':
    main()
