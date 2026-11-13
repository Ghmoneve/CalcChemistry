"""The Schur route: a certified lower bound that does not use the block constant tau.

WHY THIS ROUTE
The block certificate needs |Q_O(x)| <= tau*||x||^2 and tau >= ||O||_2 = 2.026 Ha
for water/STO-3G, while the correlation energy it must resolve is 0.049 Ha; the
bound is sound but vacuous. A subagent proved (and this project formalised the
steepest part of) that no diagonal similarity can lower that constant, because a
similarity preserves the spectrum and a symmetric matrix has 2-norm equal to its
spectral radius. So the block route is closed, not merely difficult.

THIS ROUTE
Split H around the reference determinant,

    H = [[ d,   b^T ],
         [ b,   H_rest ]],

with d = <ref|H|ref> = E_RHF, b the coupling column and H_rest the remaining
440x440 block. Then

    det(H - lam I) = (d - lam) det(Hr - lam I) - b^T adj(Hr - lam I) b,

so lam is an eigenvalue of H iff g(lam) = lam where

    g(lam) = d - b^T (Hr - lam I)^{-1} b,

provided lam < lambda_min(Hr) so the inverse exists.

WHAT IS AND IS NOT ESTABLISHED HERE
Numerically, the fixed point of g reproduces E_FCI to machine precision and
g is decreasing on the tested interval. Relaxing the assumed lower bound on
lambda_min(Hr) makes g(lam) LARGER, which is the opposite of the block route's
behaviour and is why this route might survive a crude spectral bound.

But the DIRECTION of the certified inequality is not yet settled: the fixed point
equation lam = g(lam) with g decreasing admits at most one solution, and it does
not by itself say whether g(mu) is a lower or an upper bound on lambda_1 for a
crude mu. This script therefore reports evidence toward the inequality and does NOT
assert the certificate. Any statement extending beyond the numbers below needs a
proof of the direction.

COST NOTE
Hr is 440x440 and exact diagonalisation of it is ~10^8 flop. The route is only
interesting if the bound on lambda_min(Hr) can be had more cheaply than that; the
script measures the spectrum so the comparison is possible, not as a proposed
algorithm.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from tau_measure import fci_matrix, setup

KC = 627.5094740631
OUT = Path('experiments/complement/schur_route.json')


def main():
    mol, mf = setup()
    H, _ = fci_matrix(mol, mf)
    n = H.shape[0]
    D = np.diag(H)
    iref = int(np.argmin(D))
    rest = [i for i in range(n) if i != iref]
    Hr = H[np.ix_(rest, rest)]
    b = H[rest, iref]

    w = np.linalg.eigvalsh(H)
    e0 = float(w[0])
    dr = np.linalg.eigvalsh(Hr)
    lam_min_rest = float(dr[0])
    d = float(D.min())

    res = {
        'n': n, 'd': d, 'E_FCI': e0, 'lambda_min_Hr': lam_min_rest,
        'spectrum_Hr_lowest': [float(x) for x in dr[:5]],
        'beta_row': float(np.linalg.norm(b)),
        'gap_diag': float(np.sort(D)[1] - np.sort(D)[0]),
        'correlation_energy': d - e0,
    }

    def g(lam):
        M = Hr - lam * np.eye(len(rest))
        return d - float(b @ np.linalg.solve(M, b))

    def gap_minus(lam):
        return g(lam) - lam

    # fixed point by bisection on g(lam) - lam
    lo, hi = lam_min_rest - 5.0, lam_min_rest - 1e-8
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if gap_minus(mid) > 0:
            lo = mid
        else:
            hi = mid
    fixed = 0.5 * (lo + hi)
    res['schur_fixed_point'] = fixed
    res['fixed_point_error_kcal'] = (fixed - e0) * KC

    # derivative sign over a range: is g monotone decreasing?
    eps = 1e-6
    derivs = {}
    for lam in (lam_min_rest - 2.0, lam_min_rest - 0.5, lam_min_rest - 0.1):
        dg = (g(lam + eps) - g(lam - eps)) / (2 * eps)
        derivs[f'{lam:.4f}'] = float(dg)
    res['dg_dlam'] = derivs
    res['g_decreasing'] = all(v < 0 for v in derivs.values())

    # How crude may the spectral bound be before the route degrades?
    rows = []
    print(f'n={n}  d={d:.6f}  E_FCI={e0:.6f}  corr={d-e0:.6f} Ha')
    print(f'lambda_min(Hr) = {lam_min_rest:.6f}  (= lambda_2 of H)')
    print(f'beta_row = ||b|| = {np.linalg.norm(b):.6f} Ha'
          f'  (block ||O||_2 = 2.026458, not used here)')
    print(f'\nfixed point of g = {fixed:.8f}  (E_FCI {e0:.8f},'
          f' diff {(fixed-e0)*KC:+.2e} kcal/mol)')
    print(f"g'(lam) at three probes: {derivs}")
    print(f'\n{"assumed lower bd on lambda_min(Hr)":38s}{"g(mu)":>14}'
          f'{"g(mu)-E_FCI (kcal)":>20}')
    for shift in (0.0, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0):
        mu = lam_min_rest - shift
        val = g(mu)
        rows.append({'shift': shift, 'mu': mu, 'g_mu': val,
                     'diff_kcal': (val - e0) * KC})
        print(f'  {"exact" if shift==0 else f"-{shift:.2f}":38s}{val:>14.6f}'
              f'{(val-e0)*KC:>20.2f}')
    res['spectral_bound_sensitivity'] = rows

    below = sum(1 for r in rows if r['g_mu'] < e0)
    print(f'\nof {len(rows)} probes, {below} gave g(mu) BELOW E_FCI')
    print('NOTE: the direction of the certified inequality is NOT established by')
    print('this script; see the module docstring. The numbers are evidence, not a proof.')

    OUT.write_text(json.dumps(res, indent=2))
    print(f'\nwritten {OUT}')


if __name__ == '__main__':
    main()
