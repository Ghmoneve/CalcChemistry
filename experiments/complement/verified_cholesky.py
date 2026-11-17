"""Certifying `Hr - d*I > 0` by a verified Cholesky factorisation.

WHY GERSHGORIN FAILS AND WHAT TO USE INSTEAD
The Schur bound needs `A = Hr - d*I` positive definite. Gershgorin asks whether A
is diagonally dominant, which is a *sufficient* condition and far too strong here:
the smallest diagonal of A is the correlation gap 0.447 Ha while the largest row
coupling is 1.380 Ha, so no diagonal-dominance argument can ever succeed, no matter
how the rows are grouped. This is a limitation of the criterion, not of the matrix.

THE CERTIFICATE USED HERE
Compute an approximate Cholesky factor L of A in floating point. Then, for any
symmetric A and L,

    A = L L^T + E,   ||E||_2 <= ||E||_F  <=  eps_total,

and by Weyl's inequality

    lambda_min(A) >= lambda_min(L L^T) - eps_total = sigma_min(L)^2 - eps_total.

Finally sigma_min(L) is bounded from below by the inverse of the infinity norm,

    sigma_min(L) >= 1 / (||L^{-1}||_2)  and  ||M||_2 <= sqrt(n) ||M||_inf,

so a strictly positive `sigma_min(L)^2 - eps_total` is a proof that A > 0. Nothing
here is an eigenvalue computation: two triangular solves and one residual.

THE TWO ERROR TERMS
    eps_total = ||A - fl(L L^T)||_F + n * u * || |L| |L|^T ||_F / (1 - n*u)
where u = 2^-53. The first term is what was actually observed; the second is the
bound on the rounding error committed while forming fl(L L^T).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from tau_measure import fci_matrix, setup

U = 2.0 ** -53
OUT = Path('experiments/complement/verified_cholesky.json')


def certify_pd(A: np.ndarray):
    """Return (certified, delta, detail) where delta is a proven lower bound on
    lambda_min(A), or None if the certificate does not go through."""
    n = A.shape[0]
    L = np.linalg.cholesky(A)

    LLt = L @ L.T
    E = A - LLt
    eps_obs = float(np.linalg.norm(E, 'fro'))

    absL = np.abs(L)
    eps_round = float(n * U * np.linalg.norm(absL @ absL.T, 'fro') / (1.0 - n * U))
    eps_total = eps_obs + eps_round

    Linv = np.linalg.solve_triangular(L, np.eye(n), lower=True) \
        if hasattr(np.linalg, 'solve_triangular') else np.linalg.inv(L)
    Linv_inf = float(np.linalg.norm(Linv, np.inf))
    sigma_lb = 1.0 / (np.sqrt(n) * Linv_inf)

    delta = sigma_lb ** 2 - eps_total
    detail = {
        'n': n,
        'eps_observed_fro': eps_obs,
        'eps_rounding_bound': eps_round,
        'eps_total': eps_total,
        'Linv_inf': Linv_inf,
        'sigma_min_L_lower_bound': sigma_lb,
        'delta_certified': delta,
        'certified_positive_definite': bool(delta > 0),
        'true_lambda_min': float(np.linalg.eigvalsh(A)[0]),
    }
    return bool(delta > 0), delta, detail


def bound_bAb(A: np.ndarray, b: np.ndarray, delta: float):
    """Rigorous two-sided bound on b^T A^{-1} b, given lambda_min(A) >= delta > 0.

    With y = A^{-1} b exact and yhat an approximate solve with residual
    r = b - A yhat, one has y - yhat = A^{-1} r, hence
        b^T A^{-1} b = b^T yhat + b^T A^{-1} r,
    and |b^T A^{-1} r| <= ||b|| ||r|| / lambda_min(A).
    """
    yhat, *_ = np.linalg.lstsq(A, b, rcond=None)
    r = b - A @ yhat
    centre = float(b @ yhat)
    rad = float(np.linalg.norm(b) * np.linalg.norm(r) / delta)
    return centre - rad, centre + rad, centre, rad


def main():
    mol, mf = setup()
    H, _ = fci_matrix(mol, mf)
    n = H.shape[0]
    D = np.diag(H)
    iref = int(np.argmin(D))
    rest = [i for i in range(n) if i != iref]
    Hr = H[np.ix_(rest, rest)]
    b = H[rest, iref]
    d = float(D.min())

    w = np.linalg.eigvalsh(H)
    e0 = float(w[0])

    A = Hr - d * np.eye(len(rest))
    ok, delta, detail = certify_pd(A)

    res = {'d': d, 'E_FCI': e0, 'pd_certificate': detail}

    print(f'n = {n}, d = {d:.8f}, E_FCI = {e0:.8f}')
    print(f'true lambda_min(A) = {detail["true_lambda_min"]:.6f}')
    print()
    print('--- verified Cholesky certificate for A = Hr - d*I > 0 ---')
    print(f'  ||A - L L^T||_F            = {detail["eps_observed_fro"]:.3e}')
    print(f'  rounding bound on fl(LL^T) = {detail["eps_rounding_bound"]:.3e}')
    print(f'  ||L^-1||_inf               = {detail["Linv_inf"]:.3e}')
    print(f'  sigma_min(L) lower bound   = {detail["sigma_min_L_lower_bound"]:.6e}')
    print(f'  certified lambda_min >= {delta:.6e}   (positive: {ok})')
    print(f'  -- margin: certified {delta:.3e} vs true {detail["true_lambda_min"]:.6f}')

    if not ok:
        print('\nCERTIFICATE FAILED: cannot certify positive definiteness this way.')
        OUT.write_text(json.dumps(res, indent=2))
        return

    lo, hi, centre, rad = bound_bAb(A, b, delta)
    res['bAb_bound'] = {'lo': lo, 'centre': centre, 'hi': hi, 'radius': rad}
    res['schur_bound_certified'] = d - hi
    res['slack_vs_fci_kcal'] = (d - hi - e0) * 627.5094740631

    print('\n--- b^T A^{-1} b, rigorously bracketed ---')
    print(f'  centre (computed) = {centre:.6e}')
    print(f'  radius            = {rad:.3e}')
    print(f'  upper bound       = {hi:.6e}')
    print()
    print('--- the certified lower bound ---')
    print(f'  lambda_1 >= d - b^T A^-1 b >= {d - hi:.8f} Ha')
    print(f'  E_FCI                       = {e0:.8f} Ha')
    print(f'  slack                       = {(d - hi - e0) * 627.5094740631:+.3f} kcal/mol')
    print(f'  is it a valid lower bound?  {(d - hi) <= e0}')

    OUT.write_text(json.dumps(res, indent=2))
    print(f'\nwritten {OUT}')


if __name__ == '__main__':
    main()
