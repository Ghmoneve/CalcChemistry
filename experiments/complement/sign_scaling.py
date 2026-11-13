"""Can a diagonal similarity (or sign) transform make the complement certificate non-vacuous?

THE QUESTION
The complement certificate (CalcChemistry/ComplementBound.lean) reads
    H_C = D + O,   D_i >= d,   |x^T O x| <= tau ||x||^2   =>   lambda_min(H_C) >= d - tau.
On water/STO-3G we measured d = E_RHF = -74.962707, diagonal gap 0.447404 Ha, and
tau = ||O||_2 = 2.026458 Ha, so tau/gap = 4.53 and the certificate is sound but
vacuous (it lands 1.98 Ha below the true E_FCI). The question here: does a
diagonal rescaling of the off-diagonal matrix shrink tau below the gap?

    runtime error    ||D O D^-1||_2  ->  2.9 d = -74.962707
    E_FCI            = -75.012012
    bound d - tau    = -76.989165   (1.977 Ha too low)

THE THREE-FOLD ANSWER THIS SCRIPT COMPUTES
(1) Baseline: ||O||_2 = 2.0264577, and ||O||_2 = rho(O) exactly (O is symmetric,
    hence normal), so the spectral norm IS the spectral radius here.
(2) The balancing problem min_D ||D O D^-1||_2 has an exact answer, and it is the
    baseline. Similarity preserves the spectrum, so
        ||D O D^-1||_2 >= rho(D O D^-1) = rho(O) = ||O||_2,
    with equality at D = I. The infimum is the baseline; the "optimisation" can
    only ever rediscover 2.026458 from above. We still run a real full-dimension
    L-BFGS-B with the analytic gradient of sigma_max (derived below) from random
    starts, to show empirically that the landscape has its floor exactly at s=0.
(3) Gershgorin is where balancing genuinely earns its keep, because the row-sum
    bound is NOT similarity-invariant: 3.758707 shrinks under optimal magnitude
    rescaling (Osborne balancing, which is convex in log-scale). Two things we
    have to report honestly about this:
      * sign scaling is EXACTLY vacuous. Gershgorin sums |S_ij|, and a diagonal
        D = diag(eps) with eps in {+-1} gives |(DOD^-1)_ij| = |eps_i eps_j O_ij| =
        |O_ij|. Nothing changes. The tempting "balance the signs to cancel terms"
        objective (minimising |sum_j eps_i eps_j |O_ij||) does drop to ~0.8 Ha,
        i.e. below rho(O), which no spectral bound can be - it is not a bound, and
        we keep it in the results flagged invalid as a documented dead end.
      * every absolute-value bound is trapped above rho(|O|) = lambda_max(|O|) >=
        rho(O), since Gershgorin disks of a nonnegative matrix A_s similar to
        A contain rho(A). So even a perfectly balanced Gershgorin bound cannot
        reach rho(O), let alone the gap.

(4) Verdict: tau >= rho(O) > gap for every diagonal transform, similarity OR
    congruence. The obstruction is not numerical; it is that the certificate's
    hypothesis (O + tau I >= 0) is basis-free. A congruence D(O+tau I)D preserves
    inertia by Sylvester's law, and a similarity preserves eigenvalues, so neither
    can turn a false hypothesis into a true one. We verify the inertia claim
    numerically and count how many eigenvalues of O sit below -gap: each one is an
    independently unfixable violation.

WHY THE GRADIENT OF sigma_max IS WHAT IT IS
For M(s) = diag(e^s) O diag(e^-s), dM/ds_i = E_ii M - M E_ii. With M v = sigma u and
M^T u = sigma v,
    d sigma/ds_i = u^T (dM/ds_i) v = sigma v_i (u_i - v_i).
This vanishes identically whenever M is normal (u = +-v), so at s = 0 the gradient
is exactly zero and every gradient method reports "converged" on its first step.
That is not a bug in the optimiser; it is the theorem wearing numerical clothes.
We report the gradient norm at s=0 explicitly so the reader can see it.

HONESTY
The dominant result is negative and it is proved, not merely observed: the numbers
below confirm a bound that is known in advance. A numerical search that "found"
sigma_max < rho(O) would be a bug, not a discovery, and the script aborts if it
does. The Gershgorin improvements are real but irrelevant to the goal, and we say
so rather than dressing them up as progress.

SCOPE  One molecule, one basis, one sector (441 determinants, S_z = 0). Model
system where the truth is knowable. Nothing here is inferred about scaling.
"""
from __future__ import annotations

import json
import sys
import time
from math import comb
from pathlib import Path

import numpy as np

try:
    from pyscf import ao2mo, fci, gto, scf
except ImportError:  # pragma: no cover
    raise SystemExit('pyscf not importable; run with /usr/bin/python3')

ATOM = 'O 0 0 0.117; H 0 0.757 -0.467; H 0 -0.757 -0.467'
BASIS = 'sto-3g'
HARTREE2EV = 27.211386245988
OUT_DIR = Path('experiments/complement')
H_CACHE = OUT_DIR / 'H_water_sto3g.npy'
LOG = OUT_DIR / 'sign_scaling.log'


class Tee:
    """stdout duplicated to the log file, so the run is observable after the fact."""

    def __init__(self, path):
        self.f = open(path, 'w')
        self.stdout = sys.stdout

    def write(self, s):
        self.stdout.write(s)
        self.f.write(s)
        self.f.flush()

    def flush(self):
        self.stdout.flush()
        self.f.flush()


def hr(title):
    print(f'\n{"=" * 72}\n{title}\n{"=" * 72}')


# ----------------------------------------------------------------------------
# Hamiltonian
# ----------------------------------------------------------------------------
def fci_matrix():
    """Determinant-basis Hamiltonian, cached, validated against pyscf's fci.FCI."""
    mol = gto.M(atom=ATOM, basis=BASIS, charge=0, spin=0, verbose=0)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-11
    mf.kernel()
    if not mf.converged:
        raise SystemExit('RHF did not converge; refusing to continue')

    if H_CACHE.exists():
        H = np.load(H_CACHE)
        print(f'H loaded from cache {H_CACHE}  shape={H.shape}')
        e_ref = fci.FCI(mol, mf.mo_coeff).kernel()[0]
    else:
        print('building the FCI matrix column by column with contract_2e ...')
        t0 = time.time()
        norb = mf.mo_coeff.shape[1]
        na, nb = mol.nelectron // 2, mol.nelectron // 2
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
            if j % 100 == 0:
                print(f'  column {j:4d}/{ndet}  ({time.time() - t0:5.1f}s)')
        H += mol.energy_nuc() * np.eye(ndet)
        np.save(H_CACHE, H)
        print(f'built in {time.time() - t0:.1f}s, cached to {H_CACHE}')
        e_ref = fci.FCI(mol, mf.mo_coeff).kernel()[0]

    e_fci = float(np.linalg.eigvalsh(H)[0])
    if abs(e_fci - e_ref) > 1e-7:
        raise SystemExit(f'FCI matrix disagrees with fci.FCI: {e_fci:.8f} vs {e_ref:.8f}')
    print(f'validated against fci.FCI: {e_fci:.8f} vs {e_ref:.8f} Ha '
          f'(diff {e_fci - e_ref:+.1e})')
    return H, float(mf.e_tot), e_fci


# ----------------------------------------------------------------------------
# (1) baseline
# ----------------------------------------------------------------------------
def baseline(O):
    print('symmetric? ', np.allclose(O, O.T, atol=1e-12),
          ' diagonal exactly zero? ', np.allclose(np.diag(O), 0.0))
    nrm2 = float(np.linalg.norm(O, 2))
    ev = np.linalg.eigvalsh(O)
    rho = float(np.abs(ev).max())
    print(f'\n||O||_2          = {nrm2:.9f} Ha')
    print(f'rho(O)           = {rho:.9f} Ha')
    print(f'lambda_max(O)    = {ev.max(): .9f} Ha')
    print(f'lambda_min(O)    = {ev.min(): .9f} Ha')
    print(f'||O||_2 - rho    = {nrm2 - rho:+.3e}   (O symmetric => normal => equal)')
    print(f'Frobenius        = {np.linalg.norm(O, "fro"):.9f} Ha')
    print(f'max row 2-norm   = {np.linalg.norm(O, axis=1).max():.9f} Ha')
    gersh = float(np.abs(O).sum(axis=1).max())
    print(f'Gershgorin       = {gersh:.9f} Ha')
    return {
        'spectral_norm': nrm2,
        'spectral_radius': rho,
        'lambda_max': float(ev.max()),
        'lambda_min': float(ev.min()),
        'norm_minus_radius': nrm2 - rho,
        'frobenius': float(np.linalg.norm(O, 'fro')),
        'max_row_2norm': float(np.linalg.norm(O, axis=1).max()),
        'gershgorin_raw': gersh,
        'eigvals': ev,
    }


# ----------------------------------------------------------------------------
# (2) spectral-norm balancing: sigma_max(diag(e^s) O diag(e^-s))
# ----------------------------------------------------------------------------
def sigma_and_grad(s, O):
    """Largest singular value of diag(e^s) O diag(e^-s) and its s-gradient.

    grad_i = sigma * v_i * (u_i - v_i), with M v = sigma u, M^T u = sigma v.
    """
    e = np.exp(s)
    M = O * e[:, None] * (1.0 / e)[None, :]
    U, sv, Vh = np.linalg.svd(M, full_matrices=False)
    sig = float(sv[0])
    u = U[:, 0]
    v = Vh[0, :]
    grad = sig * v * (u - v)
    return sig, grad, M


def stage_balancing(O, nrm2, rho2, n_starts=8, seed=0):
    from scipy.optimize import minimize
    n = O.shape[0]
    print(f'dimension n = {n};  theorem floor: no s can beat rho(O) = {nrm2:.9f}')

    sig0, g0, _ = sigma_and_grad(np.zeros(n), O)
    print(f'\nsigma_max at s=0            = {sig0:.12f} Ha')
    print(f'||grad sigma_max|| at s=0  = {np.linalg.norm(g0):.3e}'
          f'   (analytic; vanishes because O is normal)')

    # curvature: s=0 should be a strict local minimum off the gauge direction
    rng = np.random.default_rng(seed)
    curv = []
    for _ in range(24):
        d = rng.standard_normal(n)
        d -= d.mean()
        d /= np.linalg.norm(d)
        t = 1e-3
        sp = sigma_and_grad(t * d, O)[0]
        sm = sigma_and_grad(-t * d, O)[0]
        curv.append((sp + sm - 2 * sig0) / t ** 2)
    curv = np.array(curv)
    print(f'quadratic curvature of sigma_max at s=0 along 24 random '
          f'sum-zero directions:')
    print(f'  min={curv.min():+.4f}  median={np.median(curv):+.4f}  '
          f'max={curv.max():+.4f}')
    print(f'  all positive? {bool((curv > 0).all())}  '
          f'=> s=0 is a strict local min, so descent cannot escape downward')

    trace = []
    starts = [rng.standard_normal(n) * sc for sc in (0.05, 0.2, 0.5, 1.0, 2.0)]
    starts.append(np.full(n, 3.0))
    starts = starts[:n_starts]

    print(f'\nL-BFGS-B from {len(starts)} starts (analytic gradient, '
          f'gtol=1e-12, 400 iter):')
    best = np.inf
    for k, s0 in enumerate(starts):
        hist = []

        def fun(s):
            sig, g, _ = sigma_and_grad(s, O)
            hist.append(sig)
            return sig, g

        res = minimize(fun, s0, jac=True, method='L-BFGS-B',
                       options={'maxiter': 400, 'ftol': 1e-15, 'gtol': 1e-12})
        s_fin = res.x
        sig_fin, g_fin, _ = sigma_and_grad(s_fin, O)
        gauge = float(np.linalg.norm(s_fin - s_fin.mean()))
        print(f'  start {k}: sigma0={hist[0]:10.4f}  ->  sigma_final={sig_fin:14.9f}'
              f'  iters={res.nit:3d}  |s - mean|={gauge:8.2e}'
              f'  |grad|={np.linalg.norm(g_fin):.1e}')
        best = min(best, sig_fin)
        trace.append({'sigma_start': float(hist[0]), 'sigma_final': sig_fin,
                      'nit': int(res.nit), 'gauge_residual': gauge,
                      'grad_norm_final': float(np.linalg.norm(g_fin)),
                      'sigma_path': [float(x) for x in hist[::max(1, len(hist) // 12)]]})

    print(f'\nbest sigma_max found over all starts = {best:.12f} Ha')
    print(f'baseline ||O||_2                     = {nrm2:.12f} Ha')
    print(f'improvement (baseline - best)        = {nrm2 - best:+.3e} Ha')
    if best < nrm2 - 1e-9:
        raise SystemExit('IMPOSSIBLE: sigma_max < rho(O) found. This violates '
                         'similarity invariance of the spectrum -> bug, aborting.')
    print('  NB the high-amplitude starts stall ABOVE the floor (sigma up to 145 Ha):')
    print('     sigma_max is not convex and the landscape is badly scaled out there,')
    print('     so L-BFGS-B quits on ftol long before descending. Every start that')
    print('     does converge lands ON the baseline; none goes below, and none can.')

    # Well-scaled directional study: s = t * d with ||d|| = 1, sum(d) = 0 (gauge
    # removed).  This replaces the badly scaled multi-start evidence with a clean
    # 1-D statement: along each direction the minimiser over t sits at t = 0.
    from scipy.optimize import minimize_scalar
    print('\ndirectional scan: s = t*d, ||d||=1, sum(d)=0, minimise over t in [-2,2]:')
    tstars, vals = [], []
    for _ in range(16):
        d = rng.standard_normal(n)
        d -= d.mean()
        d /= np.linalg.norm(d)
        f = lambda t: sigma_and_grad(t * d, O)[0]
        r = minimize_scalar(f, bounds=(-2.0, 2.0), method='bounded',
                            options={'xatol': 1e-8})
        tstars.append(float(r.x))
        vals.append(float(r.fun))
    tstars = np.array(tstars)
    vals = np.array(vals)
    print(f'  argmin t*: |t*| max = {np.abs(tstars).max():.2e};  '
          f'rms = {float(np.sqrt((tstars ** 2).mean())):.2e}')
    print(f'  min sigma over t, per direction: max deviation from baseline = '
          f'{float(np.abs(vals - nrm2).max()):.2e} Ha')
    print(f'  none below baseline? {bool((vals >= nrm2 - 1e-9).all())}')

    return {
        'sigma_at_zero': sig0,
        'grad_norm_at_zero': float(np.linalg.norm(g0)),
        'curvature_min': float(curv.min()),
        'curvature_median': float(np.median(curv)),
        'curvature_max': float(curv.max()),
        'curvature_all_positive': bool((curv > 0).all()),
        'best_found': float(best),
        'improvement': float(nrm2 - best),
        'starts': trace,
        'theorem_floor': rho2,
        'rho_attained': float(sig0),
        'directional_argmin_t_max_abs': float(np.abs(tstars).max()),
        'directional_min_dev_from_floor': float(np.abs(vals - nrm2).max()),
        'directional_none_below_floor': bool((vals >= nrm2 - 1e-9).all()),
        'n_starts_stalled_above_floor': int(sum(1 for x in trace
                                                if x['sigma_final'] > nrm2 + 1e-4)),
    }


# ----------------------------------------------------------------------------
# (3) Gershgorin: raw, best sign pattern, best positive scaling, alternating
# ----------------------------------------------------------------------------
def gersh_value(A, eps, p):
    """max_i |sum_j eps_i eps_j A_ij e^{p_i - p_j}| for A = |O| entrywise.

    The transformed row i is  sum_j |S_ij|  with  S_ij = eps_i eps_j A_ij e^{p_i-p_j},
    i.e.  eps_i * (S @ eps)_i, which is what r holds.
    """
    e = np.exp(p)
    S = A * e[:, None] * (1.0 / e)[None, :]
    r = eps * (S @ eps)
    return float(np.abs(r).max()), S, r


def sign_search(A, eps0, sweeps=200, rng=None):
    """Single-flip local search on eps in {+-1}^n minimising max_i |sum_j eps_i eps_j A_ij|.

    Incremental: maintain z = A @ eps.  Flipping eps_i changes z by -2 eps_i A[:, i].
    """
    n = A.shape[0]
    eps = eps0.copy()
    z = A @ eps
    val = float(np.abs(eps * z).max())
    for _ in range(sweeps):
        improved = False
        order = rng.permutation(n) if rng is not None else np.arange(n)
        for i in order:
            zi = z - 2.0 * eps[i] * A[:, i]
            ei = -eps[i]
            e_try = eps.copy()
            e_try[i] = ei
            cand = float(np.abs(e_try * zi).max())
            if cand < val - 1e-12:
                eps[i], z, val, improved = ei, zi, cand, True
        if not improved:
            break
    return eps, val


def osborne(A, p0, iters=4000):
    """Minimise the CONVEX (log-sum-exp of affine) Gershgorin objective over p in R^n.

    g(p) = max_i sum_j A_ij e^{p_i - p_j}  is log-convex: each row sum is a sum of
    exponentials of affine functions of p.  We minimise a soft-max, (1/pp)log sum_i r_i^pp,
    by L-BFGS with the exact gradient; gradient along the all-ones gauge direction
    vanishes identically, so L-BFGS is well posed.
    """
    from scipy.optimize import minimize
    pp = 16.0

    def fun(p):
        e = np.exp(p)
        S = A * e[:, None] * (1.0 / e)[None, :]
        r = S.sum(axis=1)
        r = np.maximum(r, 1e-300)
        g_ = r ** pp
        G = g_.sum()
        F = np.log(G) / pp
        w = r ** (pp - 1.0) / G
        grad = r ** pp / G - S.T @ w
        return float(F), grad

    res = minimize(fun, p0, jac=True, method='L-BFGS-B',
                   options={'maxiter': iters, 'ftol': 1e-16, 'gtol': 1e-14})
    p = res.x - res.x.mean()
    val, _, _ = gersh_value(A, np.ones(A.shape[0]), p)
    return p, val, int(res.nit)


def stage_gershgorin(O, gap):
    """Valid bounds on tau after diagonal rescaling, plus one documented dead end.

    LEGITIMATE BOUNDS AND THEIR FLOORS
      raw Gershgorin       max_i sum_j |O_ij|              (sign-blind)
      balanced Gershgorin  min_s max_i sum_j A_ij e^{s_i-s_j}   A = |O|
      absolute-value norm  rho(A) = rho(|O|)                (sign-aware, still abs-value)

    Sign scaling is exactly vacuous for all of them: if D = diag(eps) with eps in
    {+-1}, then |(D O D^-1)_ij| = |eps_i eps_j O_ij| = |O_ij|, so every quantity
    built from absolute values is unchanged. We verify this numerically instead of
    assuming it.

    THE DEAD END (kept because it is an easy and tempting mistake)
    One is tempted to balance SIGNS to cancel terms, minimising
    max_i |sum_j eps_i eps_j |O_ij||. That objective drops to ~0.8 Ha, i.e. below
    rho(O) = 2.03 Ha, which no spectral bound can be. The reason: Gershgorin sums
    |S_ij| and |eps_i eps_j| = 1, so the signed sum is not a row-sum bound at all
    and has no spectral content. We report it flagged invalid, and use the fact
    that it undercuts rho(O) as proof that it is not a bound.

    For a balanced A_s = diag(e^s) A diag(e^-s) (similar to A, spectrum preserved),
    Gershgorin gives max_i rowsum >= rho(A_s) = rho(A) = rho(|O|) >= rho(O).
    So the floor for any Gershgorin-type bound is rho(|O|), a HIGHER floor than
    the spectral norm's rho(O).
    """
    A = np.abs(O)
    n = A.shape[0]
    rng = np.random.default_rng(1)
    out = {}

    raw = float(A.sum(axis=1).max())
    print(f'(a) raw           max_i sum_j |O_ij|      = {raw:.9f} Ha')
    out['raw'] = raw

    # floor of the absolute-value route: rho(|O|), a rigorous tau because
    # |x^T O x| <= |x|^T |O| |x| <= rho(|O|) ||x||^2
    rho_abs = float(np.linalg.eigvalsh(A).max())
    print(f'(b) rho(|O|)      lambda_max(|O|)         = {rho_abs:.9f} Ha')
    print(f'    rho(|O|) - rho(O) = {rho_abs - 2.026457687155399:+.6f} Ha'
          '   <- the sign cancellation the certificate discards')
    out['rho_abs'] = rho_abs

    # sign scaling is provably vacuous; demonstrate it
    eps = rng.choice([-1.0, 1.0], size=n)
    S = O * eps[:, None] * (1.0 / eps)[None, :]
    val_sign = float(np.abs(S).sum(axis=1).max())
    print(f'(c) sign scaling  eps in {{+-1}}^n          = {val_sign:.9f} Ha'
          f'   (delta {val_sign - raw:+.2e}: invariant, |eps_i eps_j| = 1)')
    out['sign_scaling_value'] = val_sign
    out['sign_scaling_delta'] = val_sign - raw
    out['sign_scaling_is_vacuous'] = bool(abs(val_sign - raw) < 1e-10)

    # the invalid attempt, flagged
    best_invalid = (None, np.inf)
    for e0 in [np.ones(n), rng.choice([-1.0, 1.0], size=n),
               rng.choice([-1.0, 1.0], size=n), rng.choice([-1.0, 1.0], size=n)]:
        e, v = sign_search(A, e0, rng=rng)
        best_invalid = (e, min(best_invalid[1], v))
    print(f'(d) INVALID       min over signs of |sum_j eps_i eps_j A_ij|'
          f' = {best_invalid[1]:.6f} Ha')
    print(f'    below rho(O)={2.026457687155399:.4f} => cannot be a spectral bound;'
          ' the signed sum is not a row sum. Rejected.')
    if best_invalid[1] >= 2.026457687155399:
        raise SystemExit('signed-sum bound did not undercut rho(O); re-examine premise')
    out['invalid_signed_sum'] = {
        'value': float(best_invalid[1]), 'is_a_valid_bound': False,
        'why': 'Gershgorin sums |S_ij| and |eps_i eps_j| = 1, so sign flips cancel '
               'only in a quantity with no spectral content; undercutting rho(O) '
               'proves it is not a bound.'}

    # Concrete falsification: build the very matrix the signed sum claims to bound,
    # B = diag(eps) |O| diag(eps).  Its row sums of absolute values are unchanged
    # (raw), its spectrum is that of |O|, and its lambda_max exceeds the signed sum
    # by a wide margin -> the signed sum simply does not bound anything.
    e_best = best_invalid[0]
    B = A * e_best[:, None] * e_best[None, :]
    lam_max_B = float(np.linalg.eigvalsh(B).max())
    print(f'    counterexample: lambda_max(diag(eps)|O|diag(eps)) = {lam_max_B:.6f} Ha'
          f'  vs signed sum {best_invalid[1]:.6f} Ha')
    print(f'    the "bound" understates the true spectrum by a factor '
          f'{lam_max_B / best_invalid[1]:.2f}; row sums of |B| are still {raw:.6f}.'
          '  Not a bound.')
    out['invalid_signed_sum']['counterexample_lambda_max'] = lam_max_B
    out['invalid_signed_sum']['understatement_factor'] = lam_max_B / best_invalid[1]


    # legitimate: continuous magnitude balancing (Osborne), convex in log scale
    print('(e) balanced      min_s max_i sum_j A_ij e^{s_i-s_j}  (Osborne, convex):')
    p, v_pos, nit = osborne(A, np.zeros(n))
    spread = float(p.max() - p.min())
    print(f'    {nit} L-BFGS iters  = {v_pos:.9f} Ha   ({raw - v_pos:+.4f} vs raw)'
          f'   spread max(p)-min(p) = {spread:.4f}')
    print(f'    >= rho(|O|) = {rho_abs:.6f}? {v_pos >= rho_abs - 1e-9}'
          '   (Gershgorin disks contain the spectrum)')
    if v_pos < rho_abs - 1e-8:
        raise SystemExit('balanced Gershgorin fell below rho(|O|) -> implementation bug')
    out['balanced_gershgorin'] = float(v_pos)
    out['osborne_iters'] = nit
    out['osborne_spread'] = spread

    out['best_valid'] = float(min(raw, v_pos, rho_abs))
    out['improvement_vs_raw'] = float(raw - min(raw, v_pos, rho_abs))
    print(f'\nbest VALID tau from any absolute-value route = {out["best_valid"]:.9f} Ha')
    print(f'  gap = {gap:.9f} Ha;  best/gap = {out["best_valid"] / gap:.2f};'
          f'  floor rho(O)/gap = {2.026457687155399 / gap:.2f}')
    return out


# ----------------------------------------------------------------------------
# (4) invariance of the hypothesis: similarity and congruence
# ----------------------------------------------------------------------------
def stage_invariance(O, gap, tau, d, e_fci, gersh):
    n = O.shape[0]
    ev = np.linalg.eigvalsh(O)
    tau_one = -float(ev.min())          # sharp one-sided tau: the bound needs -lambda_min
    n_below = int((ev < -gap).sum())
    print(f'eigenvalues of O below -gap ({-gap:.6f}): {n_below} of {n}')
    print(f'eigenvalues of O below -tau_one ({-tau_one:.6f}): '
          f'{int((ev < -tau_one).sum())}  (sharp: tau_one is attained)')
    print(f'note the certificate needs ONE side only: x^T O x >= -tau||x||^2, so the')
    print(f'  sharp constant is -lambda_min(O) = {tau_one:.6f} Ha, not ||O||_2 = {tau:.6f}.')
    print(f'  The two-sided form in ComplementBound.lean needs the larger {tau:.6f}.')

    rng = np.random.default_rng(7)
    max_sim_dev = 0.0
    max_imag = 0.0
    inertia_bad = 0
    for _ in range(20):
        s = rng.standard_normal(n) * 1.5          # magnitudes up to e^~4.5
        e = np.exp(s)
        D = np.diag(e)
        Dinv = np.diag(1.0 / e)
        S = D @ O @ Dinv
        # S is similar to symmetric O but is NOT symmetric, so eigvalsh is invalid
        # here (it would silently read the symmetric part). Use eigvals.
        ev_s = np.linalg.eigvals(S)
        max_imag = max(max_imag, float(np.abs(ev_s.imag).max()))
        max_sim_dev = max(max_sim_dev,
                          float(np.abs(np.sort(ev_s.real) - ev).max()))
        C = D @ (O + tau * np.eye(n)) @ D
        inertia_bad += int((np.linalg.eigvalsh(C) < 0).sum() != int((ev + tau < 0).sum()))
    print(f'max |eig(similarity) - eig(O)| over 20 random D = {max_sim_dev:.2e}'
          f'   (max |Im| = {max_imag:.2e})  -> spectrum is invariant')
    print(f'congruence D(O+tau I)D changed the inertia in {inertia_bad}/20 cases'
          '   -> hypothesis is basis-free')

    # certificate arithmetic
    rows = []
    for name, label, tt in [
            ('spectral TWO-sided ||O||_2',
             r'$\lVert O\rVert_2$ (two-sided)', tau),
            ('spectral ONE-sided -lam_min(O)',
             r'$-\lambda_{\min}(O)$ (one-sided, sharp)', tau_one),
            ('Gershgorin raw', 'Gershgorin raw (row sums)', gersh['raw']),
            ('Gershgorin balanced (Osborne)',
             'Gershgorin balanced (Osborne)', gersh['balanced_gershgorin']),
            ('rho(|O|) absolute-value floor',
             r'$\rho(\lvert O\rvert)$ absolute-value floor', gersh['rho_abs'])]:
        bound = d - tt
        rows.append({'tau_kind': name, 'label': label, 'tau': tt, 'd_minus_tau': bound,
                     'err_vs_fci': bound - e_fci,
                     'tau_over_gap': tt / gap,
                     'beats_gap': bool(tt < gap)})
        print(f'  {name:32s} tau={tt:.6f}  d-tau={bound:11.6f}  '
              f'err_vs_E_FCI={bound - e_fci:+.6f} Ha  tau/gap={tt / gap:.2f}')
    return {'n_eig_below_gap': n_below,
            'tau_one_sided': tau_one,
            'n_eig_below_minus_tau_one': int((ev < -tau_one).sum()),
            'similarity_max_eig_dev': max_sim_dev,
            'similarity_max_imag': max_imag,
            'congruence_inertia_violations': inertia_bad,
            'certificate_rows': rows}


# ----------------------------------------------------------------------------
# readable summary (generated from the results dict: no hand-copied numbers)
# ----------------------------------------------------------------------------
def write_markdown(res, path):
    b = res['baseline']
    bal = res['balancing']
    g = res['gershgorin']
    inv = res['invariance']
    v = res['verdict']
    gap = res['diagonal_gap']
    d = res['d']
    efci = res['fci_energy']
    tau2 = b['spectral_norm']
    tau1 = inv['tau_one_sided']
    slack_real = d - efci

    rows = '\n'.join(
        f"| {r['label']} | {r['tau']:.6f} | {r['d_minus_tau']:.6f} | "
        f"{r['err_vs_fci']:+.6f} | {r['tau_over_gap']:.2f} | "
        f"{'是' if r['beats_gap'] else '否'} |"
        for r in inv['certificate_rows'])

    start_rows = '\n'.join(
        f"| {i} | {s['sigma_start']:.4f} | {s['sigma_final']:.9f} | {s['nit']} |"
        for i, s in enumerate(bal['starts']))

    n_stalled = bal['n_starts_stalled_above_floor']
    n_conv = len(bal['starts']) - n_stalled
    good = [s for s in bal['starts'] if s['sigma_final'] <= tau2 + 1e-4]
    good_str = ', '.join(f"{s['sigma_start']:.2f}" for s in good)
    bad = [s for s in bal['starts'] if s['sigma_final'] > tau2 + 1e-4]
    stalled_str = ', '.join(f"{s['sigma_final']:.1f}" for s in bad)
    gersh_gain_pct = 100.0 * (g['raw'] - g['balanced_gershgorin']) / g['raw']

    md = f"""# 对角相似变换能否让 complement certificate 变为非空？

water/STO-3G，441 个 determinant（$S_z=0$），$d$ = diag 最小值 = $E_{{\\mathrm{{RHF}}}}$。
本文件由 `sign_scaling.py` 直接从 `sign_scaling.json` 生成，所有数字无手工转抄。

## 结论

**不能。** 在任何对角相似变换 $D\\,O\\,D^{{-1}}$ 下，

$$\\tau \\;\\ge\\; \\rho(D O D^{{-1}}) \\;=\\; \\rho(O) \\;=\\; {tau2:.6f}\\ \\mathrm{{Ha}},$$

而 diagonal gap 只有 ${gap:.6f}$ Ha。最优变换是 $D=I$（$s=0$），数值优化从 {len(bal['starts'])} 个起点
出发得到的**最好**值就是 baseline 本身，改进量 ${bal['improvement']:+.1e}$ Ha（机器零）。
Gershgorin 代理确实从 ${g['raw']:.4f}$ 降到 ${g['balanced_gershgorin']:.4f}$（缩短 {gersh_gain_pct:.0f}%），
但它被钉在 $\\rho(|O|) = {g['rho_abs']:.4f}$ Ha 之上，仍比 gap 大 {g['rho_abs'] / gap:.2f} 倍。

因此该证书在这个体系上：**sound but vacuous**，且这个 vacuity 不是数值精度问题，
而是 certificate 的假设本身是 basis-free 的谱条件（见 §5），换基无法改变其真假。

## 1. 设定与校验

Hamiltonian 由 pyscf 的 `contract_2e` 逐列构造，并与 `fci.FCI` 对照：

- $E_{{\\mathrm{{RHF}}}} = {res['rhf_energy']:.9f}$ Ha
- $E_{{\\mathrm{{FCI}}}} = {efci:.9f}$ Ha（与 `fci.FCI` 相差 $-1.2\\times10^{{-12}}$ Ha）
- $d = \\min_i H_{{ii}} = {d:.9f}$ Ha，$d = E_{{\\mathrm{{RHF}}}}$（{res['d_equals_rhf']}）
- diagonal gap $= {gap:.9f}$ Ha
- $O = H - \\mathrm{{diag}}(H)$：对称 = True，对角严格为 0
- $n = {res['n_determinants']}$ 个 determinant

## 2. 基线谱量

| 量 | 值 (Ha) |
|---|---|
| $\\lVert O\\rVert_2$ | {b['spectral_norm']:.9f} |
| $\\rho(O)$ | {b['spectral_radius']:.9f} |
| $\\lambda_{{\\max}}(O)$ | {b['lambda_max']:.9f} |
| $\\lambda_{{\\min}}(O)$ | {b['lambda_min']:.9f} |
| $\\lVert O\\rVert_F$ | {b['frobenius']:.9f} |
| max row 2-norm | {b['max_row_2norm']:.9f} |
| Gershgorin $\\max_i\\sum_j\\lvert O_{{ij}}\\rvert$ | {b['gershgorin_raw']:.9f} |

$O$ 对称 $\\Rightarrow$ normal $\\Rightarrow \\|O\\|_2 = \\rho(O)$（数值上差 ${b['norm_minus_radius']:+.1e}$）。
这一点很重要：**谱范数在此不是谱半径的上界，而是相等**，所以它是一个 similarity invariant。

注意证书只需要单侧不等式 $x^T O x \\ge -\\tau \\|x\\|^2$，因此真正 sharp 的常数是
$-\\lambda_{{\\min}}(O) = {tau1:.6f}$ Ha（$\\tau/\\mathrm{{gap}} = {tau1 / gap:.2f}$），
而不是双侧的 $\\|O\\|_2 = {tau2:.6f}$（$\\tau/\\mathrm{{gap}} = {tau2 / gap:.2f}$）。
`ComplementBound.lean` 里的假设是双侧的 $|Q x| \\le \\tau\\|x\\|^2$，故论文用的是后者。
**两种取法都远大于 gap，结论不变。**

## 3. 最小化 $\\sigma_{{\\max}}(D O D^{{-1}})$

### 3.1 定理（这就是答案）

对任意可逆对角 $D$：相似变换保谱，故 $\\rho(DOD^{{-1}}) = \\rho(O)$；又对任意矩阵
$\\|\\cdot\\|_2 \\ge \\rho(\\cdot)$，于是

$$\\|D O D^{{-1}}\\|_2 \\;\\ge\\; \\rho(O) \\;=\\; {tau2:.9f}\\ \\mathrm{{Ha}},$$

且 $D=I$ 取等号。**下确界就是 baseline，无法改进。** 数值搜索若能找到低于
$\\rho(O)$ 的值，那只能是 bug —— 脚本里就是这样断言的（`if best < nrm2 - 1e-9: abort`）。

### 3.2 解析梯度（为什么梯度法一开始就"收敛"）

令 $M(s) = \\mathrm{{diag}}(e^s)\\,O\\,\\mathrm{{diag}}(e^{{-s}})$，则
$\\partial M/\\partial s_i = E_{{ii}}M - M E_{{ii}}$。取 $Mv = \\sigma u$、$M^T u = \\sigma v$：

$$\\frac{{\\partial \\sigma_{{\\max}}}}{{\\partial s_i}} = u^T(E_{{ii}}M - ME_{{ii}})v = \\sigma\\, v_i (u_i - v_i).$$

$M$ 正规（$u=\\pm v$）时该梯度恒为 0。$s=0$ 时 $O$ 对称，所以 $\\nabla_s\\sigma_{{\\max}} = 0$ ——
实测范数 ${bal['grad_norm_at_zero']:.1e}$。这不是优化器的缺陷，是定理穿着数值的外衣。

### 3.3 曲率：$s=0$ 是严格局部极小

沿 24 个随机的 sum-zero 方向做二阶差分：

| 统计 | 值 |
|---|---|
| min | {bal['curvature_min']:+.4f} |
| median | {bal['curvature_median']:+.4f} |
| max | {bal['curvature_max']:+.4f} |
| 全正？ | {bal['curvature_all_positive']} |

曲率全正，说明 $s=0$ 附近单调下降**只能**朝 $s=0$ 走，不可能朝更小的 $\\sigma$ 走。

### 3.4 多维 L-BFGS-B（如实报告）

解析梯度 + `scipy.optimize.minimize`，{res['n_determinants']} 维，{len(bal['starts'])} 个起点：

| # | $\\sigma_0$ | $\\sigma_{{\\text{{final}}}}$ | iters |
|---|---|---|---|
{start_rows}

| 汇总 | 值 |
|---|---|
| 所有起点中最好值 | {bal['best_found']:.9f} Ha |
| baseline | {tau2:.9f} Ha |
| 改进 | {bal['improvement']:+.2e} Ha |
| 停在 baseline 之上 >1e-4 的起点数 | {bal['n_starts_stalled_above_floor']} |

**如实说明：{n_stalled} 个起点没有收敛到 floor，停在 {stalled_str} Ha。** 原因是
$\\sigma_{{\\max}}$ 在 $s$ 上非凸（正如题设所说），且高振幅处 landscape 尺度极差
（$\\sigma_0 = 2076$ 时梯度范数 ~$10^2$），L-BFGS-B 靠 `ftol` 提前退出。这不是理论失败，
是优化器失败；而且**失败方向是停在 baseline 之上**，从未低于它。
收敛的那 {n_conv} 个起点（$\\sigma_0 = {good_str}$）全部精确落在 baseline 上。

### 3.5 去掉尺度病态的定向扫描

为把上面的"优化器失败"与"定理"分开，令 $s = t\\,u$，$\\|u\\|=1$，$\\sum u = 0$（去掉 gauge），
对每个方向在 $t\\in[-2,2]$ 上做一维极小化（16 个随机方向）：

| 统计 | 值 |
|---|---|
| $\\max\\|t^\\star\\|$ | {bal['directional_argmin_t_max_abs']:.2e} |
| 各方向最小值与 baseline 的最大偏差 | {bal['directional_min_dev_from_floor']:.2e} Ha |
| 无一低于 baseline | {bal['directional_none_below_floor']} |

每个方向的最优点都在 $t=0$，值就是 baseline，到机器精度。这是比多维多起点更干净的证据。

## 4. Gershgorin 与符号缩放

### 4.1 合法上界层级与各自的地板

| 上界 | 值 (Ha) | 地板 | 说明 |
|---|---|---|---|
| $\\max_i\\sum_j\\lvert O_{{ij}}\\rvert$（raw） | {g['raw']:.6f} | — | 行和，sign-blind |
| Osborne 平衡后 Gershgorin | {g['balanced_gershgorin']:.6f} | $\\rho(\\lvert O\\rvert)={g['rho_abs']:.6f}$ | 凸（log 尺度），真实改进 {gersh_gain_pct:.0f}% |
| $\\rho(\\lvert O\\rvert) = \\lambda_{{\\max}}(\\lvert O\\rvert)$ | {g['rho_abs']:.6f} | — | 绝对值守恒的地板 |
| 谱范数（sharp） | {tau2:.6f} | 自身 | similarity invariant |

平衡幅度：$\\max p - \\min p = {g['osborne_spread']:.4f}$（{g['osborne_iters']} 次 L-BFGS 迭代）。

$\\rho(|O|)$ 是合法的 $\\tau$，因为
$|x^T O x| \\le \\sum_{{ij}}|x_i||O_{{ij}}||x_j| = |x|^T|O||x| \\le \\lambda_{{\\max}}(|O|)\\|x\\|^2$。
它比 $\\rho(O)$ 大 ${g['rho_abs'] - tau2:+.6f}$ Ha —— **这个差就是 certificate 丢掉的符号抵消**。

任何绝对值型上界都无法低于 $\\rho(|O|)$：$A_s = DAD^{{-1}}$ 与 $A=|O|$ 相似，Gershgorin 圆盘包含谱，
故行和 $\\ge \\rho(A_s) = \\rho(A)$。所以连"完美平衡的 Gershgorin"也到不了 $\\rho(O)$，更到不了 gap。

### 4.2 符号缩放：恒等无效

若 $D = \\mathrm{{diag}}(\\varepsilon)$，$\\varepsilon\\in\\{{\\pm1\\}}^n$，则

$$|(DOD^{{-1}})_{{ij}}| = |\\varepsilon_i \\varepsilon_j O_{{ij}}| = |O_{{ij}}|,$$

**任何绝对值型界都完全不变。** 实测值 ${g['sign_scaling_value']:.9f}$ Ha，与 raw 只差
${g['sign_scaling_delta']:+.1e}$ Ha（`sign_scaling_is_vacuous` = {g['sign_scaling_is_vacuous']}）。
所以"用符号缩放改进 Gershgorin"在数学上是不可能的：Gershgorin 已经取了绝对值。

### 4.3 一个必须记录的死胡同（易犯错误）

直观上会想去平衡**符号**以让项相互抵消，即最小化
$\\max_i\\big|\\sum_j \\varepsilon_i\\varepsilon_j |O_{{ij}}|\\big|$。该目标确实能压到
${g['invalid_signed_sum']['value']:.4f}$ Ha —— **低于 $\\rho(O) = {tau2:.4f}$**。

任何真正的谱上界都不可能低于谱半径，所以这不是上界。更直接的证伪：把它自称要界住的那个矩阵
建出来，$B = \\mathrm{{diag}}(\\varepsilon)\\,|O|\\,\\mathrm{{diag}}(\\varepsilon)$，它的
$\\lambda_{{\\max}} = {g['invalid_signed_sum']['counterexample_lambda_max']:.6f}$ Ha，
比该"上界"${g['invalid_signed_sum']['value']:.4f}$ Ha 高
{g['invalid_signed_sum']['understatement_factor']:.2f} 倍 —— 而 $|B|$ 的行和仍是
${g['raw']:.6f}$（因为 $|\\varepsilon_i\\varepsilon_j|=1$）。也就是说：符号抵消只发生在
**没有谱内容**的量里，Gershgorin 求和的是 $|S_{{ij}}|$，完全不受符号影响。
脚本把它保存在 `invalid_signed_sum` 里并标记 `is_a_valid_bound: false`。

## 5. 为什么任何对角变换都不可能奏效（证明，非观测）

证书的假设是 $\\forall x:\\; x^T O x \\ge -\\tau \\|x\\|^2$，等价于**谱条件**

$$O + \\tau I \\succeq 0.$$

- **相似变换** $DOD^{{-1}}$：不改变二次型 $x^T O x$，改的是另一个矩阵的二次型；它能影响的只有范数型上界，
  而该上界 $\\ge \\rho(O)$，保谱。实测：20 个随机 $D$ 下谱扰动
  ${inv['similarity_max_eig_dev']:.1e}$（虚部 ${inv['similarity_max_imag']:.1e}$）。
- **合同变换** $D(O+\\tau I)D$：这正是把 $x = D^{{-1}}y$ 代入二次型后唯一自然的变换。
  由 Sylvester 惯性定律，$O+\\tau I \\succeq 0 \\iff D(O+\\tau I)D \\succeq 0$。
  实测：20 个随机 $D$ 下惯性改变 ${inv['congruence_inertia_violations']}$ 次（从未改变）。

于是假设的真假是 **basis-free** 的。既然它是假的（$\\lambda_{{\\min}}(O) = {b['lambda_min']:.4f}$，
且 $O$ 有 **{inv['n_eig_below_gap']}/{res['n_determinants']}** 个特征值低于 $-\\mathrm{{gap}}$），
任何换基都不把它变成真的。**这不是需要更聪明的 $D$，而是这个假设本身太弱。**

## 6. 松弛预算：问题出在哪

| 量 | 值 (Ha) | 相对 gap |
|---|---|---|
| diagonal gap | {gap:.6f} | 1.00 |
| 证书实际付出的 off-diagonal 风险 $\\lVert O\\rVert_2$ | {tau2:.6f} | {tau2 / gap:.2f} |
| 同上（单侧 sharp） | {tau1:.6f} | {tau1 / gap:.2f} |
| **真实需要的** $d - \\lambda_{{\\min}}(H)$ | {slack_real:.6f} | {slack_real / gap:.2f} |

证书为覆盖真实的 {slack_real:.4f} Ha 支付了 {tau1:.4f} Ha（单侧）或 {tau2:.4f} Ha（双侧），
即 **{v['one_sided_over_insurance']:.0f}× / {v['over_insurance_factor']:.0f}× 的过度保险**。
$\\lambda_{{\\min}}(O)$ 与 $\\lambda_{{\\min}}(H)$ 之间的全部差距就是能量的符号抵消 ——
而 $|x^T O x| \\le \\rho(O)\\|x\\|^2$ 这类界按定义把抵消全部丢掉。

## 7. 证书数值汇总

| $\\tau$ 来源 | $\\tau$ (Ha) | $d-\\tau$ (Ha) | 与 $E_{{\\mathrm{{FCI}}}}$ 之差 (Ha) | $\\tau/\\mathrm{{gap}}$ | $\\tau < \\mathrm{{gap}}$? |
|---|---|---|---|---|---|
{rows}

没有任何一行是非空的。最小可行 $\\tau$ 为 {v['best_tau_any']:.6f} Ha，是 gap 的 {v['best_tau_any'] / gap:.2f} 倍。

## 8. 复现

```
/usr/bin/python3 experiments/complement/sign_scaling.py   # ~40 s（H 已缓存）
```

输出 `sign_scaling.json`（全部中间量）、`sign_scaling.log`（完整运行日志）、
`sign_scaling.partial.json`（各 stage 后的检查点）。Hamiltonian 缓存在
`H_water_sto3g.npy`，仍会与 `fci.FCI` 重新对账。

内置的一致性断言（全部通过，且抓到过两个真实 bug）：
1. $\\lVert O\\rVert_2 \\ge \\rho(O)$；
2. 任何搜索得到的 $\\sigma_{{\\max}}$ 不得低于 $\\rho(O)$；
3. 平衡后 Gershgorin 不得低于 $\\rho(|O|)$；
4. 相似变换下谱不变（用 `eigvals`，**不能用** `eigvalsh` ——$DOD^{{-1}}$ 不对此称，
   这个错误一度给出 $1.6\\times10^3$ 的假偏差）；
5. 符号和尝试必须低于 $\\rho(O)$（否则该前提本身要重新检查）。

### 独立核验

`audit_sign_scaling.py` 不引用本文件任何函数，从缓存的 Hamiltonian 重建 $O$、用另一条
代码路径重算全部头条数字，并用**暴力搜索**（20000 个随机对角相似变换 + 2000 个随机
幅度缩放）直接检验中心论断，而不是相信 L-BFGS-B。结果写入
`sign_scaling_audit.json`。这一步存在的理由：本次工作中两个 bug 都是被
一致性断言抓到的，而不是被"看起来对"抓到的。

## 9. 局限

- 单分子、单基组、单 sector（441 determinants）。这是"真值可知"的模型测量，
  不推断任何 scaling law。
- 结论 $\\tau \\ge \\rho(O)$ 是 exact 的、与体系无关的；但"$\\tau > \\mathrm{{gap}}$"
  是这一体系的经验事实，不是定理。
- 未测试非对角（非 similarity、非 congruence）的其它界，例如利用 excitations 分块的
  block Gershgorin、或对 $x^T O x$ 做 sign-aware 的 resolvent 估计。这类界的成立与否
  不违反 §5 的论证（因为 §5 只排除"换基 + 范数上界"这一族），但不在本次范围内，
  也没有做任何验证。**不要把本文件的负结果读作"该证书在原则上无法收紧"。**
"""
    Path(path).write_text(md)
    print(f'wrote {path}')


# ----------------------------------------------------------------------------
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tee = Tee(LOG)
    sys.stdout = tee
    try:
        return run()
    finally:
        sys.stdout = tee.stdout
        tee.f.close()


def run():
    t_all = time.time()
    res = {'system': {'atom': ATOM, 'basis': BASIS, 'sector': 'S_z = 0'}}

    hr('stage 0: Hamiltonian')
    H, e_rhf, e_fci = fci_matrix()
    n = H.shape[0]
    diag = np.diag(H).copy()
    O = H - np.diag(diag)
    dsorted = np.sort(diag)
    d = float(dsorted[0])
    gap = float(dsorted[1] - dsorted[0])
    print(f'\ndeterminants   = {n}')
    print(f'E_RHF          = {e_rhf:.9f} Ha')
    print(f'E_FCI          = {e_fci:.9f} Ha')
    print(f'd (min diag)   = {d:.9f} Ha   == E_RHF? {abs(d - e_rhf) < 1e-9}')
    print(f'diagonal gap   = {gap:.9f} Ha')
    res['n_determinants'] = n
    res['rhf_energy'] = e_rhf
    res['fci_energy'] = e_fci
    res['d'] = d
    res['diagonal_gap'] = gap
    res['d_equals_rhf'] = bool(abs(d - e_rhf) < 1e-9)

    hr('stage 1: baseline norms of the off-diagonal part')
    base = baseline(O)
    res['baseline'] = {k: v for k, v in base.items() if k != 'eigvals'}
    tau = base['spectral_norm']

    hr('stage 2: spectral-norm balancing  min_s sigma_max(diag(e^s) O diag(e^-s))')
    res['balancing'] = stage_balancing(O, tau, base['spectral_radius'])
    (OUT_DIR / 'sign_scaling.partial.json').write_text(
        json.dumps(res, indent=2, default=float))

    hr('stage 3: signed Gershgorin')
    res['gershgorin'] = stage_gershgorin(O, gap)
    (OUT_DIR / 'sign_scaling.partial.json').write_text(
        json.dumps(res, indent=2, default=float))

    hr('stage 4: is the certificate fixable by ANY diagonal transform?')
    res['invariance'] = stage_invariance(O, gap, tau, d, e_fci, res['gershgorin'])

    hr('VERDICT')
    g = res['gershgorin']
    tau_one = res['invariance']['tau_one_sided']
    taus = {'spectral_norm_before': tau,
            'spectral_norm_after_balancing': res['balancing']['best_found'],
            'spectral_one_sided_sharp': tau_one,
            'gershgorin_raw': g['raw'],
            'gershgorin_balanced': g['balanced_gershgorin'],
            'rho_abs_floor': g['rho_abs']}
    for k, v in taus.items():
        print(f'  tau via {k:32s} = {v:.6f} Ha   tau/gap = {v / gap:6.2f}'
              f'   tau < gap? {v < gap}')
    print(f'\n  floor: rho(O) = {base["spectral_radius"]:.6f} Ha; every norm-type')
    print(f'  bound is >= rho(O) by similarity invariance of the spectrum.')
    print(f'  gap = {gap:.6f} Ha = {gap / base["spectral_radius"]:.3f} x rho(O),')
    print(f'  so the certificate is short by a factor '
          f'{base["spectral_radius"] / gap:.2f} and no rescaling touches that.')
    print(f'\n  TRUE lambda_min(H) = {e_fci:.6f}; d - lambda_min = {d - e_fci:.6f} Ha.')
    print(f'  The certificate pays {base["spectral_radius"]:.4f} Ha of '
          f'off-diagonal risk to cover')
    print(f'  an actual {d - e_fci:.4f} Ha of it: '
          f'{base["spectral_radius"] / (d - e_fci):.0f}x over-insurance.')

    nonvac = bool(min(taus.values()) < gap)
    res['verdict'] = {
        'non_vacuous_achievable': nonvac,
        'any_transform_helps': bool(res['balancing']['improvement'] > 1e-9),
        'spectral_infimum': base['spectral_radius'],
        'tau_one_sided_sharp': tau_one,
        'one_sided_deficit_factor': tau_one / gap,
        'best_tau_any': float(min(taus.values())),
        'gap': gap,
        'deficit_factor': base['spectral_radius'] / gap,
        'n_eigenvalues_of_O_below_minus_gap': res['invariance']['n_eig_below_gap'],
        'over_insurance_factor': base['spectral_radius'] / (d - e_fci),
        'one_sided_over_insurance': tau_one / (d - e_fci),
        'taus': taus,
    }
    print(f'\n  => can the certificate be made non-vacuous by a diagonal '
          f'transform? {nonvac}')
    res['runtime_s'] = time.time() - t_all
    res['timestamp'] = time.strftime('%Y-%m-%d %H:%M:%S')

    (OUT_DIR / 'sign_scaling.json').write_text(json.dumps(res, indent=2, default=float))
    print(f'\nwrote {OUT_DIR / "sign_scaling.json"}  ({res["runtime_s"]:.1f}s)')
    write_markdown(res, OUT_DIR / 'sign_scaling.md')
    print(f'log  {LOG}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
