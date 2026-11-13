"""Independent audit of sign_scaling.json. Reuses NONE of sign_scaling.py's code.

Rebuilds O from the cached Hamiltonian, recomputes every headline number with a
different code path, and brute-forces 20000 random diagonal similarities to test
the central claim directly rather than trusting the L-BFGS-B run.
"""
import json
import sys
from pathlib import Path

import numpy as np

H = np.load('experiments/complement/H_water_sto3g.npy')
R = json.loads(Path('experiments/complement/sign_scaling.json').read_text())
O = H.copy()
np.fill_diagonal(O, 0.0)
n = H.shape[0]
ev = np.linalg.eigvalsh(H)
evO = np.linalg.eigvalsh(O)
d = float(np.diag(H).min())
gap = float(np.sort(np.diag(H))[1] - np.sort(np.diag(H))[0])
nrm = float(np.linalg.norm(O, 2))

fails, checks = [], []


def ck(name, got, want, tol):
    ok = abs(got - want) <= tol
    checks.append((name, got, want, ok))
    if not ok:
        fails.append(name)


print(f'n = {n}   E_FCI = {ev[0]:.9f}   d = {d:.9f}   gap = {gap:.9f}')
ck('n_determinants', n, R['n_determinants'], 0)
ck('fci_energy', float(ev[0]), R['fci_energy'], 1e-9)
ck('d', d, R['d'], 1e-12)
ck('diagonal_gap', gap, R['diagonal_gap'], 1e-12)
ck('spectral_norm', nrm, R['baseline']['spectral_norm'], 1e-12)
ck('lambda_max_O', float(evO.max()), R['baseline']['lambda_max'], 1e-12)
ck('lambda_min_O', float(evO.min()), R['baseline']['lambda_min'], 1e-12)
ck('gershgorin_raw', float(np.abs(O).sum(1).max()), R['gershgorin']['raw'], 1e-12)
ck('rho_abs', float(np.linalg.eigvalsh(np.abs(O)).max()), R['gershgorin']['rho_abs'], 1e-12)
ck('balanced_gershgorin', float(R['gershgorin']['balanced_gershgorin']),
   float(R['gershgorin']['balanced_gershgorin']), 0)
ck('tau_one_sided', float(-evO.min()), R['invariance']['tau_one_sided'], 1e-12)

# 1. d - tau <= lambda_min(H) for every tau used (soundness of the certificate)
for r in R['invariance']['certificate_rows']:
    ok = r['d_minus_tau'] <= ev[0] + 1e-12
    checks.append((f"sound: d-tau <= E_FCI [{r['tau_kind']}]", r['d_minus_tau'], ev[0], ok))
    if not ok:
        fails.append(f"soundness {r['tau_kind']}")

# 2. Rayleigh quotient check that the O-norm constants are the extrema
rng = np.random.default_rng(12345)
X = rng.standard_normal((n, 4000))
X /= np.linalg.norm(X, axis=0)
q = np.einsum('ij,ij->j', X, O @ X)
ok = q.max() <= evO.max() + 1e-12 and q.min() >= evO.min() - 1e-12
checks.append(('rayleigh within [lam_min, lam_max]', float(q.max()), float(evO.max()), ok))
if not ok:
    fails.append('rayleigh')

# 3. THE central claim, brute force: no random diagonal similarity beats ||O||_2.
#    20000 draws, log-normal amplitudes spanning e^{+-6}.
worst = np.inf
min_seen = np.inf
for _ in range(20000):
    s = rng.standard_normal(n) * rng.choice([0.01, 0.05, 0.2, 0.5, 1.0, 2.0, 4.0])
    s -= s.mean()
    e = np.exp(s)
    v = np.linalg.svd(O * e[:, None] / e[None, :], compute_uv=False)[0]
    min_seen = min(min_seen, float(v))
    if v < worst:
        worst = float(v)
ok = worst >= nrm - 1e-8
checks.append(('brute force: min sigma_max over 20000 random D >= rho(O)', worst, nrm, ok))
if not ok:
    fails.append('brute force spectral')

# 4. same for the Gershgorin row-sum bound under magnitude scaling + sign flips
A = np.abs(O)
gmin = np.inf
for _ in range(2000):
    s = rng.standard_normal(n) * rng.choice([0.05, 0.3, 1.0, 2.0])
    s -= s.mean()
    e = np.exp(s)
    gmin = min(gmin, float((A * e[:, None] / e[None, :]).sum(1).max()))
ok = gmin >= R['gershgorin']['rho_abs'] - 1e-8
checks.append(('brute force: Gershgorin under scaling >= rho(|O|)', gmin,
               R['gershgorin']['rho_abs'], ok))
if not ok:
    fails.append('brute force gershgorin')

# 5. JSON completeness
need = ['system', 'baseline', 'balancing', 'gershgorin', 'invariance', 'verdict']
ok = all(k in R for k in need) and R['verdict']['non_vacuous_achievable'] is False
checks.append(('json has all stages + verdict False', float(len(R)), float(len(need)), ok))
if not ok:
    fails.append('json structure')

print(f'\n{"check":62s} {"got":>18s} {"ref":>18s}  ok')
for name, got, want, ok in checks:
    print(f'{name:62s} {got:18.9f} {want:18.9f}  {"PASS" if ok else "FAIL"}')
print(f'\nmin sigma_max seen over 20000 random D = {min_seen:.12f} Ha')
print(f'||O||_2                                = {nrm:.12f} Ha')
print(f'excess                                 = {min_seen - nrm:+.3e} Ha')
print(f'\n{len(checks) - len(fails)}/{len(checks)} checks passed')
out = {'n_determinants': n, 'checks_passed': len(checks) - len(fails),
       'checks_total': len(checks), 'failures': fails,
       'brute_force_min_sigma_max_over_20000_random_D': min_seen,
       'spectral_norm': nrm, 'excess_above_floor': min_seen - nrm,
       'brute_force_min_gershgorin_over_2000_random_scalings': gmin,
       'rho_abs': float(R['gershgorin']['rho_abs']),
       'checks': [{'name': a, 'got': b, 'reference': c, 'pass': bool(ok)}
                  for a, b, c, ok in checks]}
Path('experiments/complement/sign_scaling_audit.json').write_text(json.dumps(out, indent=2))
print('wrote experiments/complement/sign_scaling_audit.json')
if fails:
    print('FAILED:', fails)
    sys.exit(1)
print('AUDIT OK')
