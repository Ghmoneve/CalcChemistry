"""Second system for the certificate measurement: N2 / STO-3G.

WHY A SECOND SYSTEM
Every certificate number so far comes from water / STO-3G, 441 determinants. A
reader may ask whether the 0.860 kcal/mol Schur bound and the closure of the block
route are properties of the method or of that one matrix. This repeats the
measurement on N2 in the same basis: 14400 determinants (33x water), a triple bond
rather than lone pairs, and a smaller diagonal gap.

MEMORY DISCIPLINE -- why the earlier attempt was killed
It built H, then O = H - diag(D) as a second full matrix, then A = H[rest][:,rest]
as a third, reaching ~5 GB and being OOM-killed. Nothing here holds more than one
ndet x ndet array:

  pass 1  build H (1.66 GB), read d, the reference column b, and the extremal
          eigenvalues. The off-diagonal part is applied as a LinearOperator
          v -> H@v - D*v, so O is never materialised. Then H is freed.
  pass 2  build A = Hr - dI directly into one array; H is not rebuilt.
          Factor in place with LAPACK dpotrf and solve with dpotrs, so no second
          copy of the factor exists either.

Columns are recomputed in pass 2 rather than checkpointed to disk: one contract_2e
sweep costs ~20 s and a .npy round trip costs more in complexity than it saves.

WHY ARPACK AND NOT eigvalsh
Only extremal eigenvalues are needed. Ritz values are used for reporting only; the
one quantity that must be rigorous, the sign of lambda_min(Hr - dI), is decided by
the Cholesky certificate, not by ARPACK.

SCOPE  A model measurement where the exact answer is available.
"""
from __future__ import annotations

import json
import os
import time
from math import comb
from pathlib import Path

import numpy as np
from pyscf import ao2mo, fci, gto, scf
from scipy.linalg import lapack
from scipy.sparse.linalg import LinearOperator, eigsh

KCAL = 627.5094740631
U = 2.0 ** -53
ATOM = os.environ.get('SYS_ATOM', 'N 0 0 0; N 0 0 1.0977')
BASIS = os.environ.get('SYS_BASIS', 'sto-3g')
OUT = Path(os.environ.get('SYS_OUT', 'experiments/complement/certificate_n2.json'))
LOG = Path(os.environ.get('SYS_LOG', 'experiments/complement/certificate_n2.log'))


def rss_gb():
    try:
        for line in Path('/proc/self/status').read_text().splitlines():
            if line.startswith('VmRSS:'):
                return int(line.split()[1]) / 1e6
    except OSError:
        pass
    return float('nan')


def say(msg):
    line = f'[{time.strftime("%H:%M:%S")}] {msg}  [rss {rss_gb():.2f} GB]'
    print(line, flush=True)
    with LOG.open('a') as fh:
        fh.write(line + '\n')


def main():
    LOG.write_text('')
    say('start')
    mol = gto.M(atom=ATOM, basis=BASIS, charge=0, spin=0, verbose=0)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-11
    mf.kernel()
    if not mf.converged:
        raise SystemExit('RHF not converged')
    norb = mf.mo_coeff.shape[1]
    nelec = mol.nelectron
    na, nb = nelec // 2 + nelec % 2, nelec // 2
    h1e = mf.mo_coeff.T @ mf.get_hcore() @ mf.mo_coeff
    eri = ao2mo.restore(1, ao2mo.kernel(mol, mf.mo_coeff), norb)
    h2e = fci.direct_spin1.absorb_h1e(h1e, eri, norb, (na, nb), 0.5)
    nda, ndb = comb(norb, na), comb(norb, nb)
    ndet = nda * ndb
    enuc = mol.energy_nuc()
    say(f'RHF={mf.e_tot:.9f} norb={norb} nelec={nelec} determinants={ndet}')

    def col(j):
        u = np.zeros(ndet)
        u[j] = 1.0
        return np.asarray(fci.direct_spin1.contract_2e(
            h2e, u.reshape(nda, ndb), norb, (na, nb))).ravel()

    # ---- pass 1 ----
    t0 = time.time()
    H = np.empty((ndet, ndet), order='C')
    for j in range(ndet):
        H[:, j] = col(j)
        if j % 4000 == 0:
            say(f'  H col {j}/{ndet} ({time.time()-t0:.0f}s)')
    H += enuc * np.eye(ndet)
    say(f'H built ({time.time()-t0:.0f}s)')

    e_ref = fci.FCI(mol, mf.mo_coeff).kernel()[0]
    ev = np.sort(eigsh(H, k=3, which='SA', return_eigenvectors=False))
    say(f'fci.FCI={e_ref:.9f} eigsh={ev[0]:.9f} diff={ev[0]-e_ref:+.2e}')
    if abs(ev[0] - e_ref) > 1e-6:
        raise SystemExit('eigsh disagrees with fci.FCI')

    D = np.diag(H).copy()
    iref = int(np.argmin(D))
    rest = np.array([i for i in range(ndet) if i != iref])
    b = H[rest, iref].copy()
    d = float(D[iref])
    Ds = np.sort(D)
    gap = float(Ds[1] - Ds[0])
    op = LinearOperator((ndet, ndet), matvec=lambda v: H @ v - D * v, dtype=float)
    tau = max(
        float(np.abs(eigsh(op, k=1, which='LA', return_eigenvectors=False)[0])),
        float(np.abs(eigsh(op, k=1, which='SA', return_eigenvectors=False)[0])),
    )
    del op
    say(f'd={d:.6f} gap={gap:.6f} ||O||={tau:.6f} tau/gap={tau/gap:.2f}')

    res = {
        'system': {'atom': ATOM, 'basis': BASIS},
        'n_determinants': ndet,
        'rhf_energy': float(mf.e_tot),
        'fci_energy': float(ev[0]),
        'fci_reference': float(e_ref),
        'lambda2': float(ev[1]),
        'true_gap': float(ev[1] - ev[0]),
        'diagonal_min': d,
        'diagonal_gap_to_next': gap,
        'correlation_energy': float(d - ev[0]),
        'reference_index': iref,
        'tau_spectral': tau,
        'tau_over_gap': tau / gap,
        'block_bound': d - tau,
        'block_slack_kcal': (d - tau - ev[0]) * KCAL,
    }
    OUT.write_text(json.dumps(res, indent=2))
    say(f'block d-tau={d-tau:.6f} slack {(d-tau-ev[0])*KCAL:+.1f} kcal/mol')
    H = None
    del H
    say('H released')

    # ---- pass 2 ----
    nA = rest.size
    A = np.empty((nA, nA), order='C')
    t0 = time.time()
    for jj, j in enumerate(rest):
        A[:, jj] = col(j)[rest]
        if jj % 4000 == 0:
            say(f'  A col {jj}/{nA} ({time.time()-t0:.0f}s)')
    A[np.diag_indices_from(A)] += enuc - d
    say(f'A built ({time.time()-t0:.0f}s)')

    # trace(A) = ||L||_F^2, captured before dpotrf overwrites A; this is all the
    # rounding bound needs, so |L| and |L||L^T| are never materialised.
    trace_A = float(np.trace(A))
    t0 = time.time()
    L, info = lapack.dpotrf(A, lower=1, clean=0, overwrite_a=1)
    say(f'dpotrf info={info} ({time.time()-t0:.0f}s)')
    if info != 0:
        res['certified_positive_definite'] = False
        res['dpotrf_info'] = int(info)
        OUT.write_text(json.dumps(res, indent=2))
        say('NOT positive definite')
        return

    # dpotrf(lower=1) leaves the strictly upper triangle holding the original A,
    # so every later use of L (dtrtri, the residual, np.abs) must see zeros there.
    L[np.triu_indices_from(L, 1)] = 0.0

    # || |L||L^T| ||_F <= ||L||_F^2 = trace(A) by submultiplicativity of the
    # Frobenius norm, so the rounding bound costs no extra allocation.
    eps_round = float(nA * U * trace_A / (1.0 - nA * U))

    y, info2 = lapack.dpotrs(L, b, lower=1)
    if info2 != 0:
        raise SystemExit(f'dpotrs info={info2}')
    # residual b - A y, evaluated as b - L(L^T y) because A = L L^T
    resid = b - L @ (L.T @ y)
    centre = float(b @ y)

    # Now L is no longer needed: invert it in place to bound sigma_min(L).
    Linv, info3 = lapack.dtrtri(L, lower=1, overwrite_c=1)
    if info3 != 0:
        raise SystemExit(f'dtrtri info={info3}')
    rowmax = 0.0
    for i in range(nA):
        rowmax = max(rowmax, float(np.abs(Linv[i]).sum()))
    sigma_lb = 1.0 / (np.sqrt(nA) * rowmax)
    delta = sigma_lb ** 2 - eps_round
    say(f'sigma^2>={sigma_lb**2:.4e} round<={eps_round:.2e} delta>={delta:.4e}')
    rad = float(np.linalg.norm(b) * np.linalg.norm(resid) / delta) if delta > 0 else float('inf')
    schur = d - (centre + rad)
    res.update({
        'certified_positive_definite': bool(delta > 0),
        'dpotrf_info': 0,
        'sigma_min_L_lower_bound': sigma_lb,
        'rounding_bound': eps_round,
        'delta_certified': delta,
        'schur_scalar': centre,
        'schur_scalar_radius': rad,
        'schur_bound': schur,
        'schur_slack_kcal': (schur - ev[0]) * KCAL,
        'improvement_factor': float((d - tau - ev[0]) / (schur - ev[0]))
        if schur != ev[0] else None,
    })
    OUT.write_text(json.dumps(res, indent=2))
    say(f'b^T A^-1 b={centre:.10f} +/- {rad:.2e}')
    say(f'Schur={schur:.8f} slack {(schur-ev[0])*KCAL:+.3f} kcal/mol')
    if schur != ev[0]:
        say(f'improvement {(d-tau-ev[0])/(schur-ev[0]):.0f}x')
    say('done')


if __name__ == '__main__':
    main()
