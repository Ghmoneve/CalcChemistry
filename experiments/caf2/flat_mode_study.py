"""Numerical study of the flat Ca-F-Ca bending mode: why damped BFGS stalls,
what the analytic Hessian buys, and what the basis set and grid do to it.

Sub-commands (each appends JSON so partial runs survive):
  profile   E(theta) at fixed r, plus projected Hessians along the path
  bfgs      instrumented replay of experiments/caf2/reoptimize.py:relax()
  seed      BFGS with the analytic Hessian as the initial inverse Hessian
  newton    trust-region Newton / RFO with the analytic Hessian each step
  basis     relaxed theta-scan for several basis sets and functionals
  grid      grid-level dependence of E(theta) and of the Hessian
  linear    C_infinity_v constrained optimisation + 1D bend fit
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from pyscf import dft, gto
from pyscf.hessian import rks as hrks

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'experiments/caf2/flat_mode.json'
BOHR = 0.52917721092
CONV = 5140.02            # cm^-1 per sqrt(Ha/(bohr^2 amu))
HARTREE2KCAL = 627.5094740631
ZS = [20, 9, 9]


# ---------------------------------------------------------------- geometry ---
def geom(theta_deg, r=1.9925):
    """Symmetric bent CaF2, Ca at origin, both Ca-F = r, F-Ca-F = theta."""
    th = np.radians(theta_deg)
    return np.array([[0.0, 0.0, 0.0],
                     [r * np.sin(th / 2), 0.0, r * np.cos(th / 2)],
                     [-r * np.sin(th / 2), 0.0, r * np.cos(th / 2)]])


def angle_of(x):
    v1 = x[1] - x[0]
    v2 = x[2] - x[0]
    c = float(v1 @ v2 / np.linalg.norm(v1) / np.linalg.norm(v2))
    return float(np.degrees(np.arccos(np.clip(c, -1.0, 1.0))))


def bonds_of(x):
    return [float(np.linalg.norm(x[1] - x[0])), float(np.linalg.norm(x[2] - x[0]))]


def build(x, basis='def2-svp', charge=0, spin=0):
    return gto.M(atom=[(ZS[i], tuple(x[i])) for i in range(3)],
                 basis=basis, spin=spin, charge=charge, unit='Angstrom',
                 verbose=0)


def scf(x, xc='pbe', basis='def2-svp', grid=4, conv_tol=1e-11):
    mol = build(x, basis)
    attempts = []
    for damp, init in ((0.0, None), (0.3, None), (0.0, 'atom'), (0.5, 'atom')):
        mf = dft.RKS(mol)
        mf.xc = xc
        mf.grids.level = grid
        mf.conv_tol = conv_tol
        if damp:
            mf.damp = damp
        if init:
            mf.init_guess = init
        try:
            mf.kernel()
        except Exception as exc:                        # noqa: BLE001
            attempts.append(f'{damp}/{init}:{type(exc).__name__}')
            continue
        if mf.converged:
            return float(mf.e_tot), mf
        attempts.append(f'{damp}/{init}:noconv')
    raise RuntimeError(f'SCF failed {attempts}')


def energy_grad(x, **kw):
    e, mf = scf(x, **kw)
    g = mf.nuc_grad_method().kernel().reshape(-1) / BOHR   # Ha/Angstrom
    return e, g, mf


# ------------------------------------------------- hessian + normal modes ----
def hessian(x, **kw):
    e, mf = scf(x, **kw)
    h = np.array(hrks.Hessian(mf).kernel())
    H = h.transpose(0, 2, 1, 3).reshape(9, 9)
    return 0.5 * (H + H.T), e, mf


def projector(mol, x, linear=None):
    """Projector onto the 3N-6 (3N-5 if linear) internal space, in mass-weighted
    coordinates.  Returns (P, n_null, null_vectors)."""
    masses = np.array(mol.atom_mass_list())
    sq = np.repeat(np.sqrt(masses), 3)
    if linear is None:
        linear = angle_of(x) > 179.5
    cands = []
    for a in range(3):                                   # translations
        v = np.zeros(9)
        for i in range(3):
            v[3 * i + a] = sq[3 * i]
        cands.append(v)
    for a in range(3):                                   # rotations
        e = np.zeros(3)
        e[a] = 1.0
        v = np.zeros(9)
        for i in range(3):
            v[3 * i:3 * i + 3] = sq[3 * i] * np.cross(e, x[i])
        cands.append(v)
    bv = []
    for v in cands:
        nrm = np.linalg.norm(v)
        if nrm < 1e-6:
            continue
        v = v / nrm
        for u in bv:
            v = v - float(v @ u) * u
        nrm = np.linalg.norm(v)
        if nrm > 1e-4:
            bv.append(v / nrm)
    P = np.eye(9) - sum(np.outer(u, u) for u in bv)
    return P, len(bv), bv, sq


def modes(x, **kw):
    """Physical normal modes after projecting translations+rotations."""
    H, e, mf = hessian(x, **kw)
    P, nnull, bv, sq = projector(mf.mol, x)
    Hm = H / np.outer(sq, sq)
    Hp = P @ Hm @ P
    w, vecs = np.linalg.eigh(0.5 * (Hp + Hp.T))
    keep = np.argsort(-np.abs(w))[: 9 - nnull]
    keep = keep[np.argsort(w[keep])]
    freqs = np.sign(w[keep]) * CONV * np.sqrt(np.abs(w[keep]))
    return dict(freqs=[float(f) for f in freqs],
                w=[float(v) for v in w[keep]],
                vecs=vecs[:, keep], e=e, H=H, P=P, sq=sq)


def raw_rot_eigs(x, **kw):
    """Largest mass-weighted eigenvalues of the *unprojected* Hessian: for a
    linear molecule these should be zero but the atom-centred grid breaks
    rotational symmetry."""
    H, e, mf = hessian(x, **kw)
    masses = np.array(mf.mol.atom_mass_list())
    sq = np.repeat(np.sqrt(masses), 3)
    w = np.linalg.eigvalsh(H / np.outer(sq, sq))
    return float(np.max(np.abs(w))), float(np.min(np.abs(w)))


# ------------------------------------------------------------------ utils ---
def dump(key, val):
    state = {}
    if OUT.exists():
        state = json.loads(OUT.read_text())
    state[key] = val
    OUT.write_text(json.dumps(state, indent=2, default=float))
    print(f'  -> wrote {key}')


def say(msg):
    print(f'[{time.strftime("%H:%M:%S")}] {msg}', flush=True)


# -------------------------------------------------------------- profile -----
def cmd_profile():
    basis = sys.argv[2] if len(sys.argv) > 2 else 'def2-svp'
    rows = []
    for th in [180, 175, 170, 165, 160, 155, 150, 145, 140, 135, 130, 125, 120, 110, 100]:
        x = geom(th)
        e, g, mf = energy_grad(x, basis=basis)
        rows.append(dict(theta=th, E=e, gnorm=float(np.linalg.norm(g)),
                         r=bonds_of(x), ang=angle_of(x)))
        say(f'theta {th:5.1f}  E {e:.8f}  |g| {np.linalg.norm(g):.2e}')
    dump(f'profile_{basis}', rows)

    hs = []
    for th in [180, 160, 140, 120]:
        x = geom(th)
        m = modes(x, basis=basis)
        P, nnull, bv, sq = projector(build(x, basis), x)
        # G_theta_theta = |dq/dtheta|^2 by finite difference in mass-weighted coords
        dq = ((geom(th + 1.0) - geom(th - 1.0)) / 2 * np.pi / 180).reshape(-1)
        dqm = dq * sq
        G = float(dqm @ dqm)
        hs.append(dict(theta=th, freqs=m['freqs'], w=m['w'],
                       E=m['e'], G_theta=G,
                       k_theta=float(m['w'][0] * G),
                       n_null=nnull))
        say(f'theta {th}: freqs {np.round(m["freqs"],1)}  '
            f'k_theta {m["w"][0]*G:.3e} Ha/rad^2')
    dump(f'hess_{basis}', hs)


# ----------------------------------------------------------------- bfgs -----
def cmd_bfgs():
    """Exact replay of reoptimize.py:relax, logging B_k eigenvalues."""
    x = geom(140.0)
    hinv = np.eye(9) * 0.3
    gp = xp = None
    traj = []
    e_prev = None
    for it in range(60):
        e, g, mf = energy_grad(x)
        gnorm = float(np.linalg.norm(g))
        m = modes(x)
        P, nnull, bv, sq = projector(mf.mol, x)
        # B_k eigenvalues inside the internal space
        hi = P @ hinv @ P
        wB = np.linalg.eigvalsh(hi)
        wB = wB[np.argsort(-np.abs(wB))][: 9 - nnull]
        # decomposition of g onto normal modes of the *true* Hessian
        comp = (P @ g) @ m['vecs']
        traj.append(dict(it=it, E=e, gnorm=gnorm, ang=angle_of(x),
                         r=bonds_of(x),
                         B_eig=[float(v) for v in np.sort(wB)],
                         H_eig=[float(v) for v in np.sort(m['w'])],
                         g_mode=[float(v) for v in comp],
                         g_internal=float(np.linalg.norm(P @ g)),
                         step=0.0))
        say(f'it {it:3d} E {e:.8f} |g| {gnorm:.3e} theta {angle_of(x):.2f} '
            f'B eig {np.array2string(np.sort(wB), precision=2)}')
        if gnorm < 1e-5:
            break
        if gp is not None:
            s = (x - xp).reshape(-1)
            y = g - gp
            sy = float(s @ y)
            if sy > 1e-12:
                rho = 1.0 / sy
                idn = np.eye(9)
                hinv = ((idn - rho * np.outer(s, y)) @ hinv
                        @ (idn - rho * np.outer(y, s)) + rho * np.outer(s, s))
        step = -hinv @ g
        ns = float(np.linalg.norm(step))
        if ns > 0.08:
            step *= 0.08 / ns
        traj[-1]['step'] = ns
        xp, gp = x.copy(), g.copy()
        x = x + step.reshape(3, 3)
        e_prev = e
    dump('bfgs_traj', traj)


def _bfgs_run(x0, hinv0, max_iter=60, cap=0.08, gtol=1e-5, tag='run'):
    x = np.array(x0, float)
    hinv = np.array(hinv0, float)
    gp = xp = None
    traj = []
    for it in range(max_iter):
        e, g, mf = energy_grad(x)
        gnorm = float(np.linalg.norm(g))
        traj.append(dict(it=it, E=e, gnorm=gnorm, ang=angle_of(x),
                         r=bonds_of(x)))
        say(f'{tag} it {it:3d} E {e:.9f} |g| {gnorm:.3e} '
            f'theta {angle_of(x):7.3f} r {np.mean(bonds_of(x)):.5f}')
        if gnorm < gtol:
            break
        if gp is not None:
            s = (x - xp).reshape(-1)
            y = g - gp
            sy = float(s @ y)
            if sy > 1e-14:
                rho = 1.0 / sy
                idn = np.eye(9)
                hinv = ((idn - rho * np.outer(s, y)) @ hinv
                        @ (idn - rho * np.outer(y, s)) + rho * np.outer(s, s))
        step = -hinv @ g
        ns = float(np.linalg.norm(step))
        if ns > cap:
            step *= cap / ns
        xp, gp = x.copy(), g.copy()
        x = x + step.reshape(3, 3)
    return traj, x, e, gnorm


def cmd_seed():
    """BFGS seeded with the analytic Hessian (pseudo-inverse, null space zeroed)."""
    x0 = geom(140.0)
    m = modes(x0)
    H = m['H']
    P, nnull, bv, sq = projector(build(x0), x0)
    Hm = H / np.outer(sq, sq)
    Hp = P @ Hm @ P
    w, V = np.linalg.eigh(Hp)
    keep = np.argsort(-np.abs(w))[: 9 - nnull]
    hinv = np.zeros((9, 9))
    for k in keep:
        hinv += np.outer(V[:, k], V[:, k]) / w[k]
    traj, x, e, gn = _bfgs_run(x0, hinv, tag='seedH')
    dump('bfgs_seeded', traj)
    m2 = modes(x)
    dump('bfgs_seeded_final', dict(E=e, gnorm=gn, coords=x.tolist(),
                                   ang=angle_of(x), r=bonds_of(x),
                                   freqs=m2['freqs']))


def cmd_newton():
    """Trust-region Newton / RFO with the analytic Hessian every step.

    Step solves (H + lambda I) s = -g with lambda chosen so that |s| = trust
    radius (level-shifted / Levenberg-Marquardt form of the RFO step).  Falls
    back to steepest descent along the lowest eigenvector when H has a negative
    eigenvalue and |g| is small there.
    """
    x = geom(140.0)
    traj = []
    for it in range(30):
        e, g, mf = energy_grad(x)
        gnorm = float(np.linalg.norm(g))
        m = modes(x)
        P, nnull, bv, sq = projector(mf.mol, x)
        Hm = m['H'] / np.outer(sq, sq)
        Hp = P @ Hm @ P
        w, V = np.linalg.eigh(Hp)
        keep = np.argsort(-np.abs(w))[: 9 - nnull]
        w, V = w[keep], V[:, keep]
        traj.append(dict(it=it, E=e, gnorm=gnorm, ang=angle_of(x),
                         r=bonds_of(x), H_eig=[float(v) for v in w]))
        say(f'newton it {it:2d} E {e:.9f} |g| {gnorm:.3e} '
            f'theta {angle_of(x):7.3f} r {np.mean(bonds_of(x)):.5f}')
        if gnorm < 1e-5:
            break
        gc = V.T @ (P @ g)
        lam, radius = 0.0, 0.10                     # Angstrom
        for _ in range(200):                        # bracket the shift
            s = V @ np.array([-gc[i] / (w[i] + lam) for i in range(len(w))])
            if np.linalg.norm(s) <= radius:
                break
            lam = max(1e-8, lam * 10 if lam else 1e-4)
            while np.linalg.norm(V @ np.array([-gc[i] / (w[i] + lam)
                                               for i in range(len(w))])) > radius:
                lam *= 2.0
                if lam > 1e6:
                    break
        step = s
        ns = np.linalg.norm(step)
        if ns > radius:
            step *= radius / ns
        x = x + step.reshape(3, 3)
    dump('newton_traj', traj)


# ----------------------------------------------------------------- basis ----
def relaxed_scan(basis, xc='pbe', grid=4, thetas=None):
    """For each theta, relax r by a 3-point parabola (two Newton passes)."""
    if thetas is None:
        thetas = [180, 170, 160, 150, 140, 130, 120, 110]
    rows = []
    r_guess = 2.00
    for th in thetas:
        r = r_guess
        for _ in range(3):
            es = []
            for dr in (-0.03, 0.0, 0.03):
                e, _ = scf(geom(th, r + dr), basis=basis, xc=xc, grid=grid)
                es.append(e)
            a, b, c = es
            d1 = (c - a) / 0.06
            d2 = (a - 2 * b + c) / 0.03 ** 2
            if d2 > 0:
                r = r - d1 / d2
                r = min(max(r, 1.80), 2.30)
        e, _ = scf(geom(th, r), basis=basis, xc=xc, grid=grid)
        rows.append(dict(theta=th, r=float(np.mean([r, r])), E=e))
        r_guess = r
        say(f'{basis}/{xc} theta {th:5.1f} r {r:.5f} E {e:.8f}')
    emin = min(rows, key=lambda z: z['E'])
    for z in rows:
        z['dE_kcal'] = (z['E'] - emin['E']) * HARTREE2KCAL
    return dict(basis=basis, xc=xc, grid=grid, rows=rows,
                theta_min=emin['theta'], r_at_min=emin['r'],
                bend_rise_kcal=(rows[0]['E'] - emin['E']) * HARTREE2KCAL)


def cmd_basis():
    cases = [('def2-svp', 'pbe'), ('def2-tzvp', 'pbe'), ('def2-tzvpd', 'pbe'),
             ('ma-def2-svp', 'pbe'), ('def2-svpd', 'pbe'), ('def2-tzvpp', 'pbe'),
             ('def2-svp', 'b3lyp'), ('def2-tzvpd', 'b3lyp'),
             ('def2-svp', 'pbe0'), ('def2-tzvpd', 'pbe0')]
    res = {}
    for basis, xc in cases:
        key = f'{basis}_{xc}'
        try:
            res[key] = relaxed_scan(basis, xc)
            dump(f'scan_{key}', res[key])
        except Exception as exc:                     # noqa: BLE001
            say(f'{key} FAILED {exc}')
            res[key] = {'error': str(exc)}


# ------------------------------------------------------------------ grid ----
def cmd_grid():
    rows = []
    ref = None
    for lvl in [0, 1, 2, 3, 4, 5, 6]:
        e180, _ = scf(geom(180), grid=lvl)
        e140, _ = scf(geom(140), grid=lvl)
        e160, _ = scf(geom(160), grid=lvl)
        m = modes(geom(180), grid=lvl)
        rmax, rmin = raw_rot_eigs(geom(180), grid=lvl)
        rows.append(dict(level=lvl, E180=e180, E140=e140, E160=e160,
                         dE140_180=(e140 - e180) * HARTREE2KCAL,
                         dE160_180=(e160 - e180) * HARTREE2KCAL,
                         bend_w=m['w'][0], bend_cm=m['freqs'][0],
                         raw_max_abs_eig=rmax, raw_min_abs_eig=rmin,
                         freqs=m['freqs']))
        say(f'grid {lvl}: dE(140-180) {(e140-e180)*HARTREE2KCAL:+.4f} kcal  '
            f'bend {m["freqs"][0]:+.2f} cm-1  raw|w|max {rmax:.2e}')
    # a dense, unpruned reference
    mol = build(geom(180))
    mf = dft.RKS(mol)
    mf.xc = 'pbe'
    mf.grids.atom_grid = (150, 770)
    mf.grids.prune = None
    mf.conv_tol = 1e-11
    mf.kernel()
    e_ref = float(mf.e_tot)
    h = np.array(hrks.Hessian(mf).kernel())
    H = h.transpose(0, 2, 1, 3).reshape(9, 9)
    H = 0.5 * (H + H.T)
    P, nnull, bv, sq = projector(mol, geom(180))
    w, V = np.linalg.eigh(P @ (H / np.outer(sq, sq)) @ P)
    keep = np.argsort(-np.abs(w))[: 9 - nnull]
    w = np.sort(np.array(w)[keep])
    rows.append(dict(level=-1, label='(150,770) unpruned', E180=e_ref,
                     bend_w=float(w[0]),
                     bend_cm=float(np.sign(w[0]) * CONV * np.sqrt(abs(w[0])))))
    say(f'unpruned (150,770): bend {rows[-1]["bend_cm"]:+.2f} cm-1')
    dump('grid', rows)


# ---------------------------------------------------------------- linear ----
def cmd_linear():
    """C_infinity_v: bond length only, then a 1D bend fit around 180."""
    r = 2.00
    rows = []
    for it in range(8):
        es = []
        for dr in (-0.02, 0.0, 0.02):
            e, g, mf = energy_grad(geom(180, r + dr))
            es.append((e, g))
        a, b, c = [z[0] for z in es]
        d1 = (c - a) / 0.04
        d2 = (a - 2 * b + c) / 0.02 ** 2
        rows.append(dict(it=it, r=r, E=b, dEdr=d1, d2dr2=d2))
        say(f'linear it {it}: r {r:.6f} E {b:.9f} dE/dr {d1:+.2e} k_r {d2:.5f}')
        if abs(d1) < 1e-6:
            break
        r = r - d1 / d2
    e_lin, g_lin, _ = energy_grad(geom(180, r))
    fit = []
    for th in [180, 178, 176, 174, 172, 170, 168, 165, 160]:
        e, _ = scf(geom(th, r))
        fit.append((th, e))
    th = np.array([f[0] for f in fit])
    ee = np.array([f[1] for f in fit])
    d = np.radians(180.0 - th)                       # bend angle from linear
    A2 = np.vstack([d ** 2, d ** 4]).T
    coef, *_ = np.linalg.lstsq(A2, ee - ee[0], rcond=None)
    # harmonic bend wavenumber of a 3-body bend from k_theta = 2*a2
    mF, mCa = 18.998403, 40.078
    G = 2.0 / (mF * (r / BOHR) ** 2) + 0.0        # 1/(m r^2) each, relative angle
    G = 2.0 / (mF * (r / BOHR) ** 2)
    lam = 2 * coef[0] * G                          # gaussian -> amu bohr^2
    cm = np.sign(lam) * CONV * np.sqrt(abs(lam))
    dump('linear', dict(r=r, E=e_lin, fit=[(float(a), float(b)) for a, b in fit],
                        a2=float(coef[0]), a4=float(coef[1]),
                        quartic_ratio=float(coef[1] * np.radians(20) ** 2 / coef[0]),
                        bend_cm_fit=float(cm), traj=rows))


if __name__ == '__main__':
    globals()[f'cmd_{sys.argv[1]}']()
