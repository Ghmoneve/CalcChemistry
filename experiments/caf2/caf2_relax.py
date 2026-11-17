"""Robust relaxation of a small molecule on a very flat bending potential.

This is the recipe distilled from experiments/caf2/flat_mode_study.py.  It fixes
the three things that make a damped BFGS with a fixed step cap stall on CaF2:

  1. the inverse Hessian starts from the analytic Hessian, not 0.3*I
     (0.3 is ~800x too stiff along the Ca-F-Ca bend, whose curvature is
     ~4e-4 Ha/(bohr^2 amu)),
  2. the BFGS update is only applied when the secant pair actually carries
     curvature (s.y > 1e-8 |s||y|); without that test rho = 1/(s.y) explodes at
     the convergence floor and B runs away,
  3. the gradient is taken with grid_response=True, so that it is the
     derivative of the energy that is being minimised (the atom-centred
     quadrature otherwise leaves a ~2e-4 Ha/Angstrom spurious net force), and
     convergence is judged on the translation/rotation-projected gradient.

Usage:
    python caf2_relax.py            # PBE/def2-svp, grid 5, from a 140 deg start
"""
from __future__ import annotations

import numpy as np
from pyscf import dft, gto
from pyscf.hessian import rks as hrks

BOHR = 0.52917721092
CONV = 5140.02


def build(coords, basis='def2-svp', zs=(20, 9, 9)):
    return gto.M(atom=[(zs[i], tuple(coords[i])) for i in range(len(zs))],
                 basis=basis, spin=0, charge=0, unit='Angstrom', verbose=0)


def energy_grad(mol, xc='pbe', grid=5, conv_tol=1e-11):
    mf = dft.RKS(mol)
    mf.xc = xc
    mf.grids.level = grid
    mf.conv_tol = conv_tol
    mf.kernel()
    assert mf.converged, 'SCF did not converge'
    gm = mf.nuc_grad_method()
    gm.grid_response = True          # derivative of the energy actually computed
    g = gm.kernel().reshape(-1) / BOHR
    return float(mf.e_tot), g, mf


def internal_projector(mol, coords):
    """Projector onto the 3N-6 (or 3N-5) internal space, mass-weighted."""
    sq = np.repeat(np.sqrt(np.array(mol.atom_mass_list())), 3)
    n = 3 * mol.natm
    cands = []
    for a in range(3):
        v = np.zeros(n)
        for i in range(mol.natm):
            v[3 * i + a] = sq[3 * i]
        cands.append(v)
    for a in range(3):
        e = np.zeros(3)
        e[a] = 1.0
        v = np.zeros(n)
        for i in range(mol.natm):
            v[3 * i:3 * i + 3] = sq[3 * i] * np.cross(e, coords[i])
        cands.append(v)
    basis = []
    for v in cands:
        nv = np.linalg.norm(v)
        if nv < 1e-6:
            continue
        v = v / nv
        for u in basis:
            v -= (v @ u) * u
        nv = np.linalg.norm(v)
        if nv > 1e-4:
            basis.append(v / nv)
    P = np.eye(n) - sum(np.outer(u, u) for u in basis)
    return P, sq


def hessian(mol, xc='pbe', grid=5):
    mf = dft.RKS(mol)
    mf.xc = xc
    mf.grids.level = grid
    mf.conv_tol = 1e-11
    mf.kernel()
    h = np.array(hrks.Hessian(mf).kernel())
    H = h.transpose(0, 2, 1, 3).reshape(3 * mol.natm, 3 * mol.natm)
    return 0.5 * (H + H.T), mf


def frequencies(mol, coords, xc='pbe', grid=5):
    H, mf = hessian(mol, xc, grid)
    P, sq = internal_projector(mf.mol, coords)
    Hm = P @ (H / np.outer(sq, sq)) @ P
    w, _ = np.linalg.eigh(0.5 * (Hm + Hm.T))
    keep = np.argsort(-np.abs(w))[: int(round(np.trace(P)))]
    w = np.sort(np.array(w)[keep])
    return np.sign(w) * CONV * np.sqrt(np.abs(w))


def relax(coords, basis='def2-svp', xc='pbe', grid=5, gtol=1e-5,
          maxiter=40, maxstep=0.15, log=print):
    x = np.array(coords, float)
    mol = build(x, basis)
    P, sq = internal_projector(mol, x)

    # ---- step 1: seed the inverse Hessian with the analytic Hessian --------
    H, _ = hessian(mol, xc, grid)
    Hp = P @ (H / np.outer(sq, sq)) @ P
    w, V = np.linalg.eigh(0.5 * (Hp + Hp.T))
    keep = np.argsort(-np.abs(w))[: int(round(np.trace(P)))]
    hinv = np.zeros_like(Hp)
    for k in keep:
        hinv += np.outer(V[:, k], V[:, k]) / w[k]

    gp = xp = None
    e = 0.0
    for it in range(maxiter):
        mol = build(x, basis)
        e, g, mf = energy_grad(mol, xc, grid)
        gi = float(np.linalg.norm(P @ g))
        log(f'it {it:3d}  E {e:.9f}  |g| {np.linalg.norm(g):.2e}  '
            f'|P g| {gi:.2e}  angle {angle(x):7.2f}')
        if gi < gtol:
            break
        if gp is not None:
            s = (x - xp).reshape(-1)
            y = g - gp
            sy = float(s @ y)
            # ---- step 2: relative curvature guard -------------------------
            if sy > 1e-8 * np.linalg.norm(s) * np.linalg.norm(y):
                rho = 1.0 / sy
                idn = np.eye(len(g))
                hinv = ((idn - rho * np.outer(s, y)) @ hinv
                        @ (idn - rho * np.outer(y, s)) + rho * np.outer(s, s))
        step = -hinv @ g
        ns = np.linalg.norm(step)
        if ns > maxstep:
            step *= maxstep / ns
        # ---- step 3: backtracking line search (BFGS alone guarantees nothing)
        for _ in range(8):
            xt = x + step.reshape(-1, 3)
            et, _, _ = energy_grad(build(xt, basis), xc, grid)
            if et < e - 1e-10:
                break
            step *= 0.5
        else:
            log('  line search stalled at the energy noise floor; stopping')
            break
        xp, gp = x.copy(), g.copy()
        x = xt
    return e, x, gi


def angle(x):
    v1, v2 = x[1] - x[0], x[2] - x[0]
    c = v1 @ v2 / np.linalg.norm(v1) / np.linalg.norm(v2)
    return float(np.degrees(np.arccos(np.clip(c, -1, 1))))


def bent_start(theta=140.0, r=1.9925):
    th = np.radians(theta)
    return np.array([[0.0, 0.0, 0.0],
                     [r * np.sin(th / 2), 0.0, r * np.cos(th / 2)],
                     [-r * np.sin(th / 2), 0.0, r * np.cos(th / 2)]])


if __name__ == '__main__':
    import sys
    basis = sys.argv[1] if len(sys.argv) > 1 else 'def2-svp'
    xc = sys.argv[2] if len(sys.argv) > 2 else 'pbe'
    grid = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    print(f'== {xc}/{basis}, grid level {grid} ==')
    e, x, gi = relax(bent_start(), basis=basis, xc=xc, grid=grid)
    print(f'\nE = {e:.9f} Ha   |P g| = {gi:.2e}   angle = {angle(x):.2f} deg')
    print('r(Ca-F) =', np.round(np.linalg.norm(x[1:] - x[0], axis=1), 5))
    print('frequencies cm^-1:', np.round(frequencies(build(x, basis), x, xc, grid), 1))
