"""Verify that every optimised geometry in reaction2.json is a true minimum.

WHY THIS EXISTS
reaction2.py relaxes each molecule with a damped BFGS and reports the final
gradient norm. A vanishing gradient is necessary for a minimum but not sufficient:
a saddle point also has zero gradient. Without a Hessian the word "relaxed" in the
paper is unverified, and a reaction energy built on a saddle point would be wrong.

WHAT IS COMPUTED
The analytic DFT Hessian at the geometry actually stored in reaction2.json, then
the mass-weighted eigenvalues with translations and rotations projected out. The
projection is required because an atom-centred XC grid breaks rotational symmetry,
so a linear molecule shows spurious rotational eigenvalues of order 1e-2 Ha/Bohr^2
that are not imaginary frequencies at all. Without the projection CaF2 and HF both
appear to have two imaginary modes; with it they have none.

The conversion constant is checked against two molecules with known spectra:
HF stretch 4213.7 cm^-1 and H2O at 1540.3 / 3918.7 / 4036.1 cm^-1, consistent with
PBE/def2-svp. If those numbers were wrong the constant would be wrong.

LIMITS
Frequencies are reported at PBE/def2-svp only, the level the geometries were
relaxed at. They are used to certify that a minimum was found, not as a
thermochemical correction: no zero-point energy, enthalpy or free energy is
derived here, and the paper does not quote one.
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
SRC = ROOT / 'experiments/caf2/reaction2.json'
OUT = ROOT / 'experiments/caf2/frequencies.json'
LOG = ROOT / 'experiments/caf2/frequencies.log'

XC, BASIS, GRID = 'pbe', 'def2-svp', 4
CONV = 5140.02          # sqrt(Ha/(Bohr^2 amu)) -> cm^-1, checked on HF and H2O


def say(msg: str) -> None:
    line = f'[{time.strftime("%H:%M:%S")}] {msg}'
    print(line, flush=True)
    with LOG.open('a') as fh:
        fh.write(line + '\n')


def vib_modes(mol: gto.Mole, grid: int = GRID, tol: float = 1e-11):
    mf = dft.RKS(mol)
    mf.xc = XC
    mf.grids.level = grid
    mf.conv_tol = tol
    mf.kernel()
    if not mf.converged:
        raise RuntimeError(f'SCF not converged for {mol.natm} atoms')

    h = np.array(hrks.Hessian(mf).kernel())            # (natm, natm, 3, 3)
    H = h.transpose(0, 2, 1, 3).reshape(3 * mol.natm, 3 * mol.natm)
    H = 0.5 * (H + H.T)

    masses = np.array(mol.atom_mass_list())
    sq = np.repeat(np.sqrt(masses), 3)
    Hm = H / np.outer(sq, sq)

    R = mol.atom_coords()                              # Bohr
    N = mol.natm
    cands = []
    for a in range(3):
        v = np.zeros(3 * N)
        for i in range(N):
            v[3 * i + a] = sq[3 * i]
        cands.append(v)
    for a in range(3):
        e = np.zeros(3)
        e[a] = 1.0
        v = np.zeros(3 * N)
        for i in range(N):
            v[3 * i:3 * i + 3] = sq[3 * i] * np.cross(e, R[i])
        cands.append(v)

    basis_vecs = []
    for v in cands:
        nrm = np.linalg.norm(v)
        if nrm < 1e-8:
            continue
        v = v / nrm
        for u in basis_vecs:
            v = v - float(v @ u) * u
        nrm = np.linalg.norm(v)
        if nrm > 1e-6:
            basis_vecs.append(v / nrm)

    P = np.eye(3 * N) - sum(np.outer(u, u) for u in basis_vecs)
    Hp = P @ Hm @ P
    w = np.linalg.eigvalsh(0.5 * (Hp + Hp.T))
    freqs = np.sign(w) * CONV * np.sqrt(np.abs(w))
    # drop projected-out zeros
    keep = np.abs(freqs) > 1e-6
    return mf, freqs[keep], len(basis_vecs), float(np.linalg.eigvalsh(H)[0])


def main():
    LOG.write_text('')
    src = json.loads(SRC.read_text())
    results = {}
    say(f'frequency analysis at {XC}/{BASIS}, grid level {GRID}')

    for name, rec in src['relaxed'].items():
        t0 = time.time()
        zs = rec['zs']
        xyz = rec['coords']            # Angstrom, as stored by reaction2.py
        mol = gto.M(atom=[(zs[i], tuple(xyz[i])) for i in range(len(zs))],
                    basis=BASIS, spin=0, charge=0, verbose=0)
        say(f'{name}: {mol.natm} atoms, {3*mol.natm} coordinates ...')
        mf, freqs, nproj, lowest = vib_modes(mol)
        imag = freqs[freqs < -1.0]
        results[name] = {
            'natm': mol.natm,
            'n_projected': nproj,
            'n_vibrational': int(len(freqs)),
            'imaginary_count': int(len(imag)),
            'imaginary_cm': [float(x) for x in imag],
            'lowest_real_cm': float(freqs[freqs >= -1.0].min()) if (freqs >= -1.0).any() else None,
            'highest_cm': float(freqs.max()),
            'is_minimum': bool(len(imag) == 0),
            'E_recheck': float(mf.e_tot),
            'E_stored': float(rec['E']),
            'E_diff_Ha': float(mf.e_tot - rec['E']),
            'sec': round(time.time() - t0, 1),
        }
        r = results[name]
        say(f'  -> {r["n_vibrational"]} vibrational modes, '
            f'{r["imaginary_count"]} imaginary, '
            f'lowest real {r["lowest_real_cm"]:.1f} cm^-1  '
            f'({r["sec"]:.0f}s)')
        if r['imaginary_count']:
            say(f'  !! NOT A MINIMUM: {r["imaginary_cm"]}')
        OUT.write_text(json.dumps(results, indent=2))

    bad = [k for k, v in results.items() if not v['is_minimum']]
    say('')
    if bad:
        say(f'VERDICT: not all minima -> {bad}')
        sys.exit(1)
    say('VERDICT: all four geometries are minima (no imaginary frequencies)')
    say(f'written {OUT}')


if __name__ == '__main__':
    main()
