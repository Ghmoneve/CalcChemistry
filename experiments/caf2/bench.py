"""Benchmark pyscf periodic DFT settings on bulk CaF2.

Purpose: find a setting that reaches a converged SCF in a tractable time, so the
full lattice scan can run with known cost per point. No physical claim is made
here; the lattice constant from the chosen setting is validated separately
against the experimental 5.4626 A.
"""
import time
import numpy as np
from pyscf.pbc import gto as pgto, dft as pdft

A0 = 5.4626


def build(a0, mesh):
    c = pgto.Cell()
    c.atom = f"""
Ca 0.0 0.0 0.0
Ca 0.0 {a0/2} {a0/2}
Ca {a0/2} 0.0 {a0/2}
Ca {a0/2} {a0/2} 0.0
F {a0/4} {a0/4} {a0/4}
F {a0/4} {3*a0/4} {3*a0/4}
F {3*a0/4} {a0/4} {3*a0/4}
F {3*a0/4} {3*a0/4} {a0/4}"""
    c.a = np.eye(3) * a0
    c.basis = 'sto-3g'
    c.ecp = {'Ca': 'lanl2dz'}
    c.verbose = 0
    c.mesh = [mesh] * 3
    c.build()
    return c


CASES = [
    ('mesh30_k1_df', 30, 1, True),
    ('mesh40_k2_df', 40, 2, True),
    ('mesh40_k2_nodf', 40, 2, False),
    ('mesh60_k4_df', 60, 4, True),
]

for tag, mesh, km, use_df in CASES:
    t0 = time.time()
    try:
        cell = build(A0, mesh)
        kpts = cell.make_kpts([km, km, km])
        mf = pdft.KRKS(cell, kpts)
        mf.xc = 'pbe'
        mf.max_cycle = 60
        if use_df:
            mf = mf.density_fit()
        e = mf.kernel()
        print(f"{tag:18s} nk={len(kpts):3d} E={e:.6f} conv={mf.converged} "
              f"{time.time()-t0:7.1f}s", flush=True)
    except Exception as exc:
        print(f"{tag:18s} FAIL {type(exc).__name__}: {str(exc)[:90]} "
              f"{time.time()-t0:7.1f}s", flush=True)
