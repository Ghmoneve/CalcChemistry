"""Probe molecular capabilities for the CaF2 ALD model reaction.

Target reaction (gas phase, truncated ligand model):
    Ca(formamidinate)2 + 2 HF  ->  CaF2 + 2 formamidine

The real ligand is N,N'-di-i-propylformamidinate; the model replaces the two
i-propyl groups by H. This is a model reaction and is labelled as one.

This probe only establishes what the installation can do: which basis sets are
available for Ca, whether an all-electron or ECP description is used, and
whether a geometry optimiser is present.
"""
import numpy as np
from pyscf import gto, scf, dft

print("=== geometry optimiser available? ===")
try:
    from pyscf.geomopt.geometric_solver import optimize  # noqa: F401
    print("  geomeTRIC: yes")
except Exception as exc:
    print(f"  geomeTRIC: no ({type(exc).__name__})")

print("\n=== Ca basis / ECP behaviour ===")
for basis, ecp in (('sto-3g', None), ('def2-svp', None), ('def2-svp', 'def2-svp'),
                   ('lanl2dz', 'lanl2dz')):
    try:
        m = gto.M(atom='Ca 0 0 0', basis=basis, ecp=ecp or {},
                  verbose=0, spin=0, charge=0)
        print(f"  basis={basis:10s} ecp={str(ecp):10s} nelec={m.nelectron:3d} "
              f"nao={m.nao_nr():3d}")
    except Exception as exc:
        print(f"  basis={basis:10s} ecp={str(ecp):10s} FAIL {type(exc).__name__}")

print("\n=== test molecules build and one SCF each (PBE) ===")
MOLS = {
    'HF': ('H 0 0 0; F 0 0 0.917', 0, 0),
    'CaF2': ('Ca 0 0 0; F 0 0 2.10; F 0 0 -2.10', 0, 0),
    'formamidine': ('C 0 0 0; N 0 1.35 0; N 0 -1.35 0; '
                    'H 0 0 1.10; H 0.90 1.90 0; H -0.90 1.90 0; '
                    'H 0.90 -1.90 0; H -0.90 -1.90 0', 0, 0),
    'Ca_famd2': ('Ca 0 0 0; '
                 'N 1.30 1.30 1.30; N 2.60 0.00 1.30; '
                 'C 2.60 1.30 2.30; '
                 'N -1.30 -1.30 -1.30; N -2.60 0.00 -1.30; '
                 'C -2.60 -1.30 -2.30; '
                 'H 3.50 0.00 1.30; H -3.50 0.00 -1.30; '
                 'H 3.30 1.30 2.90; H 1.90 -1.30 -2.90', 0, 0),
}
for name, (geom, spin, chg) in MOLS.items():
    try:
        m = gto.M(atom=geom, basis='def2-svp', verbose=0, spin=spin, charge=chg)
        mf = dft.RKS(m)
        mf.xc = 'pbe'
        mf.max_cycle = 200
        e = mf.kernel()
        print(f"  {name:14s} nelec={m.nelectron:3d} nao={m.nao_nr():3d} "
              f"E={e:.6f} conv={mf.converged}")
    except Exception as exc:
        print(f"  {name:14s} FAIL {type(exc).__name__}: {str(exc)[:70]}")
