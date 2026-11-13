import numpy as np
from pyscf import gto, scf, dft
from pyscf.pbc import gto as pgto

print("pyscf probe")
print("numpy", np.__version__)

# 1) HF molecule, all-electron, sanity check
mol = gto.M(atom='H 0 0 0; F 0 0 0.917', basis='def2-svp', verbose=0)
mf = scf.RHF(mol).run()
print(f"HF/def2-svp RHF energy = {mf.e_tot:.8f} Ha   (lit. ~ -100.02)")

# 2) Ca: which ECPs are available?
from pyscf.gto import basis
for name in ['lanl2dz','ccecp','def2-svp','sto-3g']:
    try:
        b = basis.load(name, 'Ca')
        print(f"Ca basis '{name}': loaded, n shells = {len(b)}")
    except Exception as e:
        print(f"Ca basis '{name}': FAIL {type(e).__name__}")
for ecp in ['lanl2dz','ccecp','def2-svp']:
    try:
        e = basis.load_ecp(ecp, 'Ca')
        print(f"Ca ECP '{ecp}': loaded, ncore = {e[0][0] if e else '?'}")
    except Exception as e:
        print(f"Ca ECP '{ecp}': FAIL {type(e).__name__}")

# 3) CaF2 bulk primitive cell, periodic, check lattice constant
a0 = 5.463  # Angstrom, experimental fluorite lattice constant
cell = pgto.Cell()
cell.atom = '''
Ca 0.0 0.0 0.0
Ca 0.0 2.7315 2.7315
Ca 2.7315 0.0 2.7315
Ca 2.7315 2.7315 0.0
F  1.36575 1.36575 1.36575
F  1.36575 3.09725 3.09725
F  3.09725 1.36575 3.09725
F  3.09725 3.09725 1.36575
'''
cell.a = np.eye(3)*a0
cell.basis = 'sto-3g'
cell.ecp = {'Ca':'lanl2dz'}
cell.verbose = 0
cell.build()
print("CaF2 cell built: natm", cell.natm, "nelectron", cell.nelectron)
