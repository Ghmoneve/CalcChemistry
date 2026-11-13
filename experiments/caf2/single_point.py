import numpy as np, time, sys
from pyscf.pbc import gto as pgto, dft as pdft

a0 = 5.463
def build(a0, mesh_pts=40):
    cell = pgto.Cell()
    cell.atom = f'''
Ca 0.0 0.0 0.0
Ca 0.0 {a0/2} {a0/2}
Ca {a0/2} 0.0 {a0/2}
Ca {a0/2} {a0/2} 0.0
F  {a0/4} {a0/4} {a0/4}
F  {a0/4} {3*a0/4} {3*a0/4}
F  {3*a0/4} {a0/4} {3*a0/4}
F  {3*a0/4} {3*a0/4} {a0/4}
'''
    cell.a = np.eye(3)*a0
    cell.basis = 'sto-3g'
    cell.ecp = {'Ca':'lanl2dz'}
    cell.verbose = 0
    cell.mesh = [mesh_pts]*3          # explicit, prevents a runaway FFT grid
    cell.build()
    return cell

t0=time.time()
cell = build(a0)
print("FFT mesh =", cell.mesh, " natm =", cell.natm, " nelec =", cell.nelectron)
kpts = cell.make_kpts([4,4,4])
mf = pdft.KRKS(cell, kpts)
mf.xc = 'pbe'
mf.max_cycle = 80
e = mf.kernel()
print(f"a0={a0}  PBE/sto-3g  E = {e:.8f} Ha   nkpts={len(kpts)}  converged={mf.converged}  t={time.time()-t0:.1f}s")
