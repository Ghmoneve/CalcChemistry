"""PBE lattice constant of bulk CaF2 (fluorite). Validation of the whole setup:
the computed a0 must be compared against experiment (5.463 A), never adjusted
to match it."""
import numpy as np, json, time, sys
from pyscf.pbc import gto as pgto, dft as pdft

MESH = [60, 60, 60]
KMESH = [4, 4, 4]
FUNC = 'pbe'
BASIS = 'sto-3g'
ECP = {'Ca': 'lanl2dz'}

def build(a0):
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
    cell.a = np.eye(3) * a0
    cell.basis = BASIS
    cell.ecp = ECP
    cell.verbose = 0
    cell.mesh = MESH
    cell.build()
    return cell

results = []
alist = [round(5.20 + 0.05*i, 3) for i in range(13)]
for a0 in alist:
    t0 = time.time()
    cell = build(a0)
    kpts = cell.make_kpts(KMESH)
    mf = pdft.KRKS(cell, kpts)
    mf.xc = FUNC
    mf.max_cycle = 100
    e = mf.kernel()
    n = cell.natm
    e_per = e / n
    results.append({'a0': a0, 'e_tot': float(e), 'e_per_atom': float(e_per),
                    'converged': bool(mf.converged), 'sec': round(time.time()-t0, 1)})
    json.dump(results, open('experiments/caf2/bulk_scan.json', 'w'), indent=2)
    print(f"a0={a0:.3f}  E={e:.8f} Ha  E/atom={e_per:.8f}  conv={mf.converged}  {time.time()-t0:.1f}s", flush=True)

# parabolic fit around the minimum
best = min(results, key=lambda r: r['e_per_atom'])
i = results.index(best)
if 0 < i < len(results)-1:
    ys = [results[i-1]['e_per_atom'], best['e_per_atom'], results[i+1]['e_per_atom']]
    xs = [results[i-1]['a0'], best['a0'], results[i+1]['a0']]
    c = np.polyfit(xs, ys, 2)
    a_fit = -c[1]/(2*c[0])
    print(f"\nparabolic minimum: a0 = {a_fit:.4f} A   (experiment 5.463 A)")
    json.dump({'a0_fit': float(a_fit), 'experiment': 5.463, 'results': results},
              open('experiments/caf2/bulk_scan.json', 'w'), indent=2)
