"""Published CaF2 ALD growth per cycle -> fractional surface coverage.

Every physical input is a literature value with its source. Nothing is fitted
and nothing is adjusted to reach a desired answer. The interval propagation is
proved in CalcChemistry/SurfaceReaction.lean (gpc_coverage_interval,
coverage_gpc_interval, max_density_from_ligand_size).

SOURCES
  [A] CaF2 fluorite lattice constant a = 5.4626 A, room-temperature XRD powder
      data: H. E. Swanson, E. Tatge, NBS Circular 539 (1953) 69-70; ICSD 52754.
  [B] CaF2(111) F-Ca-F trilayer height 3.15 A, quoted experimentally:
      Calleja et al., Surf. Sci. 582 (2005) 14-20, doi:10.1016/j.susc.2005.03.003.
      Independently derived here as a/sqrt(3) from [A].
  [C] Growth per cycle 0.4 A/cycle (175 C) and 0.3 A/cycle (250 C) for
      [Ca(amd)2]2 + anhydrous HF: J. Vac. Sci. Technol. A 41 (2023) 042405.
  [D] F-terminated, (1x1), no reconstruction: Vogt et al., Surf. Sci. 578 (2005)
      57-70, doi:10.1016/j.susc.2005.01.011; Puchin et al., J. Phys.: Condens.
      Matter 13 (2001) 2081, doi:10.1088/0953-8984/13/10/304.
"""
from fractions import Fraction as F
import math, json

A0 = F(54626, 10000)             # [A] 5.4626 Angstrom
TRILAYER = A0 / F(3).limit_denominator()  # placeholder, recomputed below
TRILAYER = A0 / F(math.isqrt(3) ** 0 + 0) if False else None
import decimal
decimal.getcontext().prec = 30
TRILAYER_F = float(A0) / math.sqrt(3.0)          # a/sqrt(3)
GPC_F = (0.3, 0.4)                               # [C]

cov_lo, cov_hi = GPC_F[0] / TRILAYER_F, GPC_F[1] / TRILAYER_F

# Ca site lattice on the F-terminated (111) plane (no reconstruction, [D]).
d_nn = float(A0) / math.sqrt(2.0)                # in-plane nearest neighbour
area_site = (math.sqrt(3.0) / 2.0) * d_nn ** 2   # area per Ca site
rho_site = 1.0 / area_site                       # sites per A^2
sigma_lo, sigma_hi = rho_site * cov_lo, rho_site * cov_hi

print("CaF2 ALD growth per cycle -> fractional layer coverage")
print(f"  lattice constant a0           = {float(A0):.4f} A            [A]")
print(f"  trilayer height a0/sqrt(3)    = {TRILAYER_F:.4f} A            [A],[B]")
print(f"  GPC (175-250 C)               = 0.3 - 0.4 A/cycle        [C]")
print(f"  fractional coverage           = {cov_lo*100:.2f}% - {cov_hi*100:.2f}% of one layer per cycle")
print()
print("Ca site geometry on CaF2(111), F-terminated, no reconstruction  [D]")
print(f"  in-plane Ca-Ca distance       = {d_nn:.4f} A   (a0/sqrt 2)")
print(f"  area per Ca site              = {area_site:.4f} A^2")
print(f"  Ca site density               = {rho_site:.5f} sites/A^2")
print(f"  Ca deposited per cycle        = {sigma_lo:.5f} - {sigma_hi:.5f} sites/A^2")
print(f"  mean area claimed per Ca       = {1/sigma_lo:.1f} - {1/sigma_hi:.1f} A^2")
print(f"  ratio to bare site area       = {(1/sigma_lo)/area_site:.1f}x - {(1/sigma_hi)/area_site:.1f}x")
print()
print("Design-rule reading (Lean: max_density_from_ligand_size, gpc_upper_from_ligand_footprint)")
print("  A ligand footprint of area A0 caps the cycle at sigma <= 1/A0.")
print(f"  Observed sigma {sigma_lo:.5f}-{sigma_hi:.5f} A^-2 therefore requires only A0 <= {1/sigma_lo:.1f} A^2.")
print("  => the observed GPC is CONSISTENT with a ligand footprint of order 100 A^2, but")
print("     consistency is not proof: kinetics could equally be the limiter.")

json.dump({'a0': float(A0), 'trilayer': TRILAYER_F, 'gpc': list(GPC_F),
           'coverage': [cov_lo, cov_hi], 'd_nn': d_nn, 'area_site': area_site,
           'rho_site': rho_site, 'sigma': [sigma_lo, sigma_hi],
           'area_per_metal': [1/sigma_lo, 1/sigma_hi]},
          open('experiments/caf2/gpc_coverage.json', 'w'), indent=2)
print("\nwritten experiments/caf2/gpc_coverage.json")
