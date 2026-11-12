"""
Yangian Y(sl(2)) Symmetry Verification for Polyacetylene Model
===============================================================

Tests Conjecture 4.3: Does the polyacetylene tight-binding Hamiltonian
possess Yangian Y(sl(2)) symmetry?

Model: H = -t Σ_j (c†_j c_{j+1} + h.c.) + U Σ_j n_{j↑} n_{j↓}
"""

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import kron, eye, csr_matrix
from itertools import product
import matplotlib.pyplot as plt

# Pauli matrices for spin-1/2
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
sigma_p = (sigma_x + 1j*sigma_y) / 2  # S^+
sigma_m = (sigma_x - 1j*sigma_y) / 2  # S^-

def fermion_operators(L):
    """
    Build fermionic creation/annihilation operators for L sites, 2 spins.
    Returns: c_dagger[site, spin], c[site, spin] as sparse matrices
    """
    # Total Hilbert space dimension: 2^(2L)
    dim = 2**(2*L)
    
    c_dag = {}
    c = {}
    
    for site in range(L):
        for spin in [0, 1]:  # 0=up, 1=down
            idx = 2*site + spin
            
            # Build Jordan-Wigner string
            phase = eye(2, format='csr')
            for k in range(idx):
                phase = kron(phase, -sigma_z, format='csr')
            
            # Apply sigma_- at site idx
            op = phase
            for k in range(idx):
                op = kron(eye(2**(k), format='csr'), op, format='csr') if k == 0 else op
            op = kron(op, csr_matrix(sigma_m), format='csr')
            for k in range(idx+1, 2*L):
                op = kron(op, eye(2, format='csr'), format='csr')
            
            c[(site, spin)] = op
            c_dag[(site, spin)] = op.conj().T
    
    return c_dag, c

def polyacetylene_hamiltonian(L, t=1.0, U=0.0, periodic=True):
    """
    Construct polyacetylene Hamiltonian in tight-binding approximation.
    
    H = -t Σ_j,σ (c†_{j,σ} c_{j+1,σ} + h.c.) + U Σ_j n_{j,↑} n_{j,↓}
    
    Parameters:
    - L: number of sites
    - t: hopping parameter
    - U: on-site Coulomb repulsion
    - periodic: use periodic boundary conditions
    """
    dim = 2**(2*L)
    H = csr_matrix((dim, dim), dtype=complex)
    
    c_dag, c = fermion_operators(L)
    
    # Hopping term
    max_j = L if periodic else L-1
    for j in range(max_j):
        j_next = (j + 1) % L
        for spin in [0, 1]:
            H += -t * (c_dag[(j, spin)] @ c[(j_next, spin)] + 
                       c[(j_next, spin)].conj().T @ c_dag[(j, spin)].conj().T)
    
    # Hubbard U term
    if U != 0:
        for j in range(L):
            n_up = c_dag[(j, 0)] @ c[(j, 0)]
            n_down = c_dag[(j, 1)] @ c[(j, 1)]
            H += U * (n_up @ n_down)
    
    return H

def yangian_sl2_generators_level0(L):
    """
    Level-0 Yangian generators: J^a = Σ_j S^a_j
    For sl(2): J^+, J^-, J^z (total spin operators)
    """
    dim = 2**(2*L)
    c_dag, c = fermion_operators(L)
    
    # S^z_j = (n_{j,↑} - n_{j,↓})/2
    J_z = csr_matrix((dim, dim), dtype=complex)
    for j in range(L):
        n_up = c_dag[(j, 0)] @ c[(j, 0)]
        n_down = c_dag[(j, 1)] @ c[(j, 1)]
        J_z += (n_up - n_down) / 2
    
    # S^+_j = c†_{j,↑} c_{j,↓}
    J_plus = csr_matrix((dim, dim), dtype=complex)
    for j in range(L):
        J_plus += c_dag[(j, 0)] @ c[(j, 1)]
    
    # S^-_j = c†_{j,↓} c_{j,↑}
    J_minus = csr_matrix((dim, dim), dtype=complex)
    for j in range(L):
        J_minus += c_dag[(j, 1)] @ c[(j, 0)]
    
    return J_plus, J_minus, J_z

def yangian_sl2_generators_level1(L):
    """
    Level-1 Yangian generators: J^a_(1) 
    
    For spin chains, these involve:
    J^z_(1) = Σ_{j<k} (S^z_j [h_k] - [h_j] S^z_k)
    
    where h_j is the local Hamiltonian term.
    
    Simplified form for nearest-neighbor:
    J^z_(1) ~ Σ_j j·S^z_j
    """
    dim = 2**(2*L)
    c_dag, c = fermion_operators(L)
    
    # Weighted sum: J^z_(1) = Σ_j j·S^z_j
    J_z1 = csr_matrix((dim, dim), dtype=complex)
    for j in range(L):
        n_up = c_dag[(j, 0)] @ c[(j, 0)]
        n_down = c_dag[(j, 1)] @ c[(j, 1)]
        J_z1 += j * (n_up - n_down) / 2
    
    # J^+_(1) = Σ_j j·S^+_j
    J_plus1 = csr_matrix((dim, dim), dtype=complex)
    for j in range(L):
        J_plus1 += j * c_dag[(j, 0)] @ c[(j, 1)]
    
    # J^-_(1) = Σ_j j·S^-_j
    J_minus1 = csr_matrix((dim, dim), dtype=complex)
    for j in range(L):
        J_minus1 += j * c_dag[(j, 1)] @ c[(j, 0)]
    
    return J_plus1, J_minus1, J_z1

def commutator(A, B):
    """Compute [A, B] = AB - BA"""
    return A @ B - B @ A

def verify_yangian_symmetry(L, t=1.0, U=0.0, periodic=True, tolerance=1e-10):
    """
    Verify [H, Q] = 0 for Yangian generators Q.
    
    Returns:
    - symmetry_holds: bool
    - commutator_norms: dict of ||[H, Q]|| values
    """
    print(f"\n{'='*70}")
    print(f"Yangian Y(sl(2)) Symmetry Check: L={L}, t={t}, U={U}, periodic={periodic}")
    print(f"{'='*70}\n")
    
    H = polyacetylene_hamiltonian(L, t, U, periodic)
    
    # Level-0 generators
    J_p0, J_m0, J_z0 = yangian_sl2_generators_level0(L)
    
    # Level-1 generators
    J_p1, J_m1, J_z1 = yangian_sl2_generators_level1(L)
    
    results = {}
    
    # Check level-0 (standard su(2) symmetry)
    print("Level-0 generators (standard SU(2)):")
    for name, J in [("J^+", J_p0), ("J^-", J_m0), ("J^z", J_z0)]:
        comm = commutator(H, J)
        norm = np.linalg.norm(comm.toarray())
        results[f"level0_{name}"] = norm
        status = "✓ PASS" if norm < tolerance else "✗ FAIL"
        print(f"  ||[H, {name:3s}]|| = {norm:.2e}  {status}")
    
    # Check level-1 (Yangian extension)
    print("\nLevel-1 generators (Yangian Y(sl(2))):")
    for name, J in [("J^+_(1)", J_p1), ("J^-_(1)", J_m1), ("J^z_(1)", J_z1)]:
        comm = commutator(H, J)
        norm = np.linalg.norm(comm.toarray())
        results[f"level1_{name}"] = norm
        status = "✓ PASS" if norm < tolerance else "✗ FAIL"
        print(f"  ||[H, {name:7s}]|| = {norm:.2e}  {status}")
    
    # Overall verdict
    level0_ok = all(results[k] < tolerance for k in results if 'level0' in k)
    level1_ok = all(results[k] < tolerance for k in results if 'level1' in k)
    
    print(f"\n{'='*70}")
    if level0_ok and level1_ok:
        print("RESULT: ✓ Yangian Y(sl(2)) symmetry VERIFIED")
    elif level0_ok and not level1_ok:
        print("RESULT: ✗ SU(2) symmetry holds, but Yangian extension FAILS")
        print("        (Level-1 generators do not commute with H)")
    else:
        print("RESULT: ✗ Neither SU(2) nor Yangian symmetry holds")
    print(f"{'='*70}\n")
    
    return (level0_ok and level1_ok), results

def exact_spectrum_bethe_ansatz(L, t=1.0, U=0.0, N_electrons=None):
    """
    Attempt exact spectrum via Bethe ansatz (if Yangian symmetry holds).
    
    For XXZ/Hubbard chain, Bethe equations:
    exp(ik_α L) = ∏_{β≠α} S(k_α - k_β)
    
    where S is the scattering matrix determined by U/t.
    """
    print(f"\nBethe Ansatz Exact Spectrum (L={L}, N={N_electrons}):")
    print("="*70)
    
    if N_electrons is None:
        N_electrons = L  # Half-filling
    
    # For U=0 (free fermions), exact solution is trivial
    if abs(U) < 1e-10:
        print("U=0: Free fermion case")
        k_values = 2*np.pi * np.arange(L) / L
        epsilon_k = -2*t * np.cos(k_values)
        
        # Fill lowest N_electrons states
        epsilon_sorted = np.sort(epsilon_k)
        E_ground = 2 * np.sum(epsilon_sorted[:N_electrons//2])  # Factor 2 for spin
        
        print(f"  Ground state energy: E_0 = {E_ground:.6f}")
        print(f"  Single-particle energies: {epsilon_sorted[:5]} ...")
        
        return E_ground, epsilon_sorted
    
    else:
        print(f"U={U}: Non-trivial interaction")
        print("  Full Bethe ansatz solution requires solving transcendental equations.")
        print("  (Not implemented here; see Richardson-Gaudin or ABACUS codes)")
        return None, None

def numerical_diagonalization(L, t=1.0, U=0.0, periodic=True, k_lowest=10):
    """
    Direct numerical diagonalization for comparison.
    """
    print(f"\nNumerical Exact Diagonalization (L={L}):")
    print("="*70)
    
    H = polyacetylene_hamiltonian(L, t, U, periodic)
    
    print(f"  Hilbert space dimension: {H.shape[0]}")
    print(f"  Matrix size: {H.shape[0]}×{H.shape[1]}")
    print(f"  Computing {k_lowest} lowest eigenvalues...")
    
    # For small systems, use dense diagonalization
    if H.shape[0] <= 256:
        eigvals, eigvecs = eigh(H.toarray())
        eigvals_lowest = eigvals[:k_lowest]
    else:
        # For larger systems, use sparse methods
        from scipy.sparse.linalg import eigsh
        eigvals_lowest, eigvecs = eigsh(H, k=k_lowest, which='SA')
        eigvals_lowest = np.sort(eigvals_lowest)
    
    print(f"\n  Lowest {k_lowest} eigenvalues:")
    for i, E in enumerate(eigvals_lowest):
        print(f"    E_{i} = {E:12.8f}")
    
    return eigvals_lowest

def analyze_failure_conditions(L_range, t=1.0, U_range=[0.0, 1.0, 2.0, 4.0]):
    """
    Scan parameter space to find where Yangian symmetry breaks.
    """
    print("\n" + "="*70)
    print("FAILURE CONDITION ANALYSIS")
    print("="*70 + "\n")
    
    results = []
    
    for L in L_range:
        for U in U_range:
            symmetry_holds, norms = verify_yangian_symmetry(L, t, U, periodic=True, tolerance=1e-8)
            
            max_level1_norm = max(norms[k] for k in norms if 'level1' in k)
            
            results.append({
                'L': L,
                'U': U,
                't': t,
                'symmetry': symmetry_holds,
                'max_violation': max_level1_norm
            })
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY: Yangian Symmetry Verification Across Parameters")
    print("="*70)
    print(f"{'L':>3s} {'U/t':>6s} {'Symmetry':>10s} {'Max ||[H,J_(1)]||':>20s}")
    print("-"*70)
    
    for r in results:
        status = "PASS" if r['symmetry'] else "FAIL"
        print(f"{r['L']:3d} {r['U']/r['t']:6.2f} {status:>10s} {r['max_violation']:20.2e}")
    
    # Identify failure regime
    passing = [r for r in results if r['symmetry']]
    failing = [r for r in results if not r['symmetry']]
    
    print("\n" + "="*70)
    if len(failing) == 0:
        print("CONCLUSION: Yangian Y(sl(2)) symmetry holds for all tested parameters")
    else:
        print("CONCLUSION: Yangian symmetry FAILS under certain conditions")
        print("\nFailure conditions:")
        for r in failing:
            print(f"  - L={r['L']}, U/t={r['U']/r['t']:.2f}")
        
        if all(r['U'] > 0 for r in failing):
            print("\nPattern: Symmetry breaks when U > 0 (interaction present)")
    print("="*70 + "\n")
    
    return results

def main():
    """Main analysis pipeline"""
    
    print("\n" + "="*70)
    print(" YANGIAN Y(SL(2)) SYMMETRY IN POLYACETYLENE MODEL")
    print(" Verification of Conjecture 4.3")
    print("="*70 + "\n")
    
    # Test 1: Small chain, no interaction (should have Yangian symmetry)
    print("\nTEST 1: L=3, U=0 (non-interacting)")
    verify_yangian_symmetry(L=3, t=1.0, U=0.0, periodic=True)
    numerical_diagonalization(L=3, t=1.0, U=0.0, periodic=True, k_lowest=6)
    
    # Test 2: Small chain, weak interaction
    print("\nTEST 2: L=3, U=1.0 (weak Hubbard interaction)")
    verify_yangian_symmetry(L=3, t=1.0, U=1.0, periodic=True)
    numerical_diagonalization(L=3, t=1.0, U=1.0, periodic=True, k_lowest=6)
    
    # Test 3: Parameter scan
    print("\n\n" + "="*70)
    print(" PARAMETER SPACE SCAN")
    print("="*70 + "\n")
    
    results = analyze_failure_conditions(
        L_range=[2, 3, 4],
        t=1.0,
        U_range=[0.0, 0.5, 1.0, 2.0]
    )
    
    # Final verdict
    print("\n\n" + "="*70)
    print(" FINAL VERDICT ON CONJECTURE 4.3")
    print("="*70 + "\n")
    
    all_pass = all(r['symmetry'] for r in results)
    
    if all_pass:
        print("✓ CONJECTURE VERIFIED:")
        print("  Polyacetylene tight-binding Hamiltonian possesses")
        print("  Yangian Y(sl(2)) symmetry under tested conditions.")
        print("\n  Implications:")
        print("  - System is integrable via Bethe ansatz")
        print("  - Exact energy spectrum can be derived algebraically")
        print("  - Richardson-Gaudin formalism applies")
    else:
        print("✗ CONJECTURE REFUTED:")
        print("  Yangian Y(sl(2)) symmetry does NOT hold in general.")
        print("\n  Failure conditions:")
        
        failing = [r for r in results if not r['symmetry']]
        if failing:
            U_fail = [r['U'] for r in failing]
            if all(u > 0 for u in U_fail):
                print("  - Symmetry breaks when U ≠ 0 (electron-electron interaction)")
                print("  - Pure tight-binding (U=0) may preserve symmetry")
            
        print("\n  Minimal correction:")
        print("  - Replace simple position-weighted generator J_(1) = Σ_j j·S^z_j")
        print("    with proper Yangian generators involving local Hamiltonian:")
        print("    J^a_(1) = Σ_{j<k} [S^a_j, h_k] where h_k is local term")
        print("  - Or restrict to U=0 (non-interacting) regime")

if __name__ == "__main__":
    main()
