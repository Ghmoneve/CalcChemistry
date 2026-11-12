# Two-level reduced pairing toy model (historical h2 names retained).
# STATUS: UNTESTED — Julia was unavailable during this correction.
# No molecular H2 accuracy claim. H_N = H_z + (eps1 + eps2) I.
# Requires Combinatorics. BigFloat checks below use scalar algebra only;
# no generic BigFloat matrix eigensolver or interval arithmetic is claimed.

using LinearAlgebra
using Printf
using Combinatorics

# ============================================================================
# 1. RICHARDSON EQUATION SOLVER (ARBITRARY PRECISION)
# ============================================================================

function solve_richardson_equation_h2(eps1::T, eps2::T, g::T) where T <: Real
    # Return paired H_N energies, not the full-sector excitation gap.
    all(isfinite, (eps1, eps2, g)) && g >= 0 ||
        throw(ArgumentError("finite real parameters and g >= 0 required"))
    radius = hypot(eps1 - eps2, g)
    center = eps1 + eps2 - g
    return center - radius, center + radius
end


function verify_richardson_equation(E::T, eps1::T, eps2::T, g::T) where T <: Real
    # H_N secular residual; nothing means inapplicable, not verified.
    isfinite(E) || throw(ArgumentError("E_pair must be finite"))
    if iszero(g) || iszero(2*eps1 - E) || iszero(2*eps2 - E)
        return nothing
    end
    lhs = one(T)
    rhs = g * (1 / (2*eps1 - E) + 1 / (2*eps2 - E))
    return abs(lhs - rhs)
end


# ============================================================================
# 2. BETHE ANSATZ WAVEFUNCTION
# ============================================================================

function construct_bethe_wavefunction(E::T, eps1::T, eps2::T, g::T;
                                     branch::Int=1) where T <: Real
    # Stable paired vector for H_N energy E; branch=1/2, including poles.
    roots = solve_richardson_equation_h2(eps1, eps2, g)
    branch in (1, 2) || throw(ArgumentError("branch must be 1 or 2"))
    scale = max(one(T), abs(E), abs(eps1-eps2), abs(g))
    isfinite(E) && abs(E-roots[branch]) <= 64*eps(float(one(T)))*scale ||
        throw(ArgumentError("E must match the selected H_N paired root"))
    d = eps1 - eps2
    if iszero(g)
        ground = d <= 0 ? [one(T), zero(T)] : [zero(T), one(T)]
    else
        radius = hypot(d, g)
        ratio = (g/radius) / (1 + abs(d)/radius)
        ground = d <= 0 ? [one(T), ratio] : [ratio, one(T)]
        ground /= norm(ground)
    end
    return branch == 1 ? ground : [-ground[2], ground[1]]
end


# ============================================================================
# 3. FULL N=2 TOY HAMILTONIAN CONSTRUCTION (H_z)
# ============================================================================

function build_fock_basis(n_orb::Int, n_elec::Int)
    # Build Fock space basis.
    n_spin_orb = 2 * n_orb
    basis = []
    
    for occupied in combinations(1:n_spin_orb, n_elec)
        state = zeros(Int, n_spin_orb)
        for i in occupied
            state[i] = 1
        end
        push!(basis, Tuple(state))
    end
    
    basis_dict = Dict(state => i for (i, state) in enumerate(basis))
    return basis, basis_dict
end


function hamiltonian_pairing_fci(eps1::T, eps2::T, g::T, 
                                  basis::Vector, basis_dict::Dict) where T <: Real
    # Build the toy H_z matrix in Fock basis.
    n_states = length(basis)
    H = zeros(T, n_states, n_states)
    
    eps = [eps1, eps2]
    
    for (i, state) in enumerate(basis)
        state_arr = collect(state)
        
        # Diagonal: S^z terms
        for p in 1:2
            n_up = state_arr[2*p - 1]
            n_down = state_arr[2*p]
            H[i, i] += 2 * eps[p] * T(0.5) * (n_up + n_down - 1)
        end
        
        # Pairing terms INCLUDING p=q: annihilate before testing creation.
        for p in 1:2
            for q in 1:2
                # Check if can apply S^-_q (pair in q)
                if state_arr[2*q - 1] == 1 && state_arr[2*q] == 1
                    new_state = copy(state_arr)
                    new_state[2*q - 1] = 0
                    new_state[2*q] = 0
                    if new_state[2*p - 1] == 0 && new_state[2*p] == 0
                        # Add pair to p
                        new_state[2*p - 1] = 1
                        new_state[2*p] = 1
                        
                        j = basis_dict[Tuple(new_state)]
                        H[j, i] += -g
                    end
                end
            end
        end
    end
    
    return Hermitian(H)
end


# ============================================================================
# 4. HIGH-PRECISION VERIFICATION
# ============================================================================

function compare_bethe_vs_fci(eps1::T, eps2::T, g::T; verbose=true) where T <: Real
    # Float64 toy-sector comparison. Source corrected; runtime untested.
    
    # Bethe ansatz solution
    E_pair_ground, E_pair_excited = solve_richardson_equation_h2(eps1, eps2, g)
    E_bethe_ground = E_pair_ground - (eps1 + eps2)
    E_bethe_excited = E_pair_excited - (eps1 + eps2)
    richardson_error = verify_richardson_equation(E_pair_ground, eps1, eps2, g)
    coeffs_bethe = construct_bethe_wavefunction(E_pair_ground, eps1, eps2, g)
    
    # FCI solution
    basis, basis_dict = build_fock_basis(2, 2)
    H_fci = hamiltonian_pairing_fci(eps1, eps2, g, basis, basis_dict)
    
    eigenvalues_fci, eigenvectors_fci = eigen(H_fci)
    
    E_fci_ground = eigenvalues_fci[1]
    wf_fci = eigenvectors_fci[:, 1]
    
    # Energy error
    energy_error = abs(E_bethe_ground - E_fci_ground)
    
    # Wavefunction comparison
    wf_bethe_full = zeros(T, length(basis))
    for (i, state) in enumerate(basis)
        if state == (1, 1, 0, 0)  # |11̄⟩
            wf_bethe_full[i] = coeffs_bethe[1]
        elseif state == (0, 0, 1, 1)  # |22̄⟩
            wf_bethe_full[i] = coeffs_bethe[2]
        end
    end
    
    # Project onto the full numerical ground eigenspace, including degeneracy.
    scale = max(one(T), opnorm(H_fci))
    ground_mask = abs.(eigenvalues_fci .- E_fci_ground) .<= 1e-12*scale
    ground_space = eigenvectors_fci[:, ground_mask]
    wf_error = norm(wf_bethe_full - ground_space*(ground_space'*wf_bethe_full))
    expected_spectrum = sort([E_bethe_ground, E_bethe_excited, zeros(T, 4)...])
    spectrum_error = maximum(abs.(expected_spectrum - eigenvalues_fci))
    residual = norm(H_fci*wf_bethe_full - E_bethe_ground*wf_bethe_full)
    normalization_error = abs(norm(wf_bethe_full)-1)
    upper_coeffs = construct_bethe_wavefunction(E_pair_excited, eps1, eps2, g; branch=2)
    upper_full = zeros(T, length(basis))
    upper_full[basis_dict[(1, 1, 0, 0)]] = upper_coeffs[1]
    upper_full[basis_dict[(0, 0, 1, 1)]] = upper_coeffs[2]
    upper_residual = norm(H_fci*upper_full - E_bethe_excited*upper_full)
    upper_norm_error = abs(norm(upper_full)-1)
    passed = all(isfinite, (energy_error, spectrum_error, residual, upper_residual,
                            wf_error, normalization_error, upper_norm_error)) &&
             max(energy_error, spectrum_error, residual, upper_residual) <= 1e-12*scale &&
             max(wf_error, normalization_error, upper_norm_error) <= 1e-12
    
    if verbose
        println("=" ^ 70)
        println("TWO-LEVEL TOY PAIRING CHECK (Julia; H_z energies)")
        println("=" ^ 70)
        @printf("Parameters: ε₁=%.4f, ε₂=%.4f, g=%.4f\n", eps1, eps2, g)
        println()
        println("BETHE ANSATZ SOLUTION:")
        @printf("  Ground state energy: E = %.18f\n", E_bethe_ground)
        println("  H_N rational residual (nothing = inapplicable): ", richardson_error)
        @printf("  Wavefunction: |Ψ⟩ = %.8f|11̄⟩ + %.8f|22̄⟩\n", coeffs_bethe[1], coeffs_bethe[2])
        println()
        println("FULL SIX-STATE TOY DIAGONALIZATION:")
        @printf("  Ground state energy: E = %.18f\n", E_fci_ground)
        println()
        println("COMPARISON:")
        @printf("  |E_Bethe - E_FCI| = %.2e model energy units\n", energy_error)
        @printf("  ||Ψ_Bethe - Ψ_FCI|| = %.2e\n", wf_error)
        println()
        if passed
            println("PASS: finite-spectrum/vector checks within tolerance")
        else
            println("✗ DISCREPANCY DETECTED")
        end
        println("=" ^ 70)
    end
    
    return Dict(
        "E_bethe" => E_bethe_ground,
        "E_pair" => E_pair_ground,
        "E_fci" => E_fci_ground,
        "energy_error" => energy_error,
        "wf_error" => wf_error,
        "richardson_error" => richardson_error,
        "spectrum_error" => spectrum_error,
        "eigenvector_residual" => max(residual, upper_residual),
        "normalization_error" => max(normalization_error, upper_norm_error),
        "passed" => passed
    )
end


# ============================================================================
# 5. EXTENDED PRECISION TEST
# ============================================================================

function test_arbitrary_precision()
    # Optional 256-bit scalar polynomial/vector checks, not matrix FCI.
    return setprecision(BigFloat, 256) do
        eps1, eps2, g = BigFloat("-1.0"), BigFloat("-0.5"), BigFloat("0.3")
        roots = solve_richardson_equation_h2(eps1, eps2, g)
        for (branch, E) in enumerate(roots)
            coeffs = construct_bethe_wavefunction(E, eps1, eps2, g; branch=branch)
            d1, d2 = 2*eps1-E, 2*eps2-E
            polynomial_error = abs(d1*d2-g*(d1+d2))
            vector_error = hypot(d1*coeffs[1]-g*sum(coeffs),
                                 d2*coeffs[2]-g*sum(coeffs))
            @assert polynomial_error < BigFloat("1e-60")
            @assert vector_error < BigFloat("1e-60")
            @assert abs(norm(coeffs)-1) < BigFloat("1e-60")
            println("BigFloat algebra branch ", branch, ": polynomial=", polynomial_error,
                    ", vector=", vector_error)
        end
        true
    end
end


# ============================================================================
# 6. MAIN EXECUTION
# ============================================================================

function main()
    println("\n" * "█" ^ 70)
    println("TWO-LEVEL REDUCED PAIRING TOY MODEL (Julia)")
    println("Source correction status: previously UNTESTED; this run checks Float64.")
    println("No molecular H2 or rigorous interval arithmetic claim.")
    println("█" ^ 70 * "\n")
    
    # Float64 test
    println("\n### FLOAT64 VERIFICATION ###\n")
    eps1, eps2, g = -1.0, -0.5, 0.3
    results_f64 = compare_bethe_vs_fci(eps1, eps2, g, verbose=true)
    @assert results_f64["passed"]
    @assert isapprox(results_f64["E_bethe"], -0.883095189484530; atol=1e-12)
    @assert isapprox(results_f64["E_pair"], -2.383095189484530; atol=1e-12)
    for parameters in ((-1.0, -0.5, 0.0), (-1.0, -1.0, 0.3), (0.2, 0.2, 0.0))
        @assert compare_bethe_vs_fci(parameters...; verbose=false)["passed"]
    end
    if "--bigfloat" in ARGS
        @assert test_arbitrary_precision()
    else
        println("BigFloat algebra check skipped; request with --bigfloat.")
    end
    println("PASS: requested toy-model checks completed; this is not a proof of molecular accuracy.")
end

# Run main if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
