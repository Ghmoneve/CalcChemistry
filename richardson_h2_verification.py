"""Numerical checks for a two-level, one-pair reduced pairing toy model.

The historical ``h2``/``fci`` callable names are retained for compatibility.
No molecular integrals, nuclear energy, or ab initio H2 accuracy are implied.
All Hamiltonians and solver energies use H_z unless explicitly called E_pair:
    H_z = sum_p eps_p (n_pu + n_pd - 1) - g sum_pq P_p^dagger P_q
    H_N = H_z + (eps1 + eps2) I,       E_pair = E_z + eps1 + eps2.
See richardson_h2_proof.md for the derivation and boundary cases.
"""

import argparse
from itertools import combinations

import numpy as np


PARAMETER_CASES = (
    (-0.5, -0.3, 0.1, "weak coupling"),
    (-1.0, -0.5, 0.5, "moderate coupling"),
    (-1.0, -0.2, 1.0, "strong coupling"),
    (-1.0, -1.0, 0.3, "degenerate orbitals"),
    (-2.0, -0.1, 0.8, "large splitting"),
    (-1.0, -0.5, 0.0, "zero coupling"),
    (0.2, 0.2, 0.0, "fully degenerate"),
    (-0.5, -1.0, 0.3, "reversed orbital order"),
    (-1.0, -0.5, 1e-14, "near a Bethe pole"),
    (0.0, 1e-14, 1.0, "nearly degenerate orbitals"),
)


def _validate_parameters(eps1, eps2, g):
    """Restrict this benchmark to finite real parameters with g >= 0."""
    if not np.all(np.isfinite([eps1, eps2, g])) or g < 0:
        raise ValueError("eps1 and eps2 must be finite; g must be finite and >= 0")


def solve_richardson_equation_h2(eps1, eps2, g):
    """Return the lower and upper Richardson pair energies for g >= 0.

    H_z's paired block is [[eps1-eps2-g, -g], [-g, eps2-eps1-g]].
    Number-convention Richardson roots are E_pair = E_z + eps1 + eps2.
    The upper paired energy need not be the first excited energy of the
    complete six-state N=2 sector, which also has four unpaired zeros.
    """
    _validate_parameters(eps1, eps2, g)
    splitting = eps1 - eps2
    radius = np.hypot(splitting, g)
    ground = eps1 + eps2 - g - radius
    excited = eps1 + eps2 - g + radius
    return ground, excited


def verify_richardson_equation(E, eps1, eps2, g):
    """Return |1 - g sum_p 1/(2 eps_p - E)| for NUMBER-convention E.

    Return None at g=0 or an exact pole: the rational equation is undefined
    there. Such cases must instead be checked by a matrix eigenvector
    residual; None is not a claim that a candidate is a valid eigenvalue.
    """
    _validate_parameters(eps1, eps2, g)
    if not np.isfinite(E):
        raise ValueError("E_pair must be finite")
    denominators = (2.0 * eps1 - E, 2.0 * eps2 - E)
    if g == 0.0 or 0.0 in denominators:
        return None
    return abs(1.0 - g / denominators[0] - g / denominators[1])


def construct_bethe_wavefunction(E, eps1, eps2, g, *, branch=None):
    """Return (N=2 wavefunction, pair coefficients) for a NUMBER-convention energy.

    g is explicit. branch=0/1 selects the lower/upper paired root; when
    both roots coincide (g=0 and eps1=eps2), these select the two paired
    determinants. Without branch the lower state is the degenerate default.
    Stable null-vector ratios avoid division by g and Bethe-root poles.
    """
    roots = solve_richardson_equation_h2(eps1, eps2, g)
    if branch is None:
        branch = int(abs(E - roots[1]) < abs(E - roots[0]))
    if branch not in (0, 1):
        raise ValueError("branch must be 0 (lower) or 1 (upper)")
    scale = max(1.0, abs(eps1 - eps2), g, abs(E))
    if not np.isfinite(E) or abs(E - roots[branch]) > 64 * np.finfo(float).eps * scale:
        raise ValueError("E must be the selected number-convention pair energy")

    splitting = eps1 - eps2
    if g == 0.0:
        ground = np.array([1.0, 0.0]) if splitting <= 0.0 else np.array([0.0, 1.0])
    else:
        radius = np.hypot(splitting, g)
        # g/(radius + |splitting|), scaled to avoid overflow/cancellation.
        ratio = (g / radius) / (1.0 + abs(splitting) / radius)
        ground = np.array([ratio, 1.0]) if splitting >= 0.0 else np.array([1.0, ratio])
        ground /= np.linalg.norm(ground)
    coeffs = ground if branch == 0 else np.array([-ground[1], ground[0]])
    basis, lookup = build_fock_basis()
    wf = np.zeros(len(basis), dtype=complex)
    wf[lookup[(1, 1, 0, 0)]] = coeffs[0]
    wf[lookup[(0, 0, 1, 1)]] = coeffs[1]
    return wf, coeffs


def build_fock_basis(n_orb=2, n_elec=2):
    """Return canonical occupation tuples and their indices at fixed N.

    Modes are (1 up, 1 down, 2 up, 2 down, ...). For K=N=2 the paired
    determinants occupy indices 0 and 5, not 0 and 3.
    """
    n_spin_orb = 2 * n_orb
    basis = []
    for occupied in combinations(range(n_spin_orb), n_elec):
        state = [0] * n_spin_orb
        for i in occupied:
            state[i] = 1
        basis.append(tuple(state))
    return basis, {state: i for i, state in enumerate(basis)}


def _apply_fermion(state, mode, create):
    """Apply one CAR operator, with parity from earlier occupied modes."""
    if state[mode] == int(create):
        return None, 0
    sign = (-1) ** sum(state[:mode])
    result = list(state)
    result[mode] = int(create)
    return tuple(result), sign


def hamiltonian_pairing_fci(eps1, eps2, g, basis, basis_idx):
    """Build the two-orbital H_z with sequential CAR operators, including p=q.

    Matrix columns are input determinants. Creation is checked AFTER
    annihilation, so P_p^dagger P_p contributes -g on an occupied pair.
    This is exact diagonalization of a toy model, not molecular FCI.
    """
    _validate_parameters(eps1, eps2, g)
    H = np.zeros((len(basis), len(basis)))
    for column, state in enumerate(basis):
        H[column, column] = sum(
            energy * (state[2*p] + state[2*p + 1] - 1)
            for p, energy in enumerate((eps1, eps2))
        )
        for p in range(2):
            for q in range(2):
                current, sign = state, 1
                # Right-to-left action of a_pu^dag a_pd^dag a_qd a_qu.
                for mode, create in ((2*q, False), (2*q + 1, False),
                                     (2*p + 1, True), (2*p, True)):
                    current, factor = _apply_fermion(current, mode, create)
                    sign *= factor
                    if current is None:
                        break
                if current is not None:
                    H[basis_idx[current], column] -= g * sign
    return H


def compare_bethe_vs_fci(eps1, eps2, g, verbose=True, *, tolerance=1e-12):
    """Check the full spectrum, both paired vectors, and equation residuals.

    Energy/residual tolerances scale with max(1, ||H_z||_2). Normalization,
    orthogonality and subspace errors use the dimensionless tolerance.
    wf_error is distance to the numerical ground eigenspace, allowing
    arbitrary phases and degenerate bases. Raw rational residuals are
    diagnostic near poles; the scaled polynomial and vector residuals
    remain checked there. No raw rational residual is silently set to zero.
    """
    pair_roots = np.array(solve_richardson_equation_h2(eps1, eps2, g))
    hz_roots = pair_roots - (eps1 + eps2)
    basis, lookup = build_fock_basis()
    H = hamiltonian_pairing_fci(eps1, eps2, g, basis, lookup)
    eigenvalues, eigenvectors = np.linalg.eigh(H)
    scale = max(1.0, float(np.linalg.norm(H, 2)))
    states, coefficients = zip(*(
        construct_bethe_wavefunction(energy, eps1, eps2, g, branch=branch)
        for branch, energy in enumerate(pair_roots)
    ))
    wavefunctions = np.column_stack(states)
    residuals = np.linalg.norm(H @ wavefunctions - wavefunctions * hz_roots, axis=0)
    normalization_error = float(np.max(np.abs(np.linalg.norm(wavefunctions, axis=0) - 1.0)))
    orthogonality_error = float(abs(np.vdot(states[0], states[1])))
    ground_space = eigenvectors[:, np.abs(eigenvalues - eigenvalues[0]) <= tolerance * scale]
    projection = ground_space @ (ground_space.conj().T @ states[0])
    wf_error = float(np.linalg.norm(states[0] - projection))
    spectrum = np.sort(np.concatenate((hz_roots, np.zeros(4))))
    energy_error = float(abs(hz_roots[0] - eigenvalues[0]))
    spectrum_error = float(np.max(np.abs(spectrum - eigenvalues)))
    hermiticity_error = float(np.linalg.norm(H - H.conj().T))

    # d_p = 2 eps_p - E_pair, computed in centered coordinates to avoid
    # subtracting large common orbital offsets. This is the determinant,
    # not the rational secular equation at its excluded singular points.
    splitting = eps1 - eps2
    d1, d2 = (splitting - hz_roots) / scale, (-splitting - hz_roots) / scale
    characteristic_error = float(np.max(np.abs(d1 * d2 - (g / scale) * (d1 + d2))))
    richardson_errors = [verify_richardson_equation(E, eps1, eps2, g) for E in pair_roots]
    secular_checked, secular_passed = [], True
    for E, error in zip(pair_roots, richardson_errors):
        denominators = np.array([2 * eps1 - E, 2 * eps2 - E])
        pole_scale = max(1.0, abs(2 * eps1), abs(2 * eps2), abs(E))
        well_conditioned = error is not None and min(abs(denominators)) > np.sqrt(np.finfo(float).eps) * pole_scale
        secular_checked.append(well_conditioned)
        if well_conditioned:
            relative_error = error / (1.0 + np.sum(np.abs(g / denominators)))
            secular_passed = secular_passed and np.isfinite(relative_error) and relative_error <= tolerance

    dimensionless_errors = (normalization_error, orthogonality_error, wf_error, characteristic_error)
    energy_errors = (energy_error, spectrum_error, float(max(residuals)), hermiticity_error)
    passed = bool(
        all(np.isfinite(error) and error <= tolerance for error in dimensionless_errors)
        and all(np.isfinite(error) and error <= tolerance * scale for error in energy_errors)
        and secular_passed
    )
    results = {
        "E_bethe": hz_roots[0], "E_bethe_excited": hz_roots[1],
        "E_pair": pair_roots[0], "pair_roots": pair_roots,
        "E_fci": eigenvalues[0], "fci_spectrum": eigenvalues,
        "energy_error": energy_error, "spectrum_error": spectrum_error,
        "wf_bethe": coefficients[0], "wf_bethe_full": states[0],
        "wf_fci": eigenvectors[:, 0], "wf_error": wf_error,
        "eigenvector_residual": float(max(residuals)), "eigenvector_residuals": residuals,
        "normalization_error": normalization_error, "orthogonality_error": orthogonality_error,
        "hermiticity_error": hermiticity_error, "characteristic_error": characteristic_error,
        "richardson_error": richardson_errors[0], "richardson_errors": richardson_errors,
        "secular_checked": secular_checked, "passed": passed,
    }
    if verbose:
        print("Toy reduced pairing benchmark (two orbitals, N=2)")
        print(f"Parameters: eps1={eps1:g}, eps2={eps2:g}, g={g:g}")
        print(f"H_z ground: analytic={hz_roots[0]:.15f}, diagonalized={eigenvalues[0]:.15f}")
        print(f"H_N ground / Richardson pair energy: {pair_roots[0]:.15f}")
        print(f"Full H_z spectrum: {eigenvalues}")
        for key in ("energy_error", "spectrum_error", "eigenvector_residual", "normalization_error",
                    "orthogonality_error", "wf_error", "characteristic_error"):
            print(f"  {key}: {results[key]:.3e}")
        print(f"Rational residuals: {richardson_errors}; well-conditioned checks: {secular_checked}")
        print("PASS" if passed else "FAIL")
    return results


def test_multiple_parameters():
    """Run the deterministic CLI cases and return whether every check passed."""
    passed_count = 0
    for eps1, eps2, g, description in PARAMETER_CASES:
        results = compare_bethe_vs_fci(eps1, eps2, g, verbose=False)
        passed_count += int(results["passed"])
        print(f"{'PASS' if results['passed'] else 'FAIL'} {description}: "
              f"spectrum={results['spectrum_error']:.2e}, "
              f"vector={results['eigenvector_residual']:.2e}, "
              f"norm={results['normalization_error']:.2e}")
    print(f"Parameter cases passed: {passed_count}/{len(PARAMETER_CASES)}")
    return passed_count == len(PARAMETER_CASES)


def visualize_wavefunction_composition(eps1, eps2, g_values, output_path=None):
    """Return a figure; save only to an explicitly supplied output_path."""
    import matplotlib.pyplot as plt

    g_values = list(g_values)
    energies, coefficients = [], []
    for g in g_values:
        pair_energy, _ = solve_richardson_equation_h2(eps1, eps2, g)
        _, coeffs = construct_bethe_wavefunction(pair_energy, eps1, eps2, g)
        energies.append(pair_energy - (eps1 + eps2))
        coefficients.append(coeffs)
    probabilities = np.abs(np.asarray(coefficients).reshape((-1, 2)))**2
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    for index in range(2):
        ax1.plot(g_values, probabilities[:, index], "o-", label=f"orbital {index + 1}")
    ax1.set(xlabel="Pairing strength g", ylabel="Pair occupation probability",
            title="Toy reduced pairing wavefunction")
    ax1.legend()
    ax2.plot(g_values, energies, "o-")
    ax2.set(xlabel="Pairing strength g", ylabel="H_z ground energy (model units)",
            title="Toy reduced pairing energy")
    fig.tight_layout()
    if output_path is not None:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
    return fig


def main(argv=None):
    """Return a nonzero exit code if any benchmark check fails."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot", metavar="PATH", help="optionally save a wavefunction plot")
    args = parser.parse_args(argv)
    print("Scope: toy reduced pairing; no ab initio H2 or molecular-accuracy claim.")
    reference = compare_bethe_vs_fci(-1.0, -0.5, 0.3)
    cases_passed = test_multiple_parameters()
    passed = reference["passed"] and cases_passed
    if args.plot:
        import matplotlib.pyplot as plt
        figure = visualize_wavefunction_composition(-1.0, -0.5, np.linspace(0.0, 1.5, 31), args.plot)
        plt.close(figure)
    print("VERIFICATION PASS" if passed else "VERIFICATION FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
