"""Independent CAR/Jordan-Wigner tests for the two-level pairing toy model.

Run with /usr/bin/python3 -B -m unittest discover -s tests -p test_pairing.py -v.
The oracle uses tensor-product operators on all 16 Fock states, independently
of the production code's sequential occupation/sign updates.
"""

from itertools import product
import unittest

import numpy as np

import richardson_h2_verification as pairing


N2_BASIS = [
    (1, 1, 0, 0),
    (1, 0, 1, 0),
    (1, 0, 0, 1),
    (0, 1, 1, 0),
    (0, 1, 0, 1),
    (0, 0, 1, 1),
]


def jordan_wigner_annihilators():
    """Return a_i = Z^(tensor i) tensor |0><1| tensor I for four modes."""
    identity = np.eye(2)
    parity = np.diag([1.0, -1.0])
    lowering = np.array([[0.0, 1.0], [0.0, 0.0]])
    operators = []
    for mode in range(4):
        operator = np.ones((1, 1))
        for site in range(4):
            factor = parity if site < mode else lowering if site == mode else identity
            operator = np.kron(operator, factor)
        operators.append(operator)
    return operators


def independent_hamiltonian(eps1, eps2, g):
    """Build H_z on all Fock states, then project to the explicit N=2 basis."""
    annihilators = jordan_wigner_annihilators()
    numbers = [a.T @ a for a in annihilators]
    pair_annihilators = [annihilators[2*p + 1] @ annihilators[2*p] for p in range(2)]
    full = np.zeros((16, 16))
    for p, energy in enumerate((eps1, eps2)):
        full += energy * (numbers[2*p] + numbers[2*p + 1] - np.eye(16))
    for left in pair_annihilators:
        for right in pair_annihilators:
            full -= g * left.T @ right
    fock_states = list(product((0, 1), repeat=4))
    selected = [fock_states.index(state) for state in N2_BASIS]
    return full[np.ix_(selected, selected)], full, sum(numbers)


class FermionicHamiltonianTests(unittest.TestCase):
    def test_oracle_satisfies_canonical_anticommutation(self):
        operators = jordan_wigner_annihilators()
        for i, left in enumerate(operators):
            for j, right in enumerate(operators):
                with self.subTest(i=i, j=j):
                    np.testing.assert_array_equal(left @ right + right @ left, np.zeros((16, 16)))
                    np.testing.assert_array_equal(
                        left @ right.T + right.T @ left,
                        np.eye(16) if i == j else np.zeros((16, 16)),
                    )

    def test_full_n2_matrix_includes_diagonal_self_pairing(self):
        expected, full, number = independent_hamiltonian(-1.0, -0.5, 0.3)
        basis, lookup = pairing.build_fock_basis()
        self.assertEqual(basis, N2_BASIS)
        actual = pairing.hamiltonian_pairing_fci(-1.0, -0.5, 0.3, basis, lookup)
        np.testing.assert_allclose(full, full.T, atol=0, rtol=0)
        np.testing.assert_allclose(full @ number, number @ full, atol=0, rtol=0)
        np.testing.assert_allclose(actual, expected, atol=1e-15, rtol=0)


class AnalyticSpectrumTests(unittest.TestCase):
    def test_paired_roots_and_four_unpaired_states_match_car_spectrum(self):
        cases = [(-1.0, -0.5, 0.3), (0.2, -0.7, 0.0),
                 (-1.0, -1.0, 0.3), (0.4, 0.4, 0.0), (1.0, 2.0, 3.0)]
        for eps1, eps2, g in cases:
            with self.subTest(eps1=eps1, eps2=eps2, g=g):
                expected, _, _ = independent_hamiltonian(eps1, eps2, g)
                roots = pairing.solve_richardson_equation_h2(eps1, eps2, g)
                np.testing.assert_allclose(
                    np.sort([roots[0] - eps1 - eps2, roots[1] - eps1 - eps2,
                             0.0, 0.0, 0.0, 0.0]),
                    np.linalg.eigvalsh(expected), atol=2e-14, rtol=0,
                )
                shift = eps1 + eps2
                np.testing.assert_allclose(
                    np.array(roots),
                    np.linalg.eigvalsh(expected[np.ix_([0, 5], [0, 5])] + shift * np.eye(2)),
                    atol=2e-14, rtol=0,
                )


class RichardsonStateTests(unittest.TestCase):
    def test_positive_secular_sign_and_number_energy_convention(self):
        for energy in (-1.5 - 0.3 - np.hypot(0.5, 0.3),
                       -1.5 - 0.3 + np.hypot(0.5, 0.3)):
            with self.subTest(energy=energy):
                expected = abs(1.0 - 0.3 / (-2.0 - energy) - 0.3 / (-1.0 - energy))
                self.assertLess(expected, 1e-13)
                self.assertAlmostEqual(
                    pairing.verify_richardson_equation(energy, -1.0, -0.5, 0.3),
                    expected, delta=1e-15,
                )
        self.assertGreater(pairing.verify_richardson_equation(-3.0, -1.0, -0.5, 0.3), 0.1)

    def test_secular_equation_marks_zero_coupling_and_poles_not_applicable(self):
        self.assertIsNone(pairing.verify_richardson_equation(-2.0, -1.0, -0.5, 0.0))
        self.assertIsNone(pairing.verify_richardson_equation(-2.0, -1.0, -1.0, 0.3))

    def test_both_wavefunctions_are_normalized_eigenvectors_without_global_g(self):
        import inspect
        self.assertIn("g", inspect.signature(pairing.construct_bethe_wavefunction).parameters)
        cases = [(-1.0, -0.5, 0.3), (-0.5, -1.0, 0.3),
                 (-1.0, -0.5, 0.0), (-0.5, -1.0, 0.0),
                 (-1.0, -1.0, 0.3), (0.2, 0.2, 0.0),
                 (-1.0, -0.5, 1e-14), (0.0, 1e-14, 1.0)]
        for eps1, eps2, g in cases:
            with self.subTest(eps1=eps1, eps2=eps2, g=g):
                expected, _, _ = independent_hamiltonian(eps1, eps2, g)
                states = []
                for branch, energy in enumerate(pairing.solve_richardson_equation_h2(eps1, eps2, g)):
                    wf, coeffs = pairing.construct_bethe_wavefunction(
                        energy, eps1, eps2, g, branch=branch,
                    )
                    self.assertTrue(np.all(np.isfinite(wf)))
                    self.assertAlmostEqual(np.linalg.norm(wf), 1.0, delta=1e-14)
                    np.testing.assert_allclose(wf[[0, 5]], coeffs, atol=0, rtol=0)
                    np.testing.assert_array_equal(wf[1:5], np.zeros(4))
                    energy_z = energy - eps1 - eps2
                    np.testing.assert_allclose(expected @ wf, energy_z * wf, atol=2e-14, rtol=0)
                    states.append(wf)
                self.assertAlmostEqual(abs(np.vdot(*states)), 0.0, delta=1e-14)


class VerificationWorkflowTests(unittest.TestCase):
    def test_comparison_checks_full_spectrum_and_both_eigenvectors(self):
        for params in [(-1.0, -0.5, 0.3), (0.0, 0.0, 0.0),
                       (-1.0, -1.0, 0.3), (1.0, 2.0, 0.0)]:
            with self.subTest(params=params):
                results = pairing.compare_bethe_vs_fci(*params, verbose=False)
                self.assertTrue(results["passed"])
                for name in ("energy_error", "spectrum_error", "eigenvector_residual",
                             "normalization_error", "orthogonality_error", "wf_error"):
                    self.assertLess(results[name], 1e-12, name)
                oracle, _, _ = independent_hamiltonian(*params)
                np.testing.assert_allclose(results["fci_spectrum"], np.linalg.eigvalsh(oracle), atol=1e-14)
                self.assertAlmostEqual(results["E_pair"], results["E_bethe"] + sum(params[:2]))

    def test_invalid_model_parameters_are_rejected(self):
        basis, lookup = pairing.build_fock_basis()
        for params in [(0.0, 1.0, -0.1), (np.nan, 1.0, 0.1), (0.0, np.inf, 0.1)]:
            with self.subTest(params=params):
                with self.assertRaises(ValueError):
                    pairing.solve_richardson_equation_h2(*params)
                with self.assertRaises(ValueError):
                    pairing.hamiltonian_pairing_fci(*params, basis, lookup)

    def test_seeded_parameters_match_independent_matrix(self):
        rng = np.random.default_rng(73129)
        basis, lookup = pairing.build_fock_basis()
        for index in range(64):
            eps1, eps2 = rng.uniform(-3.0, 3.0, 2)
            g = rng.uniform(0.01, 2.0)
            with self.subTest(index=index):
                oracle, _, _ = independent_hamiltonian(eps1, eps2, g)
                actual = pairing.hamiltonian_pairing_fci(eps1, eps2, g, basis, lookup)
                np.testing.assert_allclose(actual, oracle, atol=2e-14, rtol=0)
                self.assertTrue(pairing.compare_bethe_vs_fci(eps1, eps2, g, verbose=False)["passed"])

    def test_cli_reports_success_and_toy_scope(self):
        import pathlib
        import subprocess
        import sys
        root = pathlib.Path(__file__).resolve().parents[1]
        completed = subprocess.run(
            [sys.executable, "-B", str(root / "richardson_h2_verification.py")],
            cwd=root, capture_output=True, text=True, timeout=30,
        )
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
        self.assertIn("toy reduced pairing", completed.stdout.lower())
        self.assertIn("PASS", completed.stdout)
        self.assertNotIn("PROOF COMPLETE", completed.stdout)

    def test_cli_fails_on_wrong_eigenvector_or_normalization(self):
        import pathlib
        import subprocess
        import sys
        self.assertTrue(callable(getattr(pairing, "main", None)), "a testable CLI entry point is required")
        root = pathlib.Path(__file__).resolve().parents[1]
        for corruption in ("np.roll(wf, 1)", "2.0 * wf"):
            with self.subTest(corruption=corruption):
                code = (
                    "import numpy as np\n"
                    "import richardson_h2_verification as pairing\n"
                    "original = pairing.construct_bethe_wavefunction\n"
                    "def corrupt(*args, **kwargs):\n"
                    "    wf, coeffs = original(*args, **kwargs)\n"
                    f"    return {corruption}, coeffs\n"
                    "pairing.construct_bethe_wavefunction = corrupt\n"
                    "raise SystemExit(pairing.main([]))\n"
                )
                completed = subprocess.run([sys.executable, "-B", "-c", code], cwd=root,
                                           capture_output=True, text=True, timeout=30)
                self.assertEqual(completed.returncode, 1, completed.stdout + completed.stderr)
                self.assertIn("FAIL", completed.stdout)
                self.assertNotIn("Traceback", completed.stderr)

    def test_plot_uses_each_requested_coupling(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig = pairing.visualize_wavefunction_composition(-1.0, -0.5, [0.0, 0.3, 1.0])
        self.assertIsNotNone(fig, "return the plot without writing an unsolicited file")
        try:
            expected = []
            for g in (0.0, 0.3, 1.0):
                oracle, _, _ = independent_hamiltonian(-1.0, -0.5, g)
                _, vectors = np.linalg.eigh(oracle)
                expected.append(abs(vectors[0, 0])**2)
            np.testing.assert_allclose(fig.axes[0].lines[0].get_ydata(), expected, atol=1e-14)
        finally:
            plt.close(fig)


if __name__ == "__main__":
    unittest.main()
