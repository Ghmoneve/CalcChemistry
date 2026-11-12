"""Finite-dimensional theorem regression tests, not molecular benchmarks."""
import math
import unittest
import numpy as np


class CertificateTests(unittest.TestCase):
    def test_sharp_two_by_two_bounds_with_real_and_complex_coupling(self):
        from experiments.certificates import block_energy_bounds
        for a, c in ((-2.0, 3.0), (3.0, -2.0), (1.0, 1.0)):
            for coupling in (0.0, -0.75, 0.3 + 0.4j):
                with self.subTest(a=a, c=c, coupling=coupling):
                    H = np.array([[a, coupling], [np.conj(coupling), c]])
                    result = block_energy_bounds(H, 1)
                    reference = np.linalg.eigvalsh(H)[0]
                    self.assertAlmostEqual(result['lower'], reference, places=13)
                    self.assertEqual(result['upper'], a)
                    self.assertAlmostEqual(result['b'], abs(coupling), places=14)

    def test_two_by_two_extreme_couplings_do_not_square_first(self):
        from experiments.certificates import block_energy_bounds
        for coupling in (1e-200, 1e200):
            with self.subTest(coupling=coupling), np.errstate(all='raise'):
                result = block_energy_bounds([[0.0, coupling], [coupling, 0.0]], 1)
                self.assertTrue(math.isfinite(result['lower']))
                self.assertAlmostEqual(result['lower'] / coupling, -1.0, places=14)

    def test_seeded_hermitian_enclosures_through_dimension_ten(self):
        from experiments.certificates import block_energy_bounds
        rng = np.random.default_rng(20261112)
        for complex_entries in (False, True):
            for size in range(2, 11):
                for sample in range(4):
                    X = rng.normal(size=(size, size))
                    if complex_entries:
                        X = X + 1j * rng.normal(size=(size, size))
                    H = (X + X.conj().T) / 2
                    exact = np.linalg.eigvalsh(H)[0]
                    tolerance = 5e-12 * max(1.0, np.linalg.norm(H, 2))
                    for dimension in range(1, size):
                        with self.subTest(size=size, p=dimension, sample=sample,
                                          complex_entries=complex_entries):
                            result = block_energy_bounds(H, dimension)
                            self.assertLessEqual(result['lower'], exact + tolerance)
                            self.assertLessEqual(exact, result['upper'] + tolerance)

    def test_positive_complement_gap_controls_energy_error(self):
        from experiments.certificates import block_energy_bounds
        rng = np.random.default_rng(73)
        for size in range(2, 11):
            p = size // 2
            A = np.diag(np.linspace(-2.0, -1.0, p))
            C = np.diag(np.linspace(1.0, 2.0, size - p))
            B = 0.2 * (rng.normal(size=(p, size - p))
                       + 1j * rng.normal(size=(p, size - p)))
            H = np.block([[A, B], [B.conj().T, C]])
            result = block_energy_bounds(H, p)
            expected = np.linalg.norm(B, 2)**2 / 3.0
            with self.subTest(size=size):
                self.assertAlmostEqual(result['gap_error_bound'], expected, places=13)
                self.assertLessEqual(result['upper'] - result['exact'], expected + 5e-12)
                self.assertLessEqual(result['upper'] - result['lower'], expected + 5e-12)

    def test_omitted_low_energy_sector_is_not_certified(self):
        from experiments.certificates import block_energy_bounds
        # Active-state residual is zero, but an excluded state lies below it.
        H = np.diag([1.0, -2.0])
        result = block_energy_bounds(H, 1)
        self.assertAlmostEqual(result['lower'], -2.0)
        self.assertAlmostEqual(result['upper'], 1.0)
        self.assertIsNone(result['gap_error_bound'])
        self.assertAlmostEqual(result['exact'], -2.0)

    def test_invalid_hamiltonians_are_rejected(self):
        from experiments.certificates import block_energy_bounds
        invalid = ([1.0, 2.0], np.zeros((2, 3)), np.empty((0, 0)),
                   [[1.0, np.nan], [np.nan, 2.0]],
                   [[1.0, np.inf], [np.inf, 2.0]],
                   [[1.0, 1.0j], [1.0j, 2.0]],
                   [[1.0 + 1.0j, 0.0], [0.0, 2.0]])
        for matrix in invalid:
            with self.subTest(matrix=matrix), self.assertRaises(ValueError):
                block_energy_bounds(matrix, 1)

    def test_invalid_active_dimensions_are_rejected(self):
        from experiments.certificates import block_energy_bounds
        for dimension in (0, -1, 2, 1.5, True, None):
            with self.subTest(dimension=dimension), self.assertRaises(ValueError):
                block_energy_bounds(np.eye(2), dimension)


class BarrierIntervalTests(unittest.TestCase):
    def test_invalid_energy_intervals_are_rejected(self):
        import experiments.certificates as certificates
        self.assertTrue(callable(getattr(certificates, 'barrier_interval', None)),
                        'barrier_interval helper is missing')
        invalid = ((), (1.0,), (0.0, 1.0, 2.0), (2.0, 1.0),
                   (0.0, np.nan), (-np.inf, 1.0), (0.0, 1.0j),
                   [[0.0, 1.0]])
        for interval in invalid:
            with self.subTest(interval=interval):
                for reactant, transition in ((interval, (0.0, 1.0)),
                                             ((0.0, 1.0), interval)):
                    with self.assertRaises(ValueError):
                        certificates.barrier_interval(reactant, transition)

    def test_barrier_contains_all_endpoint_differences(self):
        import experiments.certificates as certificates
        self.assertTrue(callable(getattr(certificates, 'barrier_interval', None)),
                        'barrier_interval helper is missing')
        reactant, transition = (-3.0, -2.0), (-2.5, 0.0)
        interval = certificates.barrier_interval(reactant, transition)
        self.assertEqual(interval, (-0.5, 3.0))
        for initial in reactant:
            for final in transition:
                self.assertLessEqual(interval[0], final - initial)
                self.assertLessEqual(final - initial, interval[1])
        self.assertEqual(certificates.barrier_interval((1.0, 1.0), (3.0, 3.0)),
                         (2.0, 2.0))


class RateIntervalTests(unittest.TestCase):
    def test_invalid_rate_inputs_are_rejected(self):
        import experiments.certificates as certificates
        self.assertTrue(callable(getattr(certificates, 'rate_interval', None)),
                        'rate_interval helper is missing')
        invalid = [((2.0, 1.0), 2.0, 1.0), ((0.0, np.nan), 2.0, 1.0),
                   ((0.0, 1.0j), 2.0, 1.0), ((0.0,), 2.0, 1.0)]
        invalid += [((0.0, 1.0), A, 1.0) for A in (-1.0, np.inf, np.nan, 1.0j)]
        invalid += [((0.0, 1.0), 2.0, theta)
                    for theta in (0.0, -1.0, np.inf, np.nan, 1.0j)]
        for barrier, prefactor, thermal_energy in invalid:
            with self.subTest(parameters=(barrier, prefactor, thermal_energy)):
                with self.assertRaises(ValueError):
                    certificates.rate_interval(barrier, prefactor, thermal_energy)

    def test_rate_interval_reverses_barrier_endpoints(self):
        import experiments.certificates as certificates
        self.assertTrue(callable(getattr(certificates, 'rate_interval', None)),
                        'rate_interval helper is missing')
        interval = certificates.rate_interval((-0.5, 3.0), 7.0, 2.0)
        np.testing.assert_allclose(interval, (7.0 * math.exp(-1.5), 7.0 * math.exp(0.25)),
                                   rtol=2e-15, atol=0.0)
        for barrier in np.linspace(-0.5, 3.0, 51):
            rate = 7.0 * math.exp(-barrier / 2.0)
            self.assertLessEqual(interval[0] - 1e-14, rate)
            self.assertLessEqual(rate, interval[1] + 1e-14)
        self.assertEqual(certificates.rate_interval((1.0, 2.0), 0.0, 1.0), (0.0, 0.0))
        self.assertEqual(certificates.rate_interval((0.0, 0.0), 3.0, 1.0), (3.0, 3.0))

    def test_unrepresentable_positive_rates_are_rejected(self):
        import experiments.certificates as certificates
        self.assertTrue(callable(getattr(certificates, 'rate_interval', None)),
                        'rate_interval helper is missing')
        for barrier in ((-1000.0, -999.0), (999.0, 1000.0)):
            with self.subTest(barrier=barrier), self.assertRaises(ValueError):
                certificates.rate_interval(barrier, 1.0, 1.0)


class ExperimentReportTests(unittest.TestCase):
    def get_report(self):
        import experiments.certificates as certificates
        self.assertTrue(callable(getattr(certificates, 'run_experiments', None)),
                        'run_experiments helper is missing')
        return certificates.run_experiments(seed=20261112, samples_per_dimension=4)

    def test_seeded_report_records_every_split_and_enclosure(self):
        report = self.get_report()
        second = self.get_report()
        self.assertEqual(report['random_blocks'], second['random_blocks'])
        records = report['random_blocks']['cases']
        self.assertEqual(len(records), 2 * 4 * sum(range(1, 10)))
        self.assertEqual(report['metadata']['seed'], 20261112)
        self.assertFalse(report['metadata']['rigorous_interval_arithmetic'])
        self.assertTrue(report['passed'])
        for row in records:
            self.assertLessEqual(row['lower'] - row['exact'], row['tolerance'])
            self.assertLessEqual(row['exact'] - row['upper'], row['tolerance'])
        maximum = max(max(0.0, row['lower'] - row['exact'],
                          row['exact'] - row['upper']) for row in records)
        self.assertEqual(report['random_blocks']['max_enclosure_violation'], maximum)

    def test_fixed_ctmc_exponentials_obey_total_variation_bound(self):
        from scipy.linalg import expm
        report = self.get_report()['ctmc']
        Q, Q_tilde = np.asarray(report['Q']), np.asarray(report['Q_tilde'])
        for generator in (Q, Q_tilde):
            np.testing.assert_allclose(generator.sum(axis=1), 0.0, atol=1e-14)
            off_diagonal = generator.copy()
            np.fill_diagonal(off_diagonal, 0.0)
            self.assertGreaterEqual(off_diagonal.min(), 0.0)
        eta = 0.5 * np.max(np.sum(np.abs(Q - Q_tilde), axis=1))
        self.assertAlmostEqual(report['eta'], eta, places=14)
        self.assertEqual(len(report['cases']), 20)
        self.assertTrue(any(row['tv'] > 1e-5 for row in report['cases']))
        for row in report['cases']:
            initial, time = np.asarray(row['initial_distribution']), row['time']
            p, q = initial @ expm(time * Q), initial @ expm(time * Q_tilde)
            self.assertAlmostEqual(row['tv'], 0.5 * np.abs(p - q).sum(), places=14)
            self.assertLessEqual(row['tv'], time * eta + 5e-13)
            self.assertAlmostEqual(p.sum(), 1.0, places=12)
            self.assertAlmostEqual(q.sum(), 1.0, places=12)
            self.assertGreaterEqual(min(p.min(), q.min()), -5e-14)

    def test_synthetic_kinetics_report_has_explicit_units(self):
        report = self.get_report()['kinetics']
        self.assertTrue(report['synthetic_only'])
        self.assertEqual(report['energy_unit'], 'kJ/mol')
        self.assertEqual(report['rate_unit'], 's^-1')
        self.assertEqual(report['thermal_energy_unit'], report['energy_unit'])
        lower, upper = report['rate_interval']
        for value in report['sample_rates']:
            self.assertLessEqual(lower - 1e-8, value)
            self.assertLessEqual(value, upper + 1e-8)


if __name__ == '__main__':
    unittest.main()
