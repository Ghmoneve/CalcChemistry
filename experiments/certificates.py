"""Finite-dimensional energy enclosures; floating-point validation only.

The block routine uses dense exact-sector diagonalization as a benchmark. It is
not a scalable chemistry algorithm and does not certify floating-point rounding.
"""
import math
from typing import Any
import numpy as np


def block_energy_bounds(hamiltonian, active_dimension):
    """Benchmark a Hermitian block enclosure using dense eigensolvers.

    With a=min eig(A), c=min eig(C), b=||B||_2, the mathematical lower
    endpoint is (a+c-sqrt((a-c)^2+4b^2))/2 and the upper endpoint is a.
    The optional b^2/(c-a) error bound requires c>a. `exact` is a float64
    eigvalsh reference, not an exact arithmetic or interval certificate.
    Hermiticity is checked to absolute tolerance 1e-12.
    """
    H = np.asarray(hamiltonian, dtype=complex)
    if H.ndim != 2 or H.shape[0] != H.shape[1]:
        raise ValueError('Hamiltonian must be square')
    if not np.isfinite(H).all() or not np.allclose(H, H.conj().T, atol=1e-12, rtol=0):
        raise ValueError('Hamiltonian must be finite and Hermitian')
    p = active_dimension
    if isinstance(p, bool) or not isinstance(p, int) or not 0 < p < len(H):
        raise ValueError('Both active and complement spaces must be nonempty')
    A, B, C = H[:p, :p], H[:p, p:], H[p:, p:]
    a = float(np.linalg.eigvalsh(A)[0])
    c = float(np.linalg.eigvalsh(C)[0])
    b = float(np.linalg.norm(B, 2))
    half_gap = abs(c / 2 - a / 2)
    radius = math.hypot(half_gap, b)
    correction = 0.0 if b == 0 else b * (b / radius) / (1 + half_gap / radius)
    lower = min(a, c) - correction
    return {'lower': lower, 'upper': a, 'a': a, 'c': c, 'b': b,
            'gap_error_bound': b*(b/(c-a)) if c > a else None,
            'exact': float(np.linalg.eigvalsh(H)[0])}


def _finite_interval(interval):
    """Validate a closed, real, finite interval without silently dropping Im."""
    values = np.asarray(interval)
    if values.shape != (2,) or np.iscomplexobj(values):
        raise ValueError('An interval must contain exactly two real endpoints')
    try:
        lower, upper = map(float, values)
    except (TypeError, ValueError) as error:
        raise ValueError('Interval endpoints must be real numbers') from error
    if not math.isfinite(lower) or not math.isfinite(upper) or lower > upper:
        raise ValueError('Interval endpoints must be finite and ordered')
    return lower, upper


def barrier_interval(reactant_interval, transition_interval):
    """Bound E_transition - E_reactant in the input energy units.

    Inputs must already enclose their respective energies. Negative barriers
    are retained; these helpers do not validate a transition-state model.
    Arithmetic is ordinary float64, without directed rounding.
    """
    r_lower, r_upper = _finite_interval(reactant_interval)
    t_lower, t_upper = _finite_interval(transition_interval)
    return _finite_interval((t_lower - r_upper, t_upper - r_lower))


def rate_interval(barrier, prefactor, thermal_energy):
    """Map [L,U] to [A exp(-U/theta), A exp(-L/theta)].

    theta is k_B*T (per particle) or R*T (per mole), in the SAME units as
    the barrier. A>=0 is fixed, in inverse time units; theta>0 is fixed.
    This conditional Arrhenius model includes no prefactor/temperature
    uncertainty or dynamical validity claim. Nonfinite or underflowed
    positive rates raise ValueError instead of returning false endpoints.
    """
    lower, upper = _finite_interval(barrier)
    scalars = []
    for value in (prefactor, thermal_energy):
        if np.ndim(value) != 0 or np.iscomplexobj(value):
            raise ValueError('Prefactor and thermal energy must be real scalars')
        try:
            value = float(value)
        except (TypeError, ValueError) as error:
            raise ValueError('Prefactor and thermal energy must be real scalars') from error
        if not math.isfinite(value):
            raise ValueError('Prefactor and thermal energy must be finite')
        scalars.append(value)
    prefactor, thermal_energy = scalars
    if prefactor < 0 or thermal_energy <= 0:
        raise ValueError('Require prefactor >= 0 and thermal energy > 0')
    if prefactor == 0:
        return 0.0, 0.0
    try:
        rates = tuple(prefactor if endpoint == 0 else
                      math.exp(math.log(prefactor) - endpoint / thermal_energy)
                      for endpoint in (upper, lower))
    except OverflowError as error:
        raise ValueError('Rate is outside the supported float64 range') from error
    if not all(math.isfinite(rate) and rate > 0 for rate in rates):
        raise ValueError('Rate is outside the supported float64 range')
    return rates


def run_experiments(seed=20261112, samples_per_dimension=4) -> dict[str, Any]:
    """Execute reproducible synthetic checks; return JSON-compatible evidence."""
    import platform
    from datetime import datetime, timezone
    import scipy
    from scipy.linalg import expm

    rng = np.random.default_rng(seed)
    rows = []
    for complex_entries in (False, True):
        for size in range(2, 11):
            for sample in range(samples_per_dimension):
                X = rng.normal(size=(size, size))
                if complex_entries:
                    X = X + 1j * rng.normal(size=(size, size))
                H = (X + X.conj().T) / 2
                tolerance = 5e-12 * max(1.0, float(np.linalg.norm(H, 2)))
                for p in range(1, size):
                    result = block_energy_bounds(H, p)
                    violation = max(0.0, result['lower'] - result['exact'],
                                    result['exact'] - result['upper'])
                    gap = result['gap_error_bound']
                    gap_violation = (max(0.0, result['upper'] - result['exact'] - gap)
                                     if gap is not None else None)
                    rows.append(dict(result, dimension=size, active_dimension=p,
                                     complex_entries=complex_entries, sample=sample,
                                     tolerance=tolerance, enclosure_violation=violation,
                                     gap_violation=gap_violation))

    sharp_cases = []
    for a, c in ((-2.0, 3.0), (3.0, -2.0), (1.0, 1.0)):
        for coupling in (0.0, -0.75, 0.3 + 0.4j):
            H = [[a, coupling], [np.conj(coupling), c]]
            result = block_energy_bounds(H, 1)
            sharp_cases.append(dict(result, coupling_real=float(np.real(coupling)),
                                    coupling_imag=float(np.imag(coupling)),
                                    sharpness_error=abs(result['lower'] - result['exact'])))
    counterexample = dict(block_energy_bounds(np.diag([1.0, -2.0]), 1),
                          active_state_residual=0.0, matrix_diagonal=[1.0, -2.0])

    reactant, transition = (-100.1, -99.9), (-80.2, -79.8)
    barrier = barrier_interval(reactant, transition)
    # Synthetic molar energy example: theta=RT with R in kJ/(mol K).
    temperature, gas_constant, prefactor = 300.0, 0.00831446261815324, 1e6
    thermal_energy = gas_constant * temperature
    rates = rate_interval(barrier, prefactor, thermal_energy)
    sample_barriers = np.linspace(*barrier, 21).tolist()
    sample_rates = [prefactor * math.exp(-value / thermal_energy)
                    for value in sample_barriers]
    rate_violation = max(max(0.0, rates[0] - rate, rate - rates[1]) for rate in sample_rates)
    rate_tolerance = 5e-12 * max(1.0, rates[1])
    kinetics = dict(synthetic_only=True, energy_unit='kJ/mol', rate_unit='s^-1',
                    thermal_energy_unit='kJ/mol', temperature_unit='K',
                    reactant_interval=list(reactant), transition_interval=list(transition),
                    barrier_interval=list(barrier), rate_interval=list(rates),
                    temperature=temperature, gas_constant=gas_constant,
                    gas_constant_unit='kJ/(mol K)', thermal_energy=thermal_energy,
                    prefactor=prefactor, sample_barriers=sample_barriers,
                    sample_rates=sample_rates, max_rate_enclosure_violation=rate_violation,
                    tolerance=rate_tolerance,
                    assumption='Fixed A and T; k=A exp(-barrier/(RT)); synthetic only')

    # Row generators, with rates in s^-1. Neither is inferred from molecular data.
    Q = np.array([[-1.2, 1.0, 0.2], [0.4, -1.0, 0.6], [0.3, 0.2, -0.5]])
    Q_tilde = np.array([[-1.25, 1.05, 0.2], [0.38, -1.0, 0.62], [0.31, 0.2, -0.51]])
    eta = float(0.5 * np.max(np.sum(np.abs(Q - Q_tilde), axis=1)))
    ctmc_rows = []
    initials = [*np.eye(3), np.array([0.2, 0.3, 0.5])]
    for time in (0.0, 0.01, 0.1, 1.0, 5.0):
        P, P_tilde = expm(time * Q), expm(time * Q_tilde)
        for initial in initials:
            p, q = initial @ P, initial @ P_tilde
            tv = float(0.5 * np.sum(np.abs(p - q)))
            ctmc_rows.append(dict(time=time, initial_distribution=initial.tolist(),
                                  distribution=p.tolist(), perturbed_distribution=q.tolist(),
                                  tv=tv, bound=time * eta, violation=max(0.0, tv - time * eta),
                                  mass_error=float(max(abs(p.sum()-1), abs(q.sum()-1))),
                                  minimum_probability=float(min(p.min(), q.min()))))
    ctmc = dict(synthetic_only=True, convention='row distributions: p(t)=p(0) exp(t Q)',
                generator_unit='s^-1', time_unit='s', eta_unit='s^-1',
                Q=Q.tolist(), Q_tilde=Q_tilde.tolist(), eta=eta,
                eta_definition='0.5 * max_i sum_j abs(Q_ij - Q_tilde_ij)',
                tolerance=5e-13, cases=ctmc_rows,
                max_bound_violation=max(row['violation'] for row in ctmc_rows))
    passed = (all(row['enclosure_violation'] <= row['tolerance']
                  and (row['gap_violation'] is None or row['gap_violation'] <= row['tolerance'])
                  for row in rows)
              and all(row['sharpness_error'] <= 5e-12 for row in sharp_cases)
              and counterexample['lower'] == -2.0 and counterexample['upper'] == 1.0
              and counterexample['gap_error_bound'] is None
              and rate_violation <= rate_tolerance
              and all(row['violation'] <= 5e-13 and row['mass_error'] <= 5e-13
                      and row['minimum_probability'] >= -5e-13 for row in ctmc_rows))
    return dict(
        metadata=dict(schema_version=1, generated_at_utc=datetime.now(timezone.utc).isoformat(),
                      python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                      seed=seed, bit_generator=type(rng.bit_generator).__name__,
                      dimensions=list(range(2, 11)), samples_per_dimension=samples_per_dimension,
                      entry_distribution='independent standard normal real/imaginary; H=(X+X*)/2',
                      block_energy_unit='arbitrary model energy unit',
                      tolerance_rule='5e-12 * max(1, spectral_norm(H))',
                      rigorous_interval_arithmetic=False,
                      reference='numpy.linalg.eigvalsh float64; no outward rounding',
                      scope='Synthetic finite matrices and kinetics; no molecular H2 data',
                      reproduce='/usr/bin/python3 -B -m experiments.certificates'),
        random_blocks=dict(case_count=len(rows),
                           max_enclosure_violation=max(row['enclosure_violation'] for row in rows),
                           max_gap_violation=max(row['gap_violation'] or 0.0 for row in rows),
                           cases=rows),
        sharp_two_by_two=dict(cases=sharp_cases,
                             max_sharpness_error=max(row['sharpness_error'] for row in sharp_cases)),
        omitted_sector_counterexample=counterexample, kinetics=kinetics, ctmc=ctmc,
        passed=bool(passed))


if __name__ == '__main__':
    import json
    from pathlib import Path

    with np.errstate(all='raise'):
        report = run_experiments()
    output = Path(__file__).with_name('certificate_results.json')
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps({'output': str(output), 'passed': report['passed'],
                      'random_block_cases': report['random_blocks']['case_count'],
                      'max_enclosure_violation': report['random_blocks']['max_enclosure_violation'],
                      'ctmc_cases': len(report['ctmc']['cases']),
                      'max_ctmc_bound_violation': report['ctmc']['max_bound_violation']}, indent=2))
    raise SystemExit(0 if report['passed'] else 1)
