"""Exact counterexamples to historical general claims, using rational arithmetic.

No molecular energies or empirical DFT labels are used in these checks.
"""
from fractions import Fraction
from itertools import combinations, product
import unittest


def apply_ops(occupied, operators):
    """CAR operators listed in application order (rightmost acts first)."""
    current = set(occupied)
    sign = 1
    for mode, create in operators:
        if (mode in current) == create:
            return None, 0
        sign *= (-1) ** sum(i < mode for i in current)
        if create:
            current.add(mode)
        else:
            current.remove(mode)
    return tuple(sorted(current)), sign


def rdm_element(state_signs, creation, annihilation):
    """Equal-magnitude normalized real superposition; signs are ±1.

For D[I,J]=<a†_i1 ... a†_ip a_jp ... a_j1>, annihilation
acts in J order and creation in reversed I order.
"""
    total = Fraction(0)
    for ket, amplitude in state_signs.items():
        target, phase = apply_ops(ket, [(i, False) for i in annihilation]
                                 + [(i, True) for i in reversed(creation)])
        total += Fraction(amplitude * phase * state_signs.get(target, 0), len(state_signs))
    return total


class HistoricalCounterexamples(unittest.TestCase):
    def test_connected_component_count_is_not_independent_set_capacity(self):
        # Synthetic conflict graph on positions 0,3,6 with cutoff 3.5.
        positions = [0, 3, 6]
        edges = [(i, j) for i, j in combinations(range(3), 2)
                 if 2 * abs(positions[i] - positions[j]) <= 7]
        self.assertEqual(edges, [(0, 1), (1, 2)])
        reachable = {0}
        for _ in positions:
            for i, j in edges:
                if i in reachable or j in reachable:
                    reachable.update((i, j))
        self.assertEqual(len(reachable), 3)  # beta_0 = 1
        independent = [s for n in range(4) for s in combinations(range(3), n)
                       if all(not (i in s and j in s) for i, j in edges)]
        self.assertEqual(max(map(len, independent)), 2)

    def test_two_rdm_does_not_determine_three_rdm(self):
        plus = {(0, 1, 2): 1, (3, 4, 5): 1}
        minus = {(0, 1, 2): 1, (3, 4, 5): -1}
        # Includes repeated indices and all orders: full tensor, not a subset.
        for i, j, k, l in product(range(6), repeat=4):
            self.assertEqual(rdm_element(plus, (i, j), (k, l)),
                             rdm_element(minus, (i, j), (k, l)))
        self.assertEqual(rdm_element(plus, (0, 1, 2), (3, 4, 5)), Fraction(1, 2))
        self.assertEqual(rdm_element(minus, (0, 1, 2), (3, 4, 5)), Fraction(-1, 2))

    def test_proposed_sos_constraint_excludes_a_physical_state(self):
        psi = {(0, 1): 1, (2, 3): 1}
        d = rdm_element(psi, (0, 1), (0, 1))
        n0 = rdm_element(psi, (0,), (0,))
        n1 = rdm_element(psi, (1,), (1,))
        self.assertEqual(d*d, Fraction(1, 4))
        self.assertEqual(n0*n1*d, Fraction(1, 8))
        self.assertGreater(d*d, n0*n1*d)


if __name__ == '__main__':
    unittest.main()
