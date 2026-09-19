"""Unit tests for the parser, exact null-space solver and balancer.

Runnable with either pytest or the standard library:
    python -m unittest tests.test_chemistry -v
"""

import unittest
from fractions import Fraction

from app.chemistry import (
    BalanceError,
    FormulaError,
    balance,
    conservation_matrix,
    null_space,
    parse_formula,
    verify,
    ordered_elements,
    parse_all,
)


class TestParseFormula(unittest.TestCase):
    def test_simple(self):
        counts, charge = parse_formula("H2O")
        self.assertEqual(counts, {"H": 2, "O": 1})
        self.assertEqual(charge, 0)

    def test_nested_parentheses(self):
        counts, _ = parse_formula("Ca3(PO4)2")
        self.assertEqual(counts, {"Ca": 3, "P": 2, "O": 8})

    def test_deeply_nested_brackets(self):
        counts, _ = parse_formula("K4[Fe(CN)6]")
        self.assertEqual(counts, {"K": 4, "Fe": 1, "C": 6, "N": 6})

    def test_mixed_bracket_types(self):
        counts, _ = parse_formula("{[(NH4)2]3}2")
        self.assertEqual(counts, {"N": 12, "H": 48})

    def test_hydrate_unicode_dot(self):
        counts, _ = parse_formula("CuSO4·5H2O")
        self.assertEqual(counts, {"Cu": 1, "S": 1, "O": 9, "H": 10})

    def test_hydrate_ascii_dot(self):
        counts, _ = parse_formula("CuSO4.5H2O")
        self.assertEqual(counts, {"Cu": 1, "S": 1, "O": 9, "H": 10})

    def test_hydrate_multiple_segments(self):
        counts, _ = parse_formula("Na2CO3·10H2O")
        self.assertEqual(counts, {"Na": 2, "C": 1, "O": 13, "H": 20})

    def test_charges(self):
        self.assertEqual(parse_formula("SO4^2-")[1], -2)
        self.assertEqual(parse_formula("Fe^3+")[1], 3)
        self.assertEqual(parse_formula("H^+")[1], 1)
        self.assertEqual(parse_formula("OH^-")[1], -1)
        self.assertEqual(parse_formula("Al^{3+}")[1], 3)
        self.assertEqual(parse_formula("MnO4^-")[0], {"Mn": 1, "O": 4})

    def test_whitespace_is_ignored(self):
        counts, charge = parse_formula(" CuSO4 · 5H2O ")
        self.assertEqual(counts, {"Cu": 1, "S": 1, "O": 9, "H": 10})

    def test_invalid_formulas(self):
        for bad in ["", "Xx2", "H2O)", "(H2O", "2", "H0", "SO4^2", "h2O", "CuSO4."]:
            with self.assertRaises(FormulaError, msg=bad):
                parse_formula(bad)


class TestNullSpace(unittest.TestCase):
    def test_one_dimensional(self):
        # 2 H2 + O2 -> 2 H2O  (columns: H2, O2, -H2O)
        matrix = [
            [Fraction(2), Fraction(0), Fraction(-2)],
            [Fraction(0), Fraction(2), Fraction(-1)],
        ]
        basis = null_space(matrix)
        self.assertEqual(len(basis), 1)
        # Free variable (H2O column) fixed to 1; scaling by 2 gives 2,1,2.
        self.assertEqual(basis[0], [Fraction(1), Fraction(1, 2), Fraction(1)])

    def test_two_dimensional(self):
        matrix = [[Fraction(1), Fraction(0), Fraction(-1), Fraction(-1)]]
        self.assertEqual(len(null_space(matrix)), 3)

    def test_trivial_only(self):
        matrix = [[Fraction(1), Fraction(0)], [Fraction(0), Fraction(1)]]
        self.assertEqual(null_space(matrix), [])


class TestBalance(unittest.TestCase):
    def test_water(self):
        result = balance(["H2", "O2"], ["H2O"])
        self.assertEqual(result.reactant_coefficients, [2, 1])
        self.assertEqual(result.product_coefficients, [2])

    def test_permanganate(self):
        result = balance(["KMnO4", "HCl"], ["KCl", "MnCl2", "H2O", "Cl2"])
        self.assertEqual(result.reactant_coefficients, [2, 16])
        self.assertEqual(result.product_coefficients, [2, 2, 8, 5])

    def test_combustion(self):
        result = balance(["C7H16", "O2"], ["CO2", "H2O"])
        self.assertEqual(result.reactant_coefficients, [1, 11])
        self.assertEqual(result.product_coefficients, [7, 8])

    def test_ionic_redox(self):
        result = balance(
            ["MnO4^-", "H^+", "Fe^2+"],
            ["Mn^2+", "Fe^3+", "H2O"],
        )
        self.assertEqual(result.reactant_coefficients, [1, 8, 5])
        self.assertEqual(result.product_coefficients, [1, 5, 4])

    def test_hydrate_decomposition(self):
        result = balance(["CuSO4·5H2O"], ["CuSO4", "H2O"])
        self.assertEqual(result.reactant_coefficients, [1])
        self.assertEqual(result.product_coefficients, [1, 5])

    def test_no_solution(self):
        with self.assertRaises(BalanceError) as ctx:
            balance(["H2", "O2"], ["CO2"])
        self.assertEqual(ctx.exception.code, "no_solution")
        self.assertEqual(ctx.exception.context["null_space_dimension"], 0)

    def test_underdetermined(self):
        with self.assertRaises(BalanceError) as ctx:
            balance(["C", "O2"], ["CO", "CO2"])
        self.assertEqual(ctx.exception.code, "underdetermined")
        self.assertEqual(ctx.exception.context["null_space_dimension"], 2)

    def test_not_all_positive_zero_coefficient(self):
        # N2 is a spectator: the unique exact solution sets its coefficient to 0.
        with self.assertRaises(BalanceError) as ctx:
            balance(["H2", "O2", "N2"], ["H2O"])
        self.assertEqual(ctx.exception.code, "not_all_positive")

    def test_not_all_positive_mixed_signs(self):
        # H2O2 -> H2 + O2 balances, but as posed a species must switch sides.
        with self.assertRaises(BalanceError) as ctx:
            balance(["H2", "H2O2"], ["O2"])
        self.assertEqual(ctx.exception.code, "not_all_positive")

    def test_invalid_formula_propagates(self):
        with self.assertRaises(FormulaError):
            balance(["H2", "Xx"], ["H2O"])


class TestVerify(unittest.TestCase):
    def test_balanced(self):
        reactants = parse_all(["H2", "O2"])
        products = parse_all(["H2O"])
        elements = ordered_elements(reactants + products)
        report = verify(reactants, products, [2, 1], [2], elements)
        self.assertTrue(report["balanced"])
        self.assertTrue(all(e["balanced"] for e in report["elements"]))
        self.assertTrue(report["charge"]["balanced"])

    def test_unbalanced(self):
        reactants = parse_all(["H2", "O2"])
        products = parse_all(["H2O"])
        elements = ordered_elements(reactants + products)
        report = verify(reactants, products, [1, 1], [1], elements)
        self.assertFalse(report["balanced"])

    def test_charge_imbalance_detected(self):
        reactants = parse_all(["Fe^2+"])
        products = parse_all(["Fe^3+"])
        elements = ordered_elements(reactants + products)
        report = verify(reactants, products, [1], [1], elements)
        self.assertFalse(report["balanced"])
        self.assertFalse(report["charge"]["balanced"])
        self.assertTrue(all(e["balanced"] for e in report["elements"]))


class TestConservationMatrix(unittest.TestCase):
    def test_charge_row_and_signs(self):
        reactants = parse_all(["MnO4^-", "H^+"])
        products = parse_all(["Mn^2+", "H2O"])
        elements = ordered_elements(reactants + products)
        matrix = conservation_matrix(reactants, products, elements)
        self.assertEqual(len(matrix), len(elements) + 1)  # charge row appended
        charge_row = matrix[-1]
        self.assertEqual(charge_row[0], Fraction(-1))  # MnO4^-
        self.assertEqual(charge_row[1], Fraction(1))   # H^+
        self.assertEqual(charge_row[2], Fraction(-2))  # Mn^2+ negated (product)
        self.assertEqual(charge_row[3], Fraction(0))   # H2O


if __name__ == "__main__":
    unittest.main()
