"""Chemical formula parsing and exact equation balancing.

Balancing method
----------------
Given ``m`` reactants and ``n`` products we build an exact rational matrix
with one row per element occurring in any species plus one extra row for the
total electric charge.  Reactant columns carry a positive sign, product
columns a negative sign, so every null-space vector of the matrix is a set of
coefficients that conserves each element and the charge.

The null space is computed with Gaussian elimination over the rationals
(:class:`fractions.Fraction`).  A chemically meaningful, unique answer exists
only when the null space is one-dimensional and its basis vector can be
scaled to strictly positive integers; the coefficients are then reduced by
their greatest common divisor.

Formula syntax
--------------
* Elements: standard symbols, e.g. ``H``, ``O``, ``Fe``.
* Counts: trailing integers, e.g. ``H2O``.
* Nested groups: ``()``, ``[]`` and ``{}`` with multipliers, arbitrarily
  nested, e.g. ``Ca3(PO4)2``, ``K4[Fe(CN)6]``.
* Hydrates / dot notation: ``.``, ``·``, ``•`` or ``*`` separate segments and
  each segment may carry a leading multiplier, e.g. ``CuSO4·5H2O``.
* Ionic charge: trailing ``^`` followed by a magnitude and a sign, e.g.
  ``SO4^2-``, ``Fe^3+``, ``H^+``, ``OH^-`` (``^{2-}`` is also accepted).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
from math import gcd
from typing import Dict, List, Sequence, Tuple

# ---------------------------------------------------------------------------
# Formula parsing
# ---------------------------------------------------------------------------

ELEMENTS = frozenset(
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni "
    "Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I "
    "Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt "
    "Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr "
    "Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og".split()
)

_OPEN_TO_CLOSE = {"(": ")", "[": "]", "{": "}"}
_CLOSERS = frozenset(_OPEN_TO_CLOSE.values())
_HYDRATE_DOTS = frozenset({".", "·", "•", "*"})

_CHARGE_RE = re.compile(r"^(?:(\d+)([+-])|([+-])(\d+)?)$")
_LEADING_COEFF_RE = re.compile(r"^(\d+)(.*)$", re.S)


class FormulaError(ValueError):
    """Raised when a chemical formula cannot be parsed."""


class BalanceError(Exception):
    """Raised when an equation has no unique all-positive balancing."""

    def __init__(self, code: str, message: str, **context):
        super().__init__(message)
        self.code = code
        self.message = message
        self.context = context


def _parse_charge(text: str, source: str) -> int:
    t = text.strip()
    if t.startswith("{") and t.endswith("}"):
        t = t[1:-1]
    m = _CHARGE_RE.match(t)
    if not m:
        raise FormulaError(f"invalid charge '^{text}' in formula '{source}'")
    if m.group(2) is not None:  # e.g. "2+" / "2-"
        magnitude, sign = int(m.group(1)), m.group(2)
    else:  # e.g. "+" / "-" / "+2"
        sign = m.group(3)
        magnitude = int(m.group(4)) if m.group(4) else 1
    return magnitude if sign == "+" else -magnitude


def _read_number(s: str, i: int, source: str) -> Tuple[int, int]:
    j = i
    while j < len(s) and s[j].isdigit():
        j += 1
    if j == i:
        return 1, i
    value = int(s[i:j])
    if value == 0:
        raise FormulaError(f"multipliers must be positive in formula '{source}'")
    return value, j


def _parse_group(s: str, i: int, closer: str | None, source: str) -> Tuple[Dict[str, int], int]:
    counts: Dict[str, int] = {}
    while i < len(s):
        ch = s[i]
        if closer is not None and ch == closer:
            return counts, i + 1
        if ch in _CLOSERS:
            raise FormulaError(f"unmatched '{ch}' in formula '{source}'")
        if ch in _OPEN_TO_CLOSE:
            sub, i = _parse_group(s, i + 1, _OPEN_TO_CLOSE[ch], source)
            mult, i = _read_number(s, i, source)
            for el, n in sub.items():
                counts[el] = counts.get(el, 0) + mult * n
        elif ch.isupper():
            j = i + 1
            while j < len(s) and s[j].islower():
                j += 1
            element = s[i:j]
            if element not in ELEMENTS:
                raise FormulaError(f"unknown element '{element}' in formula '{source}'")
            number, i = _read_number(s, j, source)
            counts[element] = counts.get(element, 0) + number
        else:
            raise FormulaError(f"unexpected character '{ch}' in formula '{source}'")
    if closer is not None:
        raise FormulaError(f"missing closing '{closer}' in formula '{source}'")
    return counts, i


def _split_segments(s: str, source: str) -> List[Tuple[int, str]]:
    """Split a formula on top-level hydrate dots, keeping leading multipliers."""
    parts: List[str] = []
    depth = 0
    start = 0
    for i, ch in enumerate(s):
        if ch in _OPEN_TO_CLOSE:
            depth += 1
        elif ch in _CLOSERS:
            depth -= 1
            if depth < 0:
                raise FormulaError(f"unmatched '{ch}' in formula '{source}'")
        elif ch in _HYDRATE_DOTS and depth == 0:
            parts.append(s[start:i])
            start = i + 1
    parts.append(s[start:])

    segments: List[Tuple[int, str]] = []
    for part in parts:
        m = _LEADING_COEFF_RE.match(part)
        if m:
            rest = m.group(2)
            if not rest:
                raise FormulaError(f"dangling multiplier in formula '{source}'")
            segments.append((int(m.group(1)), rest))
        else:
            if not part:
                raise FormulaError(f"empty hydrate segment in formula '{source}'")
            segments.append((1, part))
    return segments


def parse_formula(text: str) -> Tuple[Dict[str, int], int]:
    """Parse a formula into (element counts, net charge)."""
    if not isinstance(text, str) or not text.strip():
        raise FormulaError("formula must be a non-empty string")
    s = "".join(text.split())
    charge = 0
    if "^" in s:
        if s.count("^") > 1:
            raise FormulaError(f"formula '{text}' contains more than one '^'")
        s, charge_text = s.split("^", 1)
        charge = _parse_charge(charge_text, text)
    if not s:
        raise FormulaError(f"formula '{text}' has no chemical part")
    counts: Dict[str, int] = {}
    for coeff, segment in _split_segments(s, text):
        sub, _ = _parse_group(segment, 0, None, text)
        for el, n in sub.items():
            counts[el] = counts.get(el, 0) + coeff * n
    if not counts:
        raise FormulaError(f"formula '{text}' contains no elements")
    return counts, charge


# ---------------------------------------------------------------------------
# Exact linear algebra over the rationals
# ---------------------------------------------------------------------------


def _lcm(a: int, b: int) -> int:
    return a // gcd(a, b) * b


def _rref(matrix: List[List[Fraction]]) -> Tuple[List[List[Fraction]], List[int]]:
    """Reduced row-echelon form via exact rational Gaussian elimination."""
    rows = len(matrix)
    cols = len(matrix[0]) if rows else 0
    m = [list(row) for row in matrix]
    pivots: List[int] = []
    r = 0
    for c in range(cols):
        if r == rows:
            break
        pivot = next((i for i in range(r, rows) if m[i][c] != 0), None)
        if pivot is None:
            continue
        m[r], m[pivot] = m[pivot], m[r]
        inv = 1 / m[r][c]
        m[r] = [x * inv for x in m[r]]
        for i in range(rows):
            if i != r and m[i][c] != 0:
                factor = m[i][c]
                m[i] = [a - factor * b for a, b in zip(m[i], m[r])]
        pivots.append(c)
        r += 1
    return m, pivots


def null_space(matrix: List[List[Fraction]]) -> List[List[Fraction]]:
    """Exact basis of the null space of a rational matrix."""
    if not matrix:
        return []
    cols = len(matrix[0])
    rref, pivots = _rref(matrix)
    free = [c for c in range(cols) if c not in pivots]
    basis: List[List[Fraction]] = []
    for f in free:
        vec = [Fraction(0)] * cols
        vec[f] = Fraction(1)
        for row, pc in enumerate(pivots):
            vec[pc] = -rref[row][f]
        basis.append(vec)
    return basis


# ---------------------------------------------------------------------------
# Balancing and verification
# ---------------------------------------------------------------------------


@dataclass
class ParsedSpecies:
    formula: str
    counts: Dict[str, int]
    charge: int


@dataclass
class BalanceResult:
    reactants: List[ParsedSpecies]
    products: List[ParsedSpecies]
    reactant_coefficients: List[int]
    product_coefficients: List[int]
    elements: List[str]


def parse_all(formulas: Sequence[str]) -> List[ParsedSpecies]:
    return [
        ParsedSpecies(formula=f, counts=counts, charge=charge)
        for f in formulas
        for counts, charge in [parse_formula(f)]
    ]


def ordered_elements(species: Sequence[ParsedSpecies]) -> List[str]:
    elements: List[str] = []
    seen = set()
    for sp in species:
        for el in sp.counts:
            if el not in seen:
                seen.add(el)
                elements.append(el)
    return elements


def conservation_matrix(
    reactants: Sequence[ParsedSpecies], products: Sequence[ParsedSpecies], elements: Sequence[str]
) -> List[List[Fraction]]:
    """Element-conservation matrix with a trailing charge row.

    Reactant columns are positive, product columns negative, so null-space
    vectors of this matrix are coefficient sets conserving every element and
    the total charge.
    """
    species = list(reactants) + list(products)
    n_reactants = len(reactants)
    matrix: List[List[Fraction]] = []
    for el in elements:
        matrix.append([Fraction(sp.counts.get(el, 0)) for sp in species])
    matrix.append([Fraction(sp.charge) for sp in species])
    for row in matrix:
        for j in range(n_reactants, len(species)):
            row[j] = -row[j]
    return matrix


def balance(reactant_formulas: Sequence[str], product_formulas: Sequence[str]) -> BalanceResult:
    """Balance an equation exactly; raise BalanceError when impossible."""
    reactants = parse_all(reactant_formulas)
    products = parse_all(product_formulas)
    elements = ordered_elements(list(reactants) + list(products))
    matrix = conservation_matrix(reactants, products, elements)
    basis = null_space(matrix)
    dimension = len(basis)

    if dimension == 0:
        raise BalanceError(
            "no_solution",
            "The conservation matrix only admits the trivial all-zero solution: "
            "the equation cannot be balanced.",
            null_space_dimension=0,
        )
    if dimension > 1:
        raise BalanceError(
            "underdetermined",
            f"The null space of the conservation matrix has dimension {dimension}: "
            "infinitely many independent coefficient sets exist, so no unique "
            "balancing can be returned.",
            null_space_dimension=dimension,
        )

    vector = basis[0]
    denom = reduce(_lcm, (x.denominator for x in vector), 1)
    ints = [int(x * denom) for x in vector]
    common = reduce(gcd, (abs(x) for x in ints))
    ints = [x // common for x in ints]
    if all(x < 0 for x in ints):
        ints = [-x for x in ints]

    if any(x == 0 for x in ints):
        raise BalanceError(
            "not_all_positive",
            "The one-dimensional solution assigns a zero coefficient to at least "
            "one species, so no all-positive balancing exists.",
            null_space_dimension=1,
            raw_coefficients=ints,
        )
    if not all(x > 0 for x in ints):
        raise BalanceError(
            "not_all_positive",
            "The one-dimensional solution has mixed signs: at least one species "
            "would have to be moved to the other side, so no all-positive "
            "balancing exists.",
            null_space_dimension=1,
            raw_coefficients=ints,
        )

    n = len(reactants)
    return BalanceResult(
        reactants=reactants,
        products=products,
        reactant_coefficients=ints[:n],
        product_coefficients=ints[n:],
        elements=elements,
    )


def verify(
    reactants: Sequence[ParsedSpecies],
    products: Sequence[ParsedSpecies],
    reactant_coefficients: Sequence[int],
    product_coefficients: Sequence[int],
    elements: Sequence[str],
) -> dict:
    """Per-element and charge bookkeeping for a candidate coefficient set."""
    element_checks = []
    balanced = True
    for el in elements:
        left = sum(c * sp.counts.get(el, 0) for c, sp in zip(reactant_coefficients, reactants))
        right = sum(c * sp.counts.get(el, 0) for c, sp in zip(product_coefficients, products))
        ok = left == right
        balanced = balanced and ok
        element_checks.append(
            {"element": el, "reactants": left, "products": right, "balanced": ok}
        )
    left_q = sum(c * sp.charge for c, sp in zip(reactant_coefficients, reactants))
    right_q = sum(c * sp.charge for c, sp in zip(product_coefficients, products))
    charge_ok = left_q == right_q
    balanced = balanced and charge_ok
    return {
        "elements": element_checks,
        "charge": {"reactants": left_q, "products": right_q, "balanced": charge_ok},
        "balanced": balanced,
    }
