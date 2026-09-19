"""FastAPI application exposing the chemical equation balancer."""

from __future__ import annotations

from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse

from .chemistry import BalanceError, FormulaError, balance, verify, parse_all, ordered_elements
from .schemas import (
    BalanceRequest,
    BalanceResponse,
    ChargeCheck,
    Coefficients,
    ElementCheck,
    ErrorResponse,
    SpeciesCoefficient,
    ValidateRequest,
    ValidateResponse,
    Verification,
)

app = FastAPI(
    title="Chemical Equation Balancer",
    version="1.0.0",
    description=(
        "Exact chemical equation balancing. The service builds an element-conservation "
        "matrix (including a charge row) from the submitted formulas, computes its null "
        "space with rational Gaussian elimination, and returns GCD-reduced minimal "
        "positive integer coefficients when — and only when — the null space is "
        "one-dimensional and yields an all-positive solution. Otherwise it answers "
        "HTTP 422 with the null-space dimension or the reason no all-positive "
        "solution exists."
    ),
)


class UnprocessableError(Exception):
    """Domain error reported to the client as HTTP 422."""

    def __init__(self, payload: dict):
        super().__init__(payload.get("message"))
        self.payload = payload


@app.exception_handler(UnprocessableError)
async def unprocessable_handler(_: Request, exc: UnprocessableError) -> JSONResponse:
    return JSONResponse(status_code=422, content={"detail": exc.payload})


def _fail(error: str, message: str, **extra) -> None:
    raise UnprocessableError({"error": error, "message": message, **extra})


def _build_verification(verification: dict) -> Verification:
    return Verification(
        elements=[ElementCheck(**e) for e in verification["elements"]],
        charge=ChargeCheck(**verification["charge"]),
        balanced=verification["balanced"],
    )


def _format_side(formulas, coefficients) -> str:
    parts = [f if c == 1 else f"{c} {f}" for f, c in zip(formulas, coefficients)]
    return " + ".join(parts)


@app.get("/health", tags=["meta"])
def health() -> dict:
    return {"status": "ok"}


@app.post(
    "/api/v1/balance",
    response_model=BalanceResponse,
    responses={422: {"model": ErrorResponse, "description": "Cannot be balanced uniquely"}},
    tags=["balancing"],
)
def balance_endpoint(request: BalanceRequest) -> BalanceResponse:
    """Balance a chemical equation exactly.

    Returns the minimal positive integer coefficients together with a
    per-element and charge balance sheet. Responds 422 when the null space of
    the conservation matrix is not one-dimensional or cannot be scaled to
    all-positive coefficients.
    """
    try:
        result = balance(request.reactants, request.products)
    except FormulaError as exc:
        _fail("invalid_formula", str(exc))
    except BalanceError as exc:
        _fail(exc.code, exc.message, **exc.context)

    verification = verify(
        result.reactants,
        result.products,
        result.reactant_coefficients,
        result.product_coefficients,
        result.elements,
    )
    equation = (
        f"{_format_side(request.reactants, result.reactant_coefficients)} -> "
        f"{_format_side(request.products, result.product_coefficients)}"
    )
    return BalanceResponse(
        equation=equation,
        reactants=[
            SpeciesCoefficient(formula=f, coefficient=c)
            for f, c in zip(request.reactants, result.reactant_coefficients)
        ],
        products=[
            SpeciesCoefficient(formula=f, coefficient=c)
            for f, c in zip(request.products, result.product_coefficients)
        ],
        coefficients=Coefficients(
            reactants=result.reactant_coefficients, products=result.product_coefficients
        ),
        null_space_dimension=1,
        verification=_build_verification(verification),
    )


@app.post(
    "/api/v1/validate",
    response_model=ValidateResponse,
    responses={422: {"model": ErrorResponse, "description": "Invalid input"}},
    tags=["balancing"],
)
def validate_endpoint(request: ValidateRequest) -> ValidateResponse:
    """Check a proposed set of coefficients against the conservation laws.

    Always answers 200 with a per-element and charge balance sheet;
    ``balanced`` is true only when every element and the total charge are
    conserved.
    """
    try:
        reactants = parse_all(request.reactants)
        products = parse_all(request.products)
    except FormulaError as exc:
        _fail("invalid_formula", str(exc))

    elements = ordered_elements(reactants + products)
    verification = verify(
        reactants,
        products,
        request.coefficients.reactants,
        request.coefficients.products,
        elements,
    )
    return ValidateResponse(
        balanced=verification["balanced"],
        verification=_build_verification(verification),
    )
