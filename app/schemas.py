"""Request/response models for the balancing API (Pydantic v2)."""

from __future__ import annotations

from typing import List, Optional

from pydantic import BaseModel, ConfigDict, Field, model_validator


class BalanceRequest(BaseModel):
    """Equation to balance, given as lists of formulas."""

    reactants: List[str] = Field(
        ..., min_length=1, description="Reactant formulas, e.g. [\"KMnO4\", \"HCl\"]"
    )
    products: List[str] = Field(
        ..., min_length=1, description="Product formulas, e.g. [\"KCl\", \"MnCl2\", \"H2O\", \"Cl2\"]"
    )

    model_config = ConfigDict(
        json_schema_extra={
            "examples": [
                {
                    "reactants": ["KMnO4", "HCl"],
                    "products": ["KCl", "MnCl2", "H2O", "Cl2"],
                }
            ]
        }
    )


class Coefficients(BaseModel):
    """Proposed stoichiometric coefficients, aligned with the formula lists."""

    reactants: List[int] = Field(..., min_length=1)
    products: List[int] = Field(..., min_length=1)


class ValidateRequest(BalanceRequest):
    """Equation plus a candidate coefficient set to verify."""

    coefficients: Coefficients

    @model_validator(mode="after")
    def _check_coefficients(self) -> "ValidateRequest":
        if len(self.coefficients.reactants) != len(self.reactants):
            raise ValueError(
                "coefficients.reactants must have exactly one entry per reactant "
                f"({len(self.reactants)} expected, got {len(self.coefficients.reactants)})"
            )
        if len(self.coefficients.products) != len(self.products):
            raise ValueError(
                "coefficients.products must have exactly one entry per product "
                f"({len(self.products)} expected, got {len(self.coefficients.products)})"
            )
        for value in self.coefficients.reactants + self.coefficients.products:
            if value <= 0:
                raise ValueError("coefficients must be positive integers")
        return self

    model_config = ConfigDict(
        json_schema_extra={
            "examples": [
                {
                    "reactants": ["H2", "O2"],
                    "products": ["H2O"],
                    "coefficients": {"reactants": [2, 1], "products": [2]},
                }
            ]
        }
    )


class SpeciesCoefficient(BaseModel):
    formula: str
    coefficient: int


class ElementCheck(BaseModel):
    element: str
    reactants: int
    products: int
    balanced: bool


class ChargeCheck(BaseModel):
    reactants: int
    products: int
    balanced: bool


class Verification(BaseModel):
    """Per-element and charge conservation bookkeeping."""

    elements: List[ElementCheck]
    charge: ChargeCheck
    balanced: bool


class BalanceResponse(BaseModel):
    equation: str
    reactants: List[SpeciesCoefficient]
    products: List[SpeciesCoefficient]
    coefficients: Coefficients
    null_space_dimension: int
    verification: Verification


class ValidateResponse(BaseModel):
    balanced: bool
    verification: Verification


class ErrorDetail(BaseModel):
    error: str
    message: str
    null_space_dimension: Optional[int] = None
    raw_coefficients: Optional[List[int]] = None


class ErrorResponse(BaseModel):
    detail: ErrorDetail
