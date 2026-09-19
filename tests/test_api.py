"""API tests using FastAPI's TestClient (requires fastapi + httpx)."""

import unittest

from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


class TestBalanceEndpoint(unittest.TestCase):
    def test_balance_success(self):
        resp = client.post(
            "/api/v1/balance",
            json={"reactants": ["H2", "O2"], "products": ["H2O"]},
        )
        self.assertEqual(resp.status_code, 200)
        body = resp.json()
        self.assertEqual(body["coefficients"], {"reactants": [2, 1], "products": [2]})
        self.assertEqual(body["equation"], "2 H2 + O2 -> 2 H2O")
        self.assertEqual(body["null_space_dimension"], 1)
        self.assertTrue(body["verification"]["balanced"])
        by_element = {e["element"]: e for e in body["verification"]["elements"]}
        self.assertEqual(by_element["H"]["reactants"], 4)
        self.assertEqual(by_element["H"]["products"], 4)
        self.assertEqual(by_element["O"]["reactants"], 2)
        self.assertEqual(by_element["O"]["products"], 2)
        self.assertEqual(body["verification"]["charge"]["reactants"], 0)

    def test_balance_ionic_equation(self):
        resp = client.post(
            "/api/v1/balance",
            json={
                "reactants": ["MnO4^-", "H^+", "Fe^2+"],
                "products": ["Mn^2+", "Fe^3+", "H2O"],
            },
        )
        self.assertEqual(resp.status_code, 200)
        body = resp.json()
        self.assertEqual(body["coefficients"]["reactants"], [1, 8, 5])
        self.assertEqual(body["coefficients"]["products"], [1, 5, 4])
        self.assertEqual(body["verification"]["charge"]["reactants"], 17)
        self.assertEqual(body["verification"]["charge"]["products"], 17)

    def test_balance_hydrate(self):
        resp = client.post(
            "/api/v1/balance",
            json={"reactants": ["CuSO4·5H2O"], "products": ["CuSO4", "H2O"]},
        )
        self.assertEqual(resp.status_code, 200)
        self.assertEqual(
            resp.json()["coefficients"], {"reactants": [1], "products": [1, 5]}
        )

    def test_balance_no_solution_422(self):
        resp = client.post(
            "/api/v1/balance",
            json={"reactants": ["H2", "O2"], "products": ["CO2"]},
        )
        self.assertEqual(resp.status_code, 422)
        detail = resp.json()["detail"]
        self.assertEqual(detail["error"], "no_solution")
        self.assertEqual(detail["null_space_dimension"], 0)

    def test_balance_underdetermined_422(self):
        resp = client.post(
            "/api/v1/balance",
            json={"reactants": ["C", "O2"], "products": ["CO", "CO2"]},
        )
        self.assertEqual(resp.status_code, 422)
        detail = resp.json()["detail"]
        self.assertEqual(detail["error"], "underdetermined")
        self.assertEqual(detail["null_space_dimension"], 2)

    def test_balance_not_all_positive_422(self):
        resp = client.post(
            "/api/v1/balance",
            json={"reactants": ["H2", "H2O2"], "products": ["O2"]},
        )
        self.assertEqual(resp.status_code, 422)
        self.assertEqual(resp.json()["detail"]["error"], "not_all_positive")

    def test_balance_invalid_formula_422(self):
        resp = client.post(
            "/api/v1/balance",
            json={"reactants": ["H2", "Xx9"], "products": ["H2O"]},
        )
        self.assertEqual(resp.status_code, 422)
        self.assertEqual(resp.json()["detail"]["error"], "invalid_formula")

    def test_balance_empty_list_422(self):
        resp = client.post("/api/v1/balance", json={"reactants": [], "products": ["H2O"]})
        self.assertEqual(resp.status_code, 422)


class TestValidateEndpoint(unittest.TestCase):
    def test_validate_correct(self):
        resp = client.post(
            "/api/v1/validate",
            json={
                "reactants": ["H2", "O2"],
                "products": ["H2O"],
                "coefficients": {"reactants": [2, 1], "products": [2]},
            },
        )
        self.assertEqual(resp.status_code, 200)
        self.assertTrue(resp.json()["balanced"])

    def test_validate_incorrect(self):
        resp = client.post(
            "/api/v1/validate",
            json={
                "reactants": ["H2", "O2"],
                "products": ["H2O"],
                "coefficients": {"reactants": [1, 1], "products": [1]},
            },
        )
        self.assertEqual(resp.status_code, 200)
        body = resp.json()
        self.assertFalse(body["balanced"])
        unbalanced = {e["element"] for e in body["verification"]["elements"] if not e["balanced"]}
        self.assertEqual(unbalanced, {"O"})  # H: 2 = 2 balances; O: 2 != 1 does not

    def test_validate_charge_imbalance(self):
        resp = client.post(
            "/api/v1/validate",
            json={
                "reactants": ["Fe^2+"],
                "products": ["Fe^3+"],
                "coefficients": {"reactants": [1], "products": [1]},
            },
        )
        self.assertEqual(resp.status_code, 200)
        body = resp.json()
        self.assertFalse(body["balanced"])
        self.assertFalse(body["verification"]["charge"]["balanced"])

    def test_validate_length_mismatch_422(self):
        resp = client.post(
            "/api/v1/validate",
            json={
                "reactants": ["H2", "O2"],
                "products": ["H2O"],
                "coefficients": {"reactants": [2], "products": [2]},
            },
        )
        self.assertEqual(resp.status_code, 422)

    def test_validate_nonpositive_coefficient_422(self):
        resp = client.post(
            "/api/v1/validate",
            json={
                "reactants": ["H2", "O2"],
                "products": ["H2O"],
                "coefficients": {"reactants": [2, 0], "products": [2]},
            },
        )
        self.assertEqual(resp.status_code, 422)


class TestMeta(unittest.TestCase):
    def test_health(self):
        resp = client.get("/health")
        self.assertEqual(resp.status_code, 200)
        self.assertEqual(resp.json(), {"status": "ok"})


if __name__ == "__main__":
    unittest.main()
