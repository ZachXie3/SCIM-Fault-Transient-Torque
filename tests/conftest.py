"""
Shared pytest fixtures for all test modules.

Provides session-scoped fixtures for:
  - cfg         : normalized motor parameters dict (from test fixtures).
  - motor_params: same as cfg but with derived slip applied.
  - expected_quick, expected_dq, expected_sweep: regression baseline data.
"""

import json
from pathlib import Path

import numpy as np
import pytest

from scim_calc.config import load_jsonc, normalize_params
from scim_calc.circuit import derive_slip

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture(scope="session")
def fixtures_dir():
    return FIXTURES


@pytest.fixture(scope="session")
def cfg():
    cfg = load_jsonc(FIXTURES / "input.jsonc")
    return normalize_params(cfg)


@pytest.fixture(scope="session")
def motor_params(cfg):
    p = dict(cfg)
    if p["power_provided"]:
        slip = derive_slip(
            p["V_ph"], p["omega_s"], p["Rs"], p["Rr"],
            p["Xs"], p["Xr"], p["Xm"], p["pole_pairs"], p["P_mech"],
        )
        p["slip"] = slip
    return p


@pytest.fixture(scope="session")
def expected_quick():
    return np.genfromtxt(
        FIXTURES / "expected_quick.csv",
        delimiter=",",
        names=True,
    )


@pytest.fixture(scope="session")
def expected_dq():
    return np.genfromtxt(
        FIXTURES / "expected_dq.csv",
        delimiter=",",
        names=True,
    )


@pytest.fixture(scope="session")
def expected_sweep():
    with open(FIXTURES / "expected_sweep.json") as f:
        return json.load(f)
