"""
Verify that config loading and normalization produce correct derived values
from the test fixture input.jsonc.
"""

import numpy as np
import pytest

from scim_calc.config import load_jsonc, normalize_params


def test_load_jsonc(fixtures_dir):
    """Raw fixture file should parse with correct parameter values."""
    cfg = load_jsonc(fixtures_dir / "input.jsonc")
    assert cfg["Rs"] == 0.069966
    assert cfg["Rr"] == 0.072979
    assert cfg["Xs"] == 1.272654
    assert cfg["V_LL"] == 4000.0
    assert cfg["f"] == 60.0
    assert cfg["poles"] == 2
    assert cfg["HP"] == 2000


def test_normalize_params(fixtures_dir):
    """Derived quantities (omega_s, inductances, V_ph, inertia, power) should
    be computed correctly from the raw config."""
    cfg = load_jsonc(fixtures_dir / "input.jsonc")
    p = normalize_params(cfg)
    omega_s = 2.0 * np.pi * 60.0

    assert p["omega_s"] == pytest.approx(omega_s)
    assert p["pole_pairs"] == 1
    assert p["Lls"] == pytest.approx(1.272654 / omega_s)
    assert p["Llr"] == pytest.approx(0.986283 / omega_s)
    assert p["Lm"] == pytest.approx(42.29672 / omega_s)
    assert p["V_ph"] == pytest.approx(4000.0 / np.sqrt(3.0))
    assert p["J"] == pytest.approx(40.670 / 4.0)
    assert p["P_mech"] == 2000 * 746.0
    assert p["power_provided"] is True
