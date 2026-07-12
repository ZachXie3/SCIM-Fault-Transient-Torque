"""
Test the quick analytical calculation module.

Compares full torque traces, envelope, and peak values against saved
regression baselines and expected physics-based characteristics.
"""

import numpy as np
import pytest

from scim_calc.quick import run_quick_calc


def test_quick_tnom_matches_expected(motor_params, expected_quick):
    """Nominal torque at t=0 should match the saved baseline exactly."""
    q = run_quick_calc(motor_params, slip=motor_params["slip"])
    assert q["T_nom"] == pytest.approx(expected_quick["T_fault_Nm"][0], rel=1e-10)


def test_quick_traces_match(motor_params, expected_quick):
    """The full torque trace should match the regression baseline."""
    q = run_quick_calc(motor_params, slip=motor_params["slip"])
    assert np.allclose(q["T_fault"], expected_quick["T_fault_Nm"], rtol=1e-10)


def test_quick_env_match(motor_params, expected_quick):
    """The envelope trace should match the regression baseline."""
    q = run_quick_calc(motor_params, slip=motor_params["slip"])
    assert np.allclose(q["T_env"], expected_quick["T_env_Nm"], rtol=1e-10)


def test_quick_peak_torque(motor_params):
    """Peak absolute torque should be ~376% of nominal (test-motor characteristic)."""
    q = run_quick_calc(motor_params, slip=motor_params["slip"])
    T_pct = 100.0 * np.max(np.abs(q["T_fault"])) / q["T_nom"]
    assert T_pct == pytest.approx(376.35, rel=1e-2)


def test_quick_negative_peak(motor_params):
    """Negative peak should equal positive peak in magnitude (symmetric waveform)."""
    q = run_quick_calc(motor_params, slip=motor_params["slip"])
    neg_pct = 100.0 * np.min(q["T_fault"]) / q["T_nom"]
    assert neg_pct == pytest.approx(-376.35, rel=1e-2)
