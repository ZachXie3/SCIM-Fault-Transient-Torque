"""
Test the dq model time-domain simulation module.

Compares torque traces against saved regression baselines and verifies
physics-based characteristics (asymmetry, constant speed when dynamics
are disabled).
"""

import numpy as np
import pytest

from scim_calc.dq import run_dq_simulation


def test_dq_tnom_matches_expected(motor_params, expected_dq):
    """Nominal torque at t=0 should match the saved baseline exactly."""
    d = run_dq_simulation(motor_params, slip=motor_params["slip"])
    assert d["T_nom"] == pytest.approx(expected_dq["T_fault_Nm"][0], rel=1e-10)


def test_dq_traces_match(motor_params, expected_dq):
    """The full torque trace should match the regression baseline."""
    d = run_dq_simulation(motor_params, slip=motor_params["slip"])
    assert np.allclose(d["T_fault"], expected_dq["T_fault_Nm"], rtol=1e-10)


def test_dq_speed_constant(motor_params):
    """With speed dynamics disabled, mechanical speed should remain constant."""
    d = run_dq_simulation(motor_params, slip=motor_params["slip"])
    assert np.allclose(d["omega_m_trace"], d["omega_m_trace"][0])


def test_dq_peak_abs_torque(motor_params):
    """Peak absolute torque should be ~417% of nominal (test-motor characteristic)."""
    d = run_dq_simulation(motor_params, slip=motor_params["slip"])
    T_pct = 100.0 * np.max(np.abs(d["T_fault"])) / d["T_nom"]
    assert T_pct == pytest.approx(416.84, rel=1e-2)


def test_dq_asymmetry(motor_params):
    """The dq model should produce asymmetric positive/negative peaks
    (negative peak larger in magnitude by >10% Tn)."""
    d = run_dq_simulation(motor_params, slip=motor_params["slip"])
    pos_pct = 100.0 * np.max(d["T_fault"]) / d["T_nom"]
    neg_pct = 100.0 * np.min(d["T_fault"]) / d["T_nom"]
    assert abs(neg_pct) > pos_pct
    assert (abs(neg_pct) - pos_pct) > 10.0
