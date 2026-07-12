"""
Test A: Angle invariance.

For the ideal balanced three-phase short circuit in a linear, symmetric
machine, the torque waveform should be exactly invariant to the fault
inception angle theta0. This test runs the dq simulation at multiple
angles and verifies that the torque traces are identical.
"""

import numpy as np
import pytest

from scim_calc.dq import run_dq_simulation


def test_torque_angle_invariance_0_90(motor_params):
    """Torque at theta0=0 deg and theta0=90 deg should be identical."""
    p = dict(motor_params)
    p["INITIAL_VOLTAGE_ANGLE_DEG"] = 0.0
    r0 = run_dq_simulation(p, slip=p["slip"])

    p["INITIAL_VOLTAGE_ANGLE_DEG"] = 90.0
    r90 = run_dq_simulation(p, slip=p["slip"])

    max_diff = np.max(np.abs(r0["T_fault"] - r90["T_fault"]))
    assert max_diff < 1e-8, f"Max torque difference at 90 deg: {max_diff:.2e}"


def test_torque_angle_invariance_0_45(motor_params):
    """Torque at theta0=0 deg and theta0=45 deg should be identical."""
    p = dict(motor_params)
    p["INITIAL_VOLTAGE_ANGLE_DEG"] = 0.0
    r0 = run_dq_simulation(p, slip=p["slip"])

    p["INITIAL_VOLTAGE_ANGLE_DEG"] = 45.0
    r45 = run_dq_simulation(p, slip=p["slip"])

    # Mathematically exact invariance; tiny numerical differences arise from
    # floating-point roundoff in the linear solver with rotated RHS vectors.
    max_diff = np.max(np.abs(r0["T_fault"] - r45["T_fault"]))
    assert max_diff < 1e-8, f"Max torque difference at 45 deg: {max_diff:.2e}"


def test_torque_angle_invariance_0_180(motor_params):
    """Torque at theta0=0 deg and theta0=180 deg should be identical
    (180-degree periodicity)."""
    p = dict(motor_params)
    p["INITIAL_VOLTAGE_ANGLE_DEG"] = 0.0
    r0 = run_dq_simulation(p, slip=p["slip"])

    p["INITIAL_VOLTAGE_ANGLE_DEG"] = 180.0
    r180 = run_dq_simulation(p, slip=p["slip"])

    max_diff = np.max(np.abs(r0["T_fault"] - r180["T_fault"]))
    assert max_diff < 1e-8, f"Max torque difference at 180 deg: {max_diff:.2e}"


def test_torque_angle_invariance_multiple(motor_params):
    """Torque at several angles should all be identical."""
    p = dict(motor_params)
    p["INITIAL_VOLTAGE_ANGLE_DEG"] = 0.0
    r0 = run_dq_simulation(p, slip=p["slip"])

    for angle in [15.0, 30.0, 60.0, 120.0, 135.0, 170.0]:
        p["INITIAL_VOLTAGE_ANGLE_DEG"] = angle
        r = run_dq_simulation(p, slip=p["slip"])
        max_diff = np.max(np.abs(r0["T_fault"] - r["T_fault"]))
        assert max_diff < 1e-8, f"Max torque difference at angle={angle}: {max_diff:.2e}"
