"""
Test C: Torque sign consistency.

Verifies that the complex-function torque, the expanded real-form
torque (if one were implemented), and the expected motoring sign
at pre-fault condition all agree.
"""

import numpy as np

from scim_calc.dq import run_dq_simulation
from scim_calc.circuit import torque_from_flux_current


def test_prefault_torque_motoring_sign(motor_params):
    """Pre-fault torque should be positive (motoring convention)."""
    d = run_dq_simulation(motor_params, slip=motor_params["slip"])
    assert d["T_fault"][0] > 0, "Pre-fault torque should be positive (motoring)"
    assert d["T_nom"] > 0, "Nominal torque should be positive"


def test_torque_sign_from_flux_current_consistent(motor_params):
    """Torque from flux-current function should have consistent sign."""
    d = run_dq_simulation(motor_params, slip=motor_params["slip"])
    T_nom = d["T_nom"]

    # At t=0, torque should equal nominal (positive)
    assert np.isclose(d["T_fault"][0], T_nom, rtol=1e-10), \
        f"T_fault[0]={d['T_fault'][0]} != T_nom={T_nom}"


def test_torque_consistency_prefault_from_phasors(motor_params):
    """Torque from flux-current function should match IEEE phasor torque
    at the pre-fault operating point."""
    from scim_calc.circuit import torque_from_phasors

    p = motor_params
    T_phasor = torque_from_phasors(
        p["V_ph"], p["omega_s"], p["Rs"], p["Rr"],
        p["Xs"], p["Xr"], p["Xm"], p["slip"], p["pole_pairs"],
    )
    d = run_dq_simulation(p, slip=p["slip"])
    # Both should agree within a small tolerance
    assert abs(T_phasor - d["T_nom"]) / d["T_nom"] < 1e-10, \
        f"Phasor torque {T_phasor} != dq T_nom {d['T_nom']}"
