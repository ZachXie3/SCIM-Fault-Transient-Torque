"""
Test the shared equivalent-circuit functions: flux–current inversion,
torque computation, steady-state solver, phasor torque, and slip derivation.
"""

import numpy as np
import pytest

from scim_calc.circuit import (
    currents_from_flux,
    torque_from_flux_current,
    solve_prefault,
    torque_from_phasors,
    derive_slip,
)


def test_currents_from_flux(motor_params):
    """Verify flux–current inversion: real-valued inputs produce real currents."""
    p = motor_params
    i_s, i_r = currents_from_flux(
        1.0 + 0j, 0.5 + 0j, p["Ls"], p["Lr"], p["Lm"], p["Delta"],
    )
    assert np.isreal(i_s)
    assert np.isreal(i_r)


def test_torque_from_flux_current():
    """Verify torque sign convention: aligned → zero, quadrature → positive."""
    # In-phase flux and current → zero torque
    t = torque_from_flux_current(1.0 + 0j, 1.0 + 0j, 1)
    assert abs(t) < 1e-15
    # Quadrature → positive motor torque per sign convention
    t = torque_from_flux_current(1.0 + 0j, 0.0 + 1j, 1)
    assert t == pytest.approx(1.5)


def test_solve_prefault(motor_params):
    """Pre-fault steady-state solver should produce ~3989 Nm for the test motor."""
    p = motor_params
    pf = solve_prefault(
        p["V_ph"], p["omega_s"], p["slip"],
        p["Rs"], p["Rr"], p["Ls"], p["Lr"], p["Lm"],
    )
    T = abs(torque_from_flux_current(pf["psi_s0"], pf["Is0"], p["pole_pairs"]))
    assert T == pytest.approx(3989, rel=1e-2)


def test_torque_from_phasors(motor_params):
    """Phasor-based torque should match the flux–current result (~3989 Nm)."""
    p = motor_params
    T = torque_from_phasors(
        p["V_ph"], p["omega_s"], p["Rs"], p["Rr"],
        p["Xs"], p["Xr"], p["Xm"], p["slip"], p["pole_pairs"],
    )
    assert T == pytest.approx(3989, rel=1e-2)


def test_derive_slip(motor_params):
    """Derived slip should give torque that balances the nameplate power."""
    p = motor_params
    slip = derive_slip(
        p["V_ph"], p["omega_s"], p["Rs"], p["Rr"],
        p["Xs"], p["Xr"], p["Xm"], p["pole_pairs"], p["P_mech"],
    )
    assert slip == pytest.approx(0.00779, abs=1e-4)
    # Verify power balance: T * omega_m ≈ P_mech
    T = torque_from_phasors(
        p["V_ph"], p["omega_s"], p["Rs"], p["Rr"],
        p["Xs"], p["Xr"], p["Xm"], slip, p["pole_pairs"],
    )
    omega_m = (1.0 - slip) * p["omega_s"] / p["pole_pairs"]
    P_calc = T * omega_m
    assert P_calc == pytest.approx(p["P_mech"], rel=1e-3)
