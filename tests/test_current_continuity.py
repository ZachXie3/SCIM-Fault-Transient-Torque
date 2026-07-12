"""
Test B: Current continuity at fault inception.

Verifies that stator and rotor currents immediately before and immediately
after the fault are identical (currents are continuous when flux linkages
and the inductance matrix are continuous).
"""

import numpy as np

from scim_calc.circuit import currents_from_flux, torque_from_flux_current, solve_prefault


def test_stator_current_continuous(motor_params):
    """Stator current computed from pre-fault phasors (t=0-) should equal
    current recovered from continuous flux linkages (t=0+)."""
    p = motor_params
    pf = solve_prefault(
        p["V_ph"], p["omega_s"], p["slip"],
        p["Rs"], p["Rr"], p["Ls"], p["Lr"], p["Lm"],
        angle_deg=0.0,
    )
    # Currents from the phasor solution (pre-fault, t=0-)
    Is0_phasor = pf["Is0"]
    Ir0_phasor = pf["Ir0"]

    # Currents recovered from continuous flux linkages (post-fault, t=0+)
    i_s_plus, i_r_plus = currents_from_flux(
        pf["psi_s0"], pf["psi_r0"],
        p["Ls"], p["Lr"], p["Lm"], p["Delta"],
    )

    assert np.isclose(Is0_phasor, i_s_plus, rtol=1e-12, atol=1e-12), \
        f"Stator current discontinuity: phasor={Is0_phasor}, flux_recovered={i_s_plus}"
    assert np.isclose(Ir0_phasor, i_r_plus, rtol=1e-12, atol=1e-12), \
        f"Rotor current discontinuity: phasor={Ir0_phasor}, flux_recovered={i_r_plus}"


def test_torque_continuous_across_fault(motor_params):
    """Torque should also be continuous (since it depends on continuous
    currents and fluxes)."""
    p = motor_params
    pf = solve_prefault(
        p["V_ph"], p["omega_s"], p["slip"],
        p["Rs"], p["Rr"], p["Ls"], p["Lr"], p["Lm"],
        angle_deg=0.0,
    )

    # Torque at t=0- from phasor currents
    T_pre = torque_from_flux_current(pf["psi_s0"], pf["Is0"], p["pole_pairs"])

    # Torque at t=0+ from continuous flux and recovered current
    i_s_plus, _ = currents_from_flux(
        pf["psi_s0"], pf["psi_r0"],
        p["Ls"], p["Lr"], p["Lm"], p["Delta"],
    )
    T_post = torque_from_flux_current(pf["psi_s0"], i_s_plus, p["pole_pairs"])

    assert np.isclose(T_pre, T_post, rtol=1e-12, atol=1e-12), \
        f"Torque discontinuity: pre={T_pre}, post={T_post}"
