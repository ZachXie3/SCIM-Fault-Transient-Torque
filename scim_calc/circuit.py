"""
Steady-state equivalent circuit analysis shared by quick and dq methods.

Provides the mathematical core for both calculation approaches:
  - Flux–current relations (used in the dq state-space model)
  - Phasor-based torque computation (used in the quick method and slip derivation)
  - Pre-fault steady-state solver
  - Slip derivation from nameplate power via bisection
"""

import numpy as np


def currents_from_flux(psi_s, psi_r, Ls, Lr, Lm, Delta):
    """Recover stator and rotor current space vectors from flux linkages.

    Inverts the flux–current relation:
        [psi_s]   [Ls  Lm] [i_s]
        [psi_r] = [Lm  Lr] [i_r]

    Parameters
    ----------
    psi_s, psi_r : complex — stator/rotor flux linkage space vectors (peak).
    Ls, Lr, Lm   : float — self and mutual inductances (H).
    Delta        : float — Ls*Lr - Lm^2 (leakage determinant).

    Returns
    -------
    i_s, i_r : complex — current space vectors (peak).
    """
    i_s = (Lr * psi_s - Lm * psi_r) / Delta
    i_r = (-Lm * psi_s + Ls * psi_r) / Delta
    return i_s, i_r


def torque_from_flux_current(psi_s, i_s, pole_pairs):
    """Electromagnetic torque from flux and current space vectors.

    Uses the stationary-frame dq expression:
        Te = -(3/2) * p * Im{psi_s * conj(i_s)}

    The sign convention gives positive torque for motor operation.
    """
    return -1.5 * pole_pairs * np.imag(psi_s * np.conj(i_s))


def solve_prefault(V_ph, omega_s, slip, Rs, Rr, Ls, Lr, Lm, angle_deg=0.0):
    """Solve the pre-fault steady-state phasor equations in the stationary frame.

    The 2x2 complex system encodes stator and rotor voltage equations at
    synchronous speed (omega_s for stator, slip*omega_s for rotor).

    Parameters
    ----------
    V_ph       : float — RMS line-neutral voltage (V).
    omega_s    : float — electrical synchronous frequency (rad/s).
    slip       : float — per-unit slip.
    Rs, Rr     : float — stator/rotor resistance (ohm).
    Ls, Lr, Lm : float — inductances (H).
    angle_deg  : float — initial voltage space-vector angle (deg).

    Returns
    -------
    dict with Is0, Ir0 (complex phasors, peak), psi_s0, psi_r0 (complex
    flux linkages, peak), and omega_m0 (mechanical speed in rad/s, before
    dividing by pole_pairs).
    """
    angle0 = np.deg2rad(angle_deg)
    Vs0 = np.sqrt(2.0) * V_ph * np.exp(1j * angle0)   # peak phase voltage

    # Stator–rotor coupling matrix in the stationary frame
    A = np.array(
        [
            [Rs + 1j * omega_s * Ls, 1j * omega_s * Lm],
            [1j * slip * omega_s * Lm, Rr + 1j * slip * omega_s * Lr],
        ],
        dtype=complex,
    )
    b = np.array([Vs0, 0.0 + 0.0j], dtype=complex)
    Is0, Ir0 = np.linalg.solve(A, b)

    # Initial flux linkages
    psi_s0 = Ls * Is0 + Lm * Ir0
    psi_r0 = Lm * Is0 + Lr * Ir0

    # Mechanical speed before pole-pairs division
    omega_m0 = (1.0 - slip) * omega_s

    return {
        "Is0": Is0,
        "Ir0": Ir0,
        "psi_s0": psi_s0,
        "psi_r0": psi_r0,
        "omega_m0": omega_m0,
    }


def torque_from_phasors(V_ph, omega_s, Rs, Rr, Xs, Xr, Xm, slip, pole_pairs):
    """Electromagnetic torque from the IEEE per-phase equivalent circuit.

    Uses standard steady-state induction-machine formulas with reactances
    (at frequency f) rather than inductances.

    Reference: IEEE Std 112, equivalent circuit at slip s:
        T = 3 * |Ir|^2 * (Rr/s) / omega_syn_m

    Returns torque in N·m.
    """
    # Rotor branch at slip s: Z_r = Rr/s + jXr
    Zr = Rr / slip + 1j * Xr
    # Magnetizing branch
    Zm = 1j * Xm
    # Parallel combination
    Zp = (Zm * Zr) / (Zm + Zr)
    # Stator current
    Is = V_ph / (Rs + 1j * Xs + Zp)
    # Air-gap voltage
    Eag = V_ph - Is * (Rs + 1j * Xs)
    # Rotor current
    Ir = Eag / (Rr / slip + 1j * Xr)
    # Synchronous mechanical speed
    omega_syn_m = omega_s / pole_pairs
    # Torque from rotor copper loss
    return 3.0 * abs(Ir) ** 2 * (Rr / slip) / omega_syn_m


def derive_slip(V_ph, omega_s, Rs, Rr, Xs, Xr, Xm, pole_pairs, P_mech, slip_guess=0.005):
    """Derive operating slip from nameplate power via bisection.

    Solves T_motor(s) = T_load(s) where:
        T_load(s) = P_mech / [(1-s) * omega_syn_m]

    The search brackets [1e-4, 0.1] and widens up to 0.8 if needed.
    Convergence to 1e-10 absolute tolerance in at most 100 iterations.

    Returns the slip (per-unit) that satisfies the power balance.
    """
    o_syn = omega_s / pole_pairs

    def torque(s):
        return torque_from_phasors(V_ph, omega_s, Rs, Rr, Xs, Xr, Xm, s, pole_pairs)

    def f(s):
        # Residual: motor torque − load torque
        return torque(s) - P_mech / ((1.0 - s) * o_syn)

    s_low = 1e-4
    s_high = 0.1
    f_low, f_high = f(s_low), f(s_high)

    # Widen upper bound until f changes sign
    while f_low * f_high > 0 and s_high < 0.8:
        s_high = min(s_high * 2, 0.8)
        f_high = f(s_high)

    if f_low * f_high > 0:
        raise ValueError(
            "Could not bracket slip root. "
            "Check nameplate power and motor parameters."
        )

    # Standard bisection
    for _ in range(100):
        s_mid = (s_low + s_high) / 2.0
        f_mid = f(s_mid)
        if abs(f_mid) < 1e-10:
            return s_mid
        if f_low * f_mid <= 0:
            s_high = s_mid
        else:
            s_low = s_mid
            f_low = f_mid

    return (s_low + s_high) / 2.0
