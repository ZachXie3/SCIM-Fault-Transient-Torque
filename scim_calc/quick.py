"""
Quick (closed-form analytical) SCIM short-circuit torque calculation.

Computes the transient electromagnetic torque after a balanced three-phase
terminal short circuit using approximate closed-form formulas rather than
time-domain integration.

Key assumptions
---------------
- Single-frequency torque oscillation (configurable: line or rotor frequency).
- Single exponential decay envelope.
- No speed dynamics (constant speed).
- DC offset handled approximately via an optional additive term.
- RMS phasor approximation for initial flux linkage.

Use this for screening and order-of-magnitude checks.  Use the dq model
when point-on-wave shape or detailed waveform is needed.
"""

import numpy as np

from .circuit import torque_from_phasors


def run_quick_calc(p, slip=None):
    """Run the quick analytical torque calculation.

    Parameters
    ----------
    p : dict — normalized motor parameters (from config.normalize_params).
    slip : float, optional — overrides p["slip"] if provided.

    Returns
    -------
    dict with keys:
        t                        : time array (s)
        T_fault                  : signed torque waveform (N·m)
        T_env                    : positive envelope (N·m)
        T_nom                    : pre-fault nominal torque (N·m)
        T_env0                   : initial envelope amplitude (N·m)
        envelope_initial_pu      : T_env0 / T_nom
        I_sc_ac0                 : initial symmetrical RMS current (A)
        phi                      : phase angle of torque oscillation (rad)
        omega_T                  : torque oscillation frequency (rad/s)
        T_decay                  : envelope decay time constant (s)
        T_decay_rotor            : rotor-component decay (s)
        T_decay_mixed            : stator–rotor mixed decay (s)
        Tr_sc, Ts_dc, Tr0        : derived time constants (s)
        Xs_sc, Xr_sc             : short-circuit reactances (ohm)
    """
    omega_s = p["omega_s"]
    pole_pairs = p["pole_pairs"]
    slip = slip if slip is not None else p["slip"]

    # Nominal torque from IEEE equivalent circuit (phasor-based)
    T_nom = abs(torque_from_phasors(
        p["V_ph"], omega_s, p["Rs"], p["Rr"], p["Xs"], p["Xr"], p["Xm"],
        slip, pole_pairs,
    ))
    T_base = T_nom

    # ---- Short-circuit reactances (blocked-rotor equivalents) ----
    # Stator seen from stator with rotor shorted
    Xs_sc = p["Xs"] + (p["Xm"] * p["Xr"]) / (p["Xm"] + p["Xr"])
    # Rotor seen from rotor with stator shorted
    Xr_sc = p["Xr"] + (p["Xm"] * p["Xs"]) / (p["Xm"] + p["Xs"])

    # ---- Time constants ----
    # Rotor open-circuit: flux decay when rotor is open
    Tr0 = (p["Xm"] + p["Xr"]) / (omega_s * p["Rr"])
    # Rotor short-circuit: AC component decay
    Tr_sc = Xr_sc / (omega_s * p["Rr"])
    # Stator DC component decay
    Ts_dc = Xs_sc / (omega_s * p["Rs"])
    # Torque envelope decay dominated by rotor flux
    T_decay_rotor = Tr_sc / 2.0
    # Mixed decay when DC offset is included
    T_decay_mixed = (Ts_dc * Tr_sc) / (Ts_dc + Tr_sc)
    dc_factor = p.get("DC_OFFSET_FACTOR", 0.0)
    T_decay = T_decay_mixed if dc_factor > 0 else T_decay_rotor

    # ---- Pre-fault equivalent circuit (same as circuit.py but with reactances) ----
    j = 1j
    Zr_branch = p["Rr"] / slip + j * p["Xr"]
    Zm = j * p["Xm"]
    Zp = (Zm * Zr_branch) / (Zm + Zr_branch)
    Is0 = p["V_ph"] / (p["Rs"] + j * p["Xs"] + Zp)
    Eag0 = p["V_ph"] - Is0 * (p["Rs"] + j * p["Xs"])  # air-gap voltage

    # ---- Initial short-circuit current and torque envelope ----
    # Symmetrical RMS current estimate
    I_sc_ac0 = abs(Eag0) / abs(p["Rs"] + j * Xs_sc)
    # RMS air-gap flux linkage
    Psi_m0_rms = abs(Eag0) / omega_s
    # Initial torque envelope: T_env0 = 3 * p * |Psi_m0| * I_sc_ac0
    T_env0 = 3.0 * pole_pairs * Psi_m0_rms * I_sc_ac0

    # Phase angle so the waveform starts near T_nom
    ratio_for_phi = np.clip(T_base / T_env0, -1.0, 1.0)
    phi = np.arccos(ratio_for_phi)

    # Torque oscillation frequency
    freq_mode = p.get("TORQUE_FREQUENCY_MODE", "rotor").lower()
    if freq_mode == "line":
        omega_T = omega_s
    else:
        omega_T = (1.0 - slip) * omega_s  # rotor electrical frequency

    fault_angle = np.deg2rad(p.get("FAULT_ANGLE_DEG", 0.0))

    # ---- Construct waveform ----
    t = np.linspace(0.0, p["T_END"], p["N_POINTS"])

    # Main AC component
    T_ac = T_env0 * np.exp(-t / T_decay_rotor) * np.cos(omega_T * t + phi)

    # Optional DC-offset sensitivity term (same frequency, different decay)
    T_dc = (
        dc_factor
        * T_env0
        * np.exp(-t / T_decay_mixed)
        * np.cos(omega_T * t + fault_angle)
    )
    T_fault = T_ac + T_dc

    # Positive envelope (magnitude only)
    T_env = T_env0 * np.exp(-t / T_decay)

    return {
        "t": t,
        "T_fault": T_fault,
        "T_env": T_env,
        "T_nom": T_nom,
        "T_env0": T_env0,
        "envelope_initial_pu": T_env0 / T_base,
        "I_sc_ac0": I_sc_ac0,
        "phi": phi,
        "omega_T": omega_T,
        "T_decay": T_decay,
        "T_decay_rotor": T_decay_rotor,
        "T_decay_mixed": T_decay_mixed,
        "Tr_sc": Tr_sc,
        "Ts_dc": Ts_dc,
        "Tr0": Tr0,
        "Xs_sc": Xs_sc,
        "Xr_sc": Xr_sc,
    }
