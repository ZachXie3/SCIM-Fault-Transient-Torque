"""
Stationary-frame (alpha-beta / Clarke) space-vector SCIM short-circuit simulation.

Technically an alpha-beta model (stationary frame, omega=0), but referred to
as "dq" in the code following industry convention (two-axis state-space
machine models are universally called dq models in the literature).

Integrates the flux-linkage differential equations using RK4 after a
balanced three-phase bolted terminal short circuit (v_s = 0 for t >= 0).

Compared to the quick method, this naturally captures:
  - Full multi-component torque spectrum (not just a single frequency).
  - DC offset via flux continuity at fault inception.
  - Point-on-wave dependence (configurable inception angle).
  - Optional speed dynamics.
"""

import numpy as np

from .circuit import currents_from_flux, torque_from_flux_current, solve_prefault


def run_dq_simulation(p, slip=None):
    """Run the dq-model time-domain simulation.

    Parameters
    ----------
    p : dict — normalized motor parameters (from config.normalize_params).
    slip : float, optional — overrides p["slip"] if provided.

    Returns
    -------
    dict with keys:
        t                : time array (s)
        T_fault          : signed torque waveform (N·m)
        T_nom            : pre-fault nominal torque (N·m)
        omega_m_trace    : mechanical speed at each step (rad/s)
        Is_abs           : peak stator current magnitude at each step (A)
        psi_s0           : initial stator flux linkage (complex, peak)
        Is0              : initial stator current (complex, peak)
    """
    omega_s = p["omega_s"]
    pole_pairs = p["pole_pairs"]
    slip = slip if slip is not None else p["slip"]
    angle_deg = p.get("INITIAL_VOLTAGE_ANGLE_DEG", 0.0)
    use_speed_dynamics = p.get("USE_SPEED_DYNAMICS", False)
    j_total = p.get("J", 0.0)  # total inertia for speed dynamics

    # ---- Pre-fault initial condition ----
    pf = solve_prefault(
        p["V_ph"], omega_s, slip, p["Rs"], p["Rr"],
        p["Ls"], p["Lr"], p["Lm"], angle_deg,
    )
    Is0 = pf["Is0"]
    Ir0 = pf["Ir0"]
    psi_s0 = pf["psi_s0"]
    psi_r0 = pf["psi_r0"]

    # Mechanical speed at t=0 (pole-pairs division applied here)
    omega_m0 = pf["omega_m0"] / pole_pairs

    # Nominal torque and sign convention (positive = motor)
    T_nom_raw = torque_from_flux_current(psi_s0, Is0, pole_pairs)
    torque_sign = 1.0 if T_nom_raw >= 0.0 else -1.0
    T_nom = abs(T_nom_raw)
    T_load_raw = T_nom_raw  # load torque equals pre-fault motor torque

    # State vector: [Re(psi_s), Im(psi_s), Re(psi_r), Im(psi_r), omega_m]
    y = np.array([psi_s0.real, psi_s0.imag, psi_r0.real, psi_r0.imag, omega_m0])

    # Local copies for speed inside nested functions
    Ls, Lr, Lm, Delta = p["Ls"], p["Lr"], p["Lm"], p["Delta"]
    Rs, Rr = p["Rs"], p["Rr"]

    # ---- Right-hand side of the ODE system ----
    def rhs(state):
        ps = state[0] + 1j * state[1]          # psi_s (complex)
        pr = state[2] + 1j * state[3]          # psi_r (complex)
        om = state[4] if use_speed_dynamics else omega_m0
        ore = pole_pairs * om                   # rotor electrical speed
        i_s, i_r = currents_from_flux(ps, pr, Ls, Lr, Lm, Delta)

        # Stator: v_s=0 after short circuit -> dpsi_s/dt = -Rs * i_s
        dpsi_s = -Rs * i_s
        # Rotor (squirrel-cage, v_r=0): dpsi_r/dt = j*omega_re*psi_r - Rr*i_r
        dpsi_r = 1j * ore * pr - Rr * i_r

        dom = 0.0
        if use_speed_dynamics:
            Te = torque_from_flux_current(ps, i_s, pole_pairs)
            dom = (Te - T_load_raw) / j_total

        return np.array([dpsi_s.real, dpsi_s.imag, dpsi_r.real, dpsi_r.imag, dom])

    # ---- Classical RK4 integrator ----
    def rk4_step(state, dt):
        k1 = rhs(state)
        k2 = rhs(state + 0.5 * dt * k1)
        k3 = rhs(state + 0.5 * dt * k2)
        k4 = rhs(state + dt * k3)
        return state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    # ---- Time integration loop ----
    t = np.linspace(0.0, p["T_END"], p["N_POINTS"])
    dt = t[1] - t[0]
    T_fault = np.zeros(p["N_POINTS"])
    omega_m_trace = np.zeros(p["N_POINTS"])
    Is_abs = np.zeros(p["N_POINTS"])

    state = y.copy()
    for k in range(p["N_POINTS"]):
        ps = state[0] + 1j * state[1]
        pr = state[2] + 1j * state[3]
        i_s, i_r = currents_from_flux(ps, pr, Ls, Lr, Lm, Delta)
        Te_raw = torque_from_flux_current(ps, i_s, pole_pairs)
        T_fault[k] = torque_sign * Te_raw
        omega_m_trace[k] = state[4]
        Is_abs[k] = abs(i_s)
        if k < p["N_POINTS"] - 1:
            state = rk4_step(state, dt)

    return {
        "t": t,
        "T_fault": T_fault,
        "T_nom": T_nom,
        "omega_m_trace": omega_m_trace,
        "Is_abs": Is_abs,
        "psi_s0": psi_s0,
        "Is0": Is0,
    }
