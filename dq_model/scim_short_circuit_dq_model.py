# scim_short_circuit_dq_model_fixed.py
"""
Stationary-frame dq/space-vector SCIM three-phase terminal short-circuit model.

Replace the USER PARAMETERS block with your motor data. Parameters are per-phase,
stator-referred, and reactances are at frequency f.
"""

import json
import re
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

input_path = Path(__file__).resolve().parent.parent / "input.jsonc"
text = input_path.read_text(encoding="utf-8")
text = re.sub(r"//.*", "", text)
text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
cfg = json.loads(text)

Rs = cfg["Rs"]
Rr = cfg["Rr"]
J = cfg["J"]
slip = cfg["slip"]
V_LL = cfg["V_LL"]
f = cfg["f"]
pole_pairs = cfg["pole_pairs"]

# Accept either reactances (Xs, Xr, Xm) or inductances (Ls, Lr, Lm)
omega_s = 2.0 * np.pi * f
if "Ls" in cfg and "Lr" in cfg and "Lm" in cfg:
    Lls = cfg["Ls"]
    Llr = cfg["Lr"]
    Lm = cfg["Lm"]
else:
    Xs = cfg["Xs"]
    Xr = cfg["Xr"]
    Xm = cfg["Xm"]
    Lls = Xs / omega_s
    Llr = Xr / omega_s
    Lm = Xm / omega_s

if min(Lls, Llr, Lm) <= 0:
    raise ValueError("Ls, Lr, and Lm must be positive.")

T_END = cfg["T_END"]
N_POINTS = cfg["N_POINTS"]
INITIAL_VOLTAGE_ANGLE_DEG = cfg["INITIAL_VOLTAGE_ANGLE_DEG"]
USE_SPEED_DYNAMICS = cfg["USE_SPEED_DYNAMICS"]

# =========================
# HELPER FUNCTIONS
# =========================
def validate_inputs():
    if not (0.0 < slip < 1.0):
        raise ValueError("slip must be between 0 and 1 pu for motor operation.")
    if min(Rs, Rr, V_LL, f, pole_pairs) <= 0:
        raise ValueError("Rs, Rr, V_LL, f, and pole_pairs must be positive.")
    if USE_SPEED_DYNAMICS and J <= 0:
        raise ValueError("J must be positive when speed dynamics are enabled.")


def currents_from_flux(psi_s, psi_r, Ls, Lr, Lm, Delta):
    i_s = (Lr * psi_s - Lm * psi_r) / Delta
    i_r = (-Lm * psi_s + Ls * psi_r) / Delta
    return i_s, i_r


def torque_from_flux_current(psi_s, i_s, p):
    # Te = (3/2)*p*(psi_ds*i_qs - psi_qs*i_ds)
    # For psi = psi_ds + j psi_qs and i = i_ds + j i_qs:
    return -1.5 * p * np.imag(psi_s * np.conj(i_s))


# =========================
# INITIALIZATION
# =========================
validate_inputs()

Ls = Lls + Lm
Lr = Llr + Lm
Delta = Ls * Lr - Lm ** 2
if Delta <= 0:
    raise ValueError("Invalid inductance set: Delta = Ls*Lr - Lm^2 must be positive.")

if cfg.get("CONNECTION", "wye").lower() == "delta":
    V_ph = V_LL
else:
    V_ph = V_LL / np.sqrt(3.0)
angle0 = np.deg2rad(INITIAL_VOLTAGE_ANGLE_DEG)
Vs0 = np.sqrt(2.0) * V_ph * np.exp(1j * angle0)

# Synchronous steady-state phasor equations on peak basis.
A = np.array(
    [
        [Rs + 1j * omega_s * Ls, 1j * omega_s * Lm],
        [1j * slip * omega_s * Lm, Rr + 1j * slip * omega_s * Lr],
    ],
    dtype=complex,
)
b = np.array([Vs0, 0.0 + 0.0j], dtype=complex)
Is0, Ir0 = np.linalg.solve(A, b)
psi_s0 = Ls * Is0 + Lm * Ir0
psi_r0 = Lm * Is0 + Lr * Ir0

omega_m0 = (1.0 - slip) * omega_s / pole_pairs
omega_re0 = pole_pairs * omega_m0

T_nom_raw = torque_from_flux_current(psi_s0, Is0, pole_pairs)
torque_sign = 1.0 if T_nom_raw >= 0.0 else -1.0
T_nom = abs(T_nom_raw)
if T_nom <= 0:
    raise ValueError("Computed nominal torque is zero; check input parameters and slip.")
T_load_raw = T_nom_raw  # same sign convention as raw torque for speed dynamics

# State vector: [Re(psi_s), Im(psi_s), Re(psi_r), Im(psi_r), omega_m]
y = np.array([psi_s0.real, psi_s0.imag, psi_r0.real, psi_r0.imag, omega_m0], dtype=float)


def rhs(state):
    psi_s = state[0] + 1j * state[1]
    psi_r = state[2] + 1j * state[3]
    omega_m = state[4] if USE_SPEED_DYNAMICS else omega_m0
    omega_re = pole_pairs * omega_m

    i_s, i_r = currents_from_flux(psi_s, psi_r, Ls, Lr, Lm, Delta)

    # Terminal short circuit: v_s = 0
    dpsi_s = -Rs * i_s
    dpsi_r = 1j * omega_re * psi_r - Rr * i_r

    Te_raw = torque_from_flux_current(psi_s, i_s, pole_pairs)
    if USE_SPEED_DYNAMICS:
        domega_m = (Te_raw - T_load_raw) / J
    else:
        domega_m = 0.0

    return np.array([dpsi_s.real, dpsi_s.imag, dpsi_r.real, dpsi_r.imag, domega_m], dtype=float)


def rk4_step(state, dt):
    k1 = rhs(state)
    k2 = rhs(state + 0.5 * dt * k1)
    k3 = rhs(state + 0.5 * dt * k2)
    k4 = rhs(state + dt * k3)
    return state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

# =========================
# TIME INTEGRATION
# =========================
t = np.linspace(0.0, T_END, N_POINTS)
dt = t[1] - t[0]
T_fault = np.zeros_like(t)
omega_m_trace = np.zeros_like(t)
Is_abs = np.zeros_like(t)

state = y.copy()
for k, tk in enumerate(t):
    psi_s = state[0] + 1j * state[1]
    psi_r = state[2] + 1j * state[3]
    i_s, i_r = currents_from_flux(psi_s, psi_r, Ls, Lr, Lm, Delta)
    Te_raw = torque_from_flux_current(psi_s, i_s, pole_pairs)
    T_fault[k] = torque_sign * Te_raw
    omega_m_trace[k] = state[4]
    Is_abs[k] = abs(i_s)
    if k < len(t) - 1:
        state = rk4_step(state, dt)

T_percent = 100.0 * T_fault / T_nom

# =========================
# OUTPUT
# =========================
outdir = Path(__file__).resolve().parent.parent / "data"
csv_path = outdir / "dq_scim_short_circuit.csv"
png_path = outdir / "dq_scim_short_circuit.png"

np.savetxt(
    csv_path,
    np.column_stack([t, T_fault, T_percent, omega_m_trace, Is_abs]),
    delimiter=",",
    header="time_s,T_fault_Nm,T_fault_over_T_nom_percent,omega_m_rad_per_s,abs_i_s_peak_A",
    comments="",
)

plt.figure(figsize=(9, 5))
plt.plot(t, T_percent, label="dq model torque")
plt.axhline(100.0, color="black", linewidth=0.8, linestyle=":", label="Pre-fault nominal torque")
plt.xlabel("Time (s)")
plt.ylabel("T_fault / T_nominal (%)")
plt.title("SCIM 3-Phase Terminal Short Circuit — dq Model")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(png_path, dpi=160)
plt.show()

print("dq SCIM short-circuit calculation complete")
print(f"T_nom = {T_nom:.6g} N*m")
print(f"Initial stator current peak magnitude = {abs(Is0):.6g} A")
print(f"Initial rotor current peak magnitude = {abs(Ir0):.6g} A")
print(f"Initial mechanical speed = {omega_m0:.6g} rad/s")
print(f"Peak plotted ratio = {np.max(np.abs(T_percent)):.6g} %")
print(f"Saved CSV: {csv_path}")
print(f"Saved plot: {png_path}")
