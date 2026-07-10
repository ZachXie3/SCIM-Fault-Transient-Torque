# scim_short_circuit_quick_calculation_fixed.py
"""
Quick SCIM three-phase terminal short-circuit torque estimate.

Replace the USER PARAMETERS block with your motor data. Parameters are per-phase,
stator-referred, and reactances are at frequency f.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# =========================
# USER PARAMETERS
# =========================
Rs = 0.006       # ohm/phase, stator resistance
Rr = 0.020       # ohm/phase, rotor resistance referred to stator; do NOT divide by slip
Xs = 0.340       # ohm/phase, stator leakage reactance at f
Xr = 0.330       # ohm/phase, rotor leakage reactance at f
Xm = 10.600      # ohm/phase, magnetizing reactance at f
J = 1.0          # kg*m^2, not used in quick constant-speed estimate
slip = 0.020     # pu pre-fault slip
V_LL = 230.0     # V RMS line-line pre-fault voltage
f = 60.0         # Hz
pole_pairs = 3   # p

T_END = 0.200    # seconds
N_POINTS = 5000

# Optional first-cycle asymmetry sensitivity term.
# 0.0 = no extra dc/asymmetry term; 1.0 = strong first-cycle sensitivity term.
DC_OFFSET_FACTOR = 0.0
FAULT_ANGLE_DEG = 0.0

# Torque frequency mode: "rotor" uses (1-slip)*omega_s; "line" uses omega_s.
TORQUE_FREQUENCY_MODE = "line"

# =========================
# CALCULATION
# =========================
if not (0.0 < slip < 1.0):
    raise ValueError("slip must be between 0 and 1 pu for motor operation.")
if min(Rs, Rr, Xs, Xr, Xm, V_LL, f, pole_pairs) <= 0:
    raise ValueError("Rs, Rr, Xs, Xr, Xm, V_LL, f, and pole_pairs must be positive.")

j = 1j
omega_s = 2.0 * np.pi * f
V_phi = V_LL / np.sqrt(3.0)  # RMS phase voltage

# Initial steady-state equivalent circuit
Zr = Rr / slip + j * Xr
Zm = j * Xm
Zp = (Zm * Zr) / (Zm + Zr)
Is0 = V_phi / (Rs + j * Xs + Zp)
Eag0 = V_phi - Is0 * (Rs + j * Xs)
Ir0 = Eag0 / (Rr / slip + j * Xr)

omega_syn_m = omega_s / pole_pairs
T_nom = 3.0 * abs(Ir0) ** 2 * (Rr / slip) / omega_syn_m
T_base = abs(T_nom)
if T_base <= 0:
    raise ValueError("Computed nominal torque is zero; check input parameters and slip.")

# Short-circuit reactances and time constants
Xs_sc = Xs + (Xm * Xr) / (Xm + Xr)
Xr_sc = Xr + (Xm * Xs) / (Xm + Xs)
Tr0 = (Xm + Xr) / (omega_s * Rr)
Tr_sc = Xr_sc / (omega_s * Rr)
Ts_dc = Xs_sc / (omega_s * Rs)
T_decay_rotor = Tr_sc / 2.0
T_decay_mixed = (Ts_dc * Tr_sc) / (Ts_dc + Tr_sc)
T_decay = T_decay_mixed if DC_OFFSET_FACTOR > 0 else T_decay_rotor

# Initial short-circuit current and torque envelope
I_sc_ac0 = abs(Eag0) / abs(Rs + j * Xs_sc)  # RMS symmetrical estimate
Psi_m0_rms = abs(Eag0) / omega_s
T_env0 = 3.0 * pole_pairs * Psi_m0_rms * I_sc_ac0
if T_env0 <= 0:
    raise ValueError("Computed torque envelope is zero; check input parameters.")

ratio_for_phi = np.clip(T_nom / T_env0, -1.0, 1.0)
phi = np.arccos(ratio_for_phi)

if TORQUE_FREQUENCY_MODE.lower() == "line":
    omega_T = omega_s
else:
    omega_T = (1.0 - slip) * omega_s

fault_angle = np.deg2rad(FAULT_ANGLE_DEG)
t = np.linspace(0.0, T_END, N_POINTS)

T_ac = T_env0 * np.exp(-t / T_decay_rotor) * np.cos(omega_T * t + phi)
T_dc_sensitivity = (
    DC_OFFSET_FACTOR
    * T_env0
    * np.exp(-t / T_decay_mixed)
    * np.cos(omega_T * t + fault_angle)
)
T_fault = T_ac + T_dc_sensitivity
T_env = T_env0 * np.exp(-t / T_decay)

T_percent = 100.0 * T_fault / T_base
T_env_percent = 100.0 * T_env / T_base

# =========================
# OUTPUT
# =========================
outdir = Path(__file__).resolve().parent
csv_path = outdir / "quick_scim_short_circuit_fixed.csv"
png_path = outdir / "quick_scim_short_circuit_fixed.png"

np.savetxt(
    csv_path,
    np.column_stack([t, T_fault, T_percent, T_env, T_env_percent]),
    delimiter=",",
    header="time_s,T_fault_Nm,T_fault_over_T_nom_percent,T_envelope_Nm,T_envelope_over_T_nom_percent",
    comments="",
)

plt.figure(figsize=(9, 5))
plt.plot(t, T_percent, label="Quick signed torque estimate")
plt.plot(t, T_env_percent, "--", label="Positive envelope")
plt.plot(t, -T_env_percent, "--", color="gray", linewidth=0.9, label="Negative envelope")
plt.axhline(100.0, color="black", linewidth=0.8, linestyle=":", label="Pre-fault nominal torque")
plt.xlabel("Time (s)")
plt.ylabel("T_fault / T_nominal (%)")
plt.title("SCIM 3-Phase Terminal Short Circuit — Quick Calculation")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(png_path, dpi=160)
plt.show()

print("Quick SCIM short-circuit calculation complete")
print(f"T_nom = {T_nom:.6g} N*m")
print(f"Xs_sc = {Xs_sc:.6g} ohm, Xr_sc = {Xr_sc:.6g} ohm")
print(f"Tr0 = {Tr0:.6g} s, Tr_sc = {Tr_sc:.6g} s, Ts_dc = {Ts_dc:.6g} s")
print(f"I_sc_ac0 = {I_sc_ac0:.6g} A RMS")
print(f"T_env0/T_nom = {T_env0 / T_base:.6g} pu")
print(f"Peak plotted ratio = {np.max(np.abs(T_percent)):.6g} %")
print(f"Saved CSV: {csv_path}")
print(f"Saved plot: {png_path}")
