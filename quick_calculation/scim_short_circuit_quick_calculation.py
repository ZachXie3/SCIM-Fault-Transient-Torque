# scim_short_circuit_quick_calculation_fixed.py
"""
Quick SCIM three-phase terminal short-circuit torque estimate.

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
slip = cfg["slip"]
V_LL = cfg["V_LL"]
f = cfg["f"]
pole_pairs = cfg["poles"] // 2

# Accept either reactances (Xs, Xr, Xm) or inductances (Ls, Lr, Lm)
omega_s = 2.0 * np.pi * f
if "Ls" in cfg and "Lr" in cfg and "Lm" in cfg:
    Xs = omega_s * cfg["Ls"]
    Xr = omega_s * cfg["Lr"]
    Xm = omega_s * cfg["Lm"]
else:
    Xs = cfg["Xs"]
    Xr = cfg["Xr"]
    Xm = cfg["Xm"]

T_END = cfg["T_END"]
N_POINTS = cfg["N_POINTS"]
DC_OFFSET_FACTOR = cfg["DC_OFFSET_FACTOR"]
FAULT_ANGLE_DEG = cfg["FAULT_ANGLE_DEG"]
TORQUE_FREQUENCY_MODE = cfg["TORQUE_FREQUENCY_MODE"]

# =========================
# CALCULATION
# =========================
if cfg.get("CONNECTION", "wye").lower() == "delta":
    V_ph = V_LL
else:
    V_ph = V_LL / np.sqrt(3.0)

# Optional: derive slip from nameplate power when HP or kW is provided
if "HP" in cfg and "kW" in cfg:
    raise ValueError("Provide EITHER 'HP' or 'kW', not both.")
power_provided = False
P_mech_W = None
if "HP" in cfg:
    P_mech_W = cfg["HP"] * 746.0
    power_provided = True
elif "kW" in cfg:
    P_mech_W = cfg["kW"] * 1000.0
    power_provided = True

if power_provided:
    def _torque_from_slip(s):
        Zr = Rr / s + 1j * Xr
        Zm = 1j * Xm
        Zp = (Zm * Zr) / (Zm + Zr)
        Is = V_ph / (Rs + 1j * Xs + Zp)
        Eag = V_ph - Is * (Rs + 1j * Xs)
        Ir = Eag / (Rr / s + 1j * Xr)
        o_syn = omega_s / pole_pairs
        return 3.0 * abs(Ir) ** 2 * (Rr / s) / o_syn

    o_syn = omega_s / pole_pairs
    def _f(s):
        return _torque_from_slip(s) - P_mech_W / ((1.0 - s) * o_syn)

    s_low = 1e-4
    s_high = 0.1
    f_low, f_high = _f(s_low), _f(s_high)
    # If root not bracketed (torque may cross target twice due to breakdown peak),
    # widen upper bound until opposite sign is found
    while f_low * f_high > 0 and s_high < 0.8:
        s_high = min(s_high * 2, 0.8)
        f_high = _f(s_high)
    if f_low * f_high > 0:
        print("Warning: slip derivation did not bracket root; using provided slip.")
    else:
        for _ in range(100):
            s_mid = (s_low + s_high) / 2.0
            f_mid = _f(s_mid)
            if abs(f_mid) < 1e-10:
                break
            if f_low * f_mid <= 0:
                s_high = s_mid
            else:
                s_low = s_mid
                f_low = f_mid
        s_derived = (s_low + s_high) / 2.0
        print(f"Derived slip from nameplate power: {s_derived:.6g} pu (was {slip:.6g})")
        slip = s_derived

if not (0.0 < slip < 1.0):
    raise ValueError("slip must be between 0 and 1 pu for motor operation.")
if min(Rs, Rr, V_LL, f, pole_pairs) <= 0:
    raise ValueError("Rs, Rr, V_LL, f, and pole_pairs must be positive.")

j = 1j

# Initial steady-state equivalent circuit
Zr = Rr / slip + j * Xr
Zm = j * Xm
Zp = (Zm * Zr) / (Zm + Zr)
Is0 = V_ph / (Rs + j * Xs + Zp)
Eag0 = V_ph - Is0 * (Rs + j * Xs)
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
outdir = Path(__file__).resolve().parent.parent / "data"
csv_path = outdir / "quick_scim_short_circuit.csv"
png_path = outdir / "quick_scim_short_circuit.png"

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
print(f"T_nom (from equivalent circuit) = {T_nom:.6g} N*m")
if power_provided:
    omega_m_mech = (1.0 - slip) * omega_s / pole_pairs
    T_full_load = P_mech_W / omega_m_mech
    pct_diff = (T_nom - T_full_load) / T_full_load * 100.0
    print(f"Speed = {omega_m_mech:.4g} rad/s  ({omega_m_mech * 60 / (2*np.pi):.4g} RPM)")
    print(f"T_full_load (from nameplate power) = {T_full_load:.6g} N*m")
    print(f"T_nom vs T_full_load: {pct_diff:+.4g} % difference ({(T_nom/T_full_load - 1)*100:.4g} % error)")
print(f"Xs_sc = {Xs_sc:.6g} ohm, Xr_sc = {Xr_sc:.6g} ohm")
print(f"Tr0 = {Tr0:.6g} s, Tr_sc = {Tr_sc:.6g} s, Ts_dc = {Ts_dc:.6g} s")
print(f"I_sc_ac0 = {I_sc_ac0:.6g} A RMS")
print(f"T_env0/T_nom = {T_env0 / T_base:.6g} pu")
print(f"Peak plotted ratio = {np.max(np.abs(T_percent)):.6g} %")
print(f"Saved CSV: {csv_path}")
print(f"Saved plot: {png_path}")
