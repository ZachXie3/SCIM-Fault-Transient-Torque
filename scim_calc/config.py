"""
Load and normalize motor parameters from a JSONC input file.

Reads the shared input.jsonc, strips JSON-style comments, and produces a flat
params dict with all derived quantities (inductances, phase voltage, inertia,
simulation options) ready for passing to calculation functions.
"""

import json
import re
from pathlib import Path

import numpy as np


def load_jsonc(path):
    """Load a JSON-with-comments file, return raw Python dict.
    Strips // line comments and /* block comments before parsing.
    """
    text = Path(path).read_text(encoding="utf-8")
    text = re.sub(r"//.*", "", text)
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
    return json.loads(text)


def normalize_params(cfg):
    """Normalize raw config dict into a flat params dict with derived quantities.

    The returned dict contains all original keys plus computed values:
    omega_s, pole_pairs, inductances (Lls, Llr, Lm, Ls, Lr, Delta),
    phase voltage V_ph, total inertia J, nameplate power P_mech,
    and simulation options with defaults.

    Parameters are stored in both reactance (Xs, Xr, Xm) and inductance
    form so that both quick (reactance-based) and dq (inductance-based)
    methods can use the same params dict.
    """
    p = dict(cfg)

    # Electrical angular frequency and pole-pair count
    p["omega_s"] = 2.0 * np.pi * cfg["f"]
    p["pole_pairs"] = cfg["poles"] // 2

    # Inductance / reactance dual representation
    if "Lls" in cfg and "Llr" in cfg and "Lm" in cfg:
        p["Xs"] = p["omega_s"] * cfg["Lls"]
        p["Xr"] = p["omega_s"] * cfg["Llr"]
        p["Xm"] = p["omega_s"] * cfg["Lm"]
        p["Lls"] = cfg["Lls"]
        p["Llr"] = cfg["Llr"]
        p["Lm"] = cfg["Lm"]
    else:
        p["Xs"] = cfg["Xs"]
        p["Xr"] = cfg["Xr"]
        p["Xm"] = cfg["Xm"]
        p["Lls"] = cfg["Xs"] / p["omega_s"]
        p["Llr"] = cfg["Xr"] / p["omega_s"]
        p["Lm"] = cfg["Xm"] / p["omega_s"]

    # Total self-inductances and leakage determinant
    p["Ls"] = p["Lls"] + p["Lm"]
    p["Lr"] = p["Llr"] + p["Lm"]
    p["Delta"] = p["Ls"] * p["Lr"] - p["Lm"] ** 2

    # Phase voltage (line-neutral RMS) from line-line value
    if cfg.get("CONNECTION", "wye").lower() == "delta":
        p["V_ph"] = cfg["V_LL"]
    else:
        p["V_ph"] = cfg["V_LL"] / np.sqrt(3.0)

    # Total rotating inertia (rotor + optional load)
    p["J"] = 0.0
    if "ROTOR_GD2" in cfg:
        p["J"] += cfg["ROTOR_GD2"] / 4.0       # GD^2 / 4 = J (kg·m^2)
    elif "ROTOR_WK2" in cfg:
        p["J"] += cfg["ROTOR_WK2"] * 0.042140   # WK^2 → J conversion (lb·ft^2 → kg·m^2)
    if "LOAD_GD2" in cfg:
        p["J"] += cfg["LOAD_GD2"] / 4.0
    elif "LOAD_WK2" in cfg:
        p["J"] += cfg["LOAD_WK2"] * 0.042140

    # Nameplate mechanical output power
    if "HP" in cfg and "kW" in cfg:
        raise ValueError("Provide EITHER HP or kW, not both.")
    if "HP" in cfg:
        p["P_mech"] = cfg["HP"] * 746.0          # 1 HP = 746 W
        p["power_provided"] = True
    elif "kW" in cfg:
        p["P_mech"] = cfg["kW"] * 1000.0
        p["power_provided"] = True
    else:
        p["P_mech"] = None
        p["power_provided"] = False

    # Simulation options with safe defaults
    p["T_END"] = cfg.get("T_END", 0.2)
    p["N_POINTS"] = cfg.get("N_POINTS", 5000)
    p["INITIAL_VOLTAGE_ANGLE_DEG"] = cfg.get("INITIAL_VOLTAGE_ANGLE_DEG", 0.0)
    p["USE_SPEED_DYNAMICS"] = cfg.get("USE_SPEED_DYNAMICS", False)
    p["DC_OFFSET_FACTOR"] = cfg.get("DC_OFFSET_FACTOR", 0.0)
    p["FAULT_ANGLE_DEG"] = cfg.get("FAULT_ANGLE_DEG", 0.0)
    p["TORQUE_FREQUENCY_MODE"] = cfg.get("TORQUE_FREQUENCY_MODE", "rotor")

    return p
