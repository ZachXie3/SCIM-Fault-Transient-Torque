"""
SCIM Fault Transient Torque — main entry point.

Loads motor parameters from ``input.jsonc`` (JSON with comments), runs both
the quick analytical calculation and the dq model time-domain simulation,
compares their results, and saves CSVs and overlay plots to ``data/``.

Usage
-----
    python run.py

All outputs are written to the ``data/`` directory (gitignored).
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from scim_calc.config import load_jsonc, normalize_params
from scim_calc.circuit import derive_slip
from scim_calc.quick import run_quick_calc
from scim_calc.dq import run_dq_simulation
from scim_calc.compare import compute_comparison, print_comparison

INPUT_FILE = Path(__file__).parent / "input.jsonc"
DATA_DIR = Path(__file__).parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)


def main():
    # ---- Load and normalize motor parameters ----
    print(f"Loading: {INPUT_FILE}")
    cfg = load_jsonc(str(INPUT_FILE))
    p = normalize_params(cfg)
    slip = cfg["slip"]

    if p["power_provided"]:
        slip = derive_slip(
            p["V_ph"], p["omega_s"], p["Rs"], p["Rr"],
            p["Xs"], p["Xr"], p["Xm"], p["pole_pairs"], p["P_mech"],
        )
        print(f"Derived slip from nameplate power: {slip:.6g} pu")
    else:
        print(f"Using slip from config: {slip}")

    p["slip"] = slip

    # ---- Quick calculation ----
    print("\n--- Quick Calculation ---")
    q = run_quick_calc(p, slip=slip)
    _save_quick(q)
    _print_quick_results(q)

    # ---- DQ model simulation ----
    print("\n--- DQ Model Simulation ---")
    d = run_dq_simulation(p, slip=slip)
    _save_dq(d)
    _print_dq_results(d)

    # ---- Comparison & export ----
    print("\n--- Comparison ---")
    cmp = compute_comparison(q, d)
    print_comparison(cmp)
    _save_comparison(cmp)

    # ---- Overlay plot ----
    _plot_overlay(q, d)

    print("\nDone. Results saved to:", DATA_DIR.resolve())


# ---------------------------------------------------------------------------
# Internal helpers — save CSV, print results, plot overlay
# ---------------------------------------------------------------------------

def _save_quick(q):
    csv = DATA_DIR / "quick_scim_short_circuit.csv"
    T_pct = 100.0 * q["T_fault"] / q["T_nom"]
    T_env_pct = 100.0 * q["T_env"] / q["T_nom"]
    np.savetxt(
        csv,
        np.column_stack([q["t"], q["T_fault"], T_pct, q["T_env"], T_env_pct]),
        delimiter=",",
        header=(
            "time_s,T_fault_Nm,T_fault_over_T_nom_percent,"
            "T_envelope_Nm,T_envelope_over_T_nom_percent"
        ),
        comments="",
    )
    print(f"  Saved: {csv.name}")


def _save_dq(d):
    csv = DATA_DIR / "dq_scim_short_circuit.csv"
    T_pct = 100.0 * d["T_fault"] / d["T_nom"]
    np.savetxt(
        csv,
        np.column_stack([d["t"], d["T_fault"], T_pct, d["omega_m_trace"], d["Is_abs"]]),
        delimiter=",",
        header="time_s,T_fault_Nm,T_fault_over_T_nom_percent,omega_m_rad_per_s,abs_i_s_peak_A",
        comments="",
    )
    print(f"  Saved: {csv.name}")


def _save_comparison(cmp):
    csv = DATA_DIR / "comparison_summary.csv"
    with open(csv, "w") as f:
        f.write("metric,value\n")
        for k, v in cmp.items():
            if v is not None:
                f.write(f"{k},{v}\n")
    print(f"  Saved: {csv.name}")


def _print_quick_results(q):
    print(f"  T_nom           = {q['T_nom']:.4g} Nm")
    print(f"  Peak |T|/T_nom  = {100 * np.max(np.abs(q['T_fault'])) / q['T_nom']:.2f} %")
    print(f"  Envelope init   = {100 * q['T_env0'] / q['T_nom']:.2f} %")
    print(f"  I_sc_ac0        = {q['I_sc_ac0']:.4g} A RMS")
    print(f"  Xs_sc           = {q['Xs_sc']:.4g} ohm,  Xr_sc = {q['Xr_sc']:.4g} ohm")
    print(f"  Tr0             = {q['Tr0']:.4g} s,  Tr_sc = {q['Tr_sc']:.4g} s,  Ts_dc = {q['Ts_dc']:.4g} s")


def _print_dq_results(d):
    print(f"  T_nom           = {d['T_nom']:.4g} Nm")
    print(f"  Peak |T|/T_nom  = {100 * np.max(np.abs(d['T_fault'])) / d['T_nom']:.2f} %")
    print(f"  Speed           = {d['omega_m_trace'][0]:.4g} rad/s")
    print(f"  Initial Is peak = {abs(d['Is0']):.4g} A")


def _plot_overlay(q, d):
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(q["t"], 100 * q["T_fault"] / q["T_nom"],
            label="Quick calculation", lw=1.2, ls="--")
    ax.plot(d["t"], 100 * d["T_fault"] / d["T_nom"],
            label="dq model", lw=1.2)
    ax.axhline(100, color="black", lw=0.8, ls=":", label="Nominal torque")
    ax.set(
        xlabel="Time (s)",
        ylabel="T / T_nom (%)",
        title="SCIM Short-Circuit Transient Torque — Quick vs dq Model",
    )
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    png = DATA_DIR / "scim_short_circuit_overlay.png"
    fig.savefig(png, dpi=160)
    print(f"  Saved: {png.name}")


if __name__ == "__main__":
    main()
