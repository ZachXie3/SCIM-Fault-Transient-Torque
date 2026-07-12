"""
Compare quick calculation and dq model results.

Computes point-wise difference metrics, peak comparison, and time-segmented
RMS errors to quantify where and why the two methods diverge.
"""

import numpy as np


def compute_comparison(quick_res, dq_res):
    """Compute comparison metrics between quick and dq results.

    Parameters
    ----------
    quick_res : dict from run_quick_calc
    dq_res : dict from run_dq_simulation

    Returns
    -------
    dict with keys:
        T_nom, quick_pos_peak, quick_neg_peak,
        dq_pos_peak, dq_neg_peak,
        diff_rms, diff_rms_pct, diff_max, diff_max_pct,
        diff_early_rms_pct (t < 50 ms), diff_late_rms_pct (t > 100 ms)
    """
    T_nom = quick_res["T_nom"]
    Tq = quick_res["T_fault"]
    Td = dq_res["T_fault"]
    t = quick_res["t"]

    # If time vectors differ, interpolate dq onto the quick grid
    if len(Td) != len(t):
        Td = np.interp(t, dq_res["t"], Td)

    diff = Td - Tq
    abs_diff = np.abs(diff)

    return {
        "T_nom": T_nom,
        "quick_pos_peak": np.max(Tq),
        "quick_neg_peak": np.min(Tq),
        "dq_pos_peak": np.max(Td),
        "dq_neg_peak": np.min(Td),
        "diff_rms": np.sqrt(np.mean(diff ** 2)),
        "diff_rms_pct": 100.0 * np.sqrt(np.mean(diff ** 2)) / T_nom,
        "diff_max": np.max(abs_diff),
        "diff_max_pct": 100.0 * np.max(abs_diff) / T_nom,
        "diff_late_rms_pct": (
            100.0 * np.sqrt(np.mean(diff[t >= 0.1] ** 2)) / T_nom
            if np.any(t >= 0.1)
            else None
        ),
        "diff_early_rms_pct": (
            100.0 * np.sqrt(np.mean(diff[t <= 0.05] ** 2)) / T_nom
            if np.any(t <= 0.05)
            else None
        ),
    }


def print_comparison(cmp):
    """Pretty-print the comparison metrics to the console."""
    print("=" * 60)
    print("          Quick Calc   vs   DQ Model   Comparison")
    print("=" * 60)
    print(f"  T_nom                  = {cmp['T_nom']:.4g} Nm")
    print(f"  Quick positive peak    = {100 * cmp['quick_pos_peak'] / cmp['T_nom']:.2f} % Tn")
    print(f"  DQ    positive peak    = {100 * cmp['dq_pos_peak'] / cmp['T_nom']:.2f} % Tn")
    print(f"  Quick negative peak    = {100 * cmp['quick_neg_peak'] / cmp['T_nom']:.2f} % Tn")
    print(f"  DQ    negative peak    = {100 * cmp['dq_neg_peak'] / cmp['T_nom']:.2f} % Tn")
    print(f"  RMS difference (full)  = {cmp['diff_rms_pct']:.2f} % Tn")
    print(f"  RMS diff (t < 50ms)    = {cmp['diff_early_rms_pct']:.2f} % Tn")
    print(f"  RMS diff (t > 100ms)   = {cmp['diff_late_rms_pct']:.2f} % Tn")
    print(f"  Max point-wise diff    = {cmp['diff_max_pct']:.2f} % Tn")
    print("=" * 60)
