"""
Sweep fault inception angle to find the worst-case transient torque.

Uses a two-phase approach for efficiency:
  1. Coarse sweep over [0, 180) degrees at a user-specified step.
  2. Fine refinement around the coarse worst-case angle.

The peak absolute torque is invariant with inception angle for this motor
(classical result for balanced faults), but the positive/negative asymmetry
varies. This sweep quantifies the asymmetry range.
"""

import numpy as np

from .dq import run_dq_simulation


def angle_sweep(p, slip=None, coarse_step=2.0, refine_width=4.0, refine_step=0.5):
    """Sweep INITIAL_VOLTAGE_ANGLE_DEG and find the worst-case torque.

    Parameters
    ----------
    p : dict — normalized motor parameters.
    slip : float, optional.
    coarse_step : float — step size (deg) for initial coarse sweep (default 2).
    refine_width : float — half-width (deg) for refinement (default 4).
    refine_step : float — step size (deg) for refinement (default 0.5).

    Returns
    -------
    dict with keys:
        angles, pos_peaks, neg_peaks, abs_peaks, T_nom,
        worst_angle_deg, worst_abs_Nm, worst_abs_pct,
        worst_neg_angle, worst_neg_Nm, worst_neg_pct,
        worst_pos_angle, worst_pos_Nm, worst_pos_pct
    """
    # Reduce temporal resolution during sweep for speed
    sweep_n = max(500, p["N_POINTS"] // 10)
    sweep_t = min(0.1, p["T_END"])

    angles0 = np.arange(0.0, 180.0, coarse_step)
    T_nom = None
    pos0, neg0, abs0 = [], [], []

    orig_n = p["N_POINTS"]
    orig_t = p["T_END"]
    p["N_POINTS"] = sweep_n
    p["T_END"] = sweep_t

    # ---- Phase 1: coarse sweep ----
    for a in angles0:
        res = run_dq_simulation(p, slip=slip)
        if T_nom is None:
            T_nom = res["T_nom"]
        Tf = res["T_fault"]
        pos0.append(np.max(Tf))
        neg0.append(np.min(Tf))
        abs0.append(np.max(np.abs(Tf)))

    p["N_POINTS"] = orig_n
    p["T_END"] = orig_t

    ci = int(np.argmax(abs0))
    ca = angles0[ci]

    # ---- Phase 2: refinement around coarse worst ----
    ra = np.arange(
        max(0, ca - refine_width),
        min(180, ca + refine_width + refine_step),
        refine_step,
    )

    # Merge coarse and refined results
    pmap = dict(zip(angles0, pos0))
    nmap = dict(zip(angles0, neg0))
    amap = dict(zip(angles0, abs0))
    for a in ra:
        res = run_dq_simulation(p, slip=slip)
        Tf = res["T_fault"]
        pmap[a] = np.max(Tf)
        nmap[a] = np.min(Tf)
        amap[a] = np.max(np.abs(Tf))

    all_a = np.array(sorted(set(angles0) | set(ra)))
    pos = np.array([pmap[x] for x in all_a])
    neg = np.array([nmap[x] for x in all_a])
    ab = np.array([amap[x] for x in all_a])

    wi = int(np.argmax(ab))
    wni = int(np.argmin(neg))
    wpi = int(np.argmax(pos))

    return {
        "angles": all_a,
        "pos_peaks": pos,
        "neg_peaks": neg,
        "abs_peaks": ab,
        "T_nom": T_nom,
        "worst_angle_deg": all_a[wi],
        "worst_abs_Nm": ab[wi],
        "worst_abs_pct": 100.0 * ab[wi] / T_nom,
        "worst_neg_angle": all_a[wni],
        "worst_neg_Nm": neg[wni],
        "worst_neg_pct": 100.0 * neg[wni] / T_nom,
        "worst_pos_angle": all_a[wpi],
        "worst_pos_Nm": pos[wpi],
        "worst_pos_pct": 100.0 * pos[wpi] / T_nom,
    }
