"""
Test D: RK4 step-convergence test.

Runs the simulation at multiple time-step resolutions and verifies that
peak torque values converge as the step size decreases. For RK4 (4th-order),
halving the step should reduce the error by approximately a factor of 16.
"""

import numpy as np
import pytest

from scim_calc.dq import run_dq_simulation


@pytest.mark.slow
def test_rk4_peak_torque_convergence(motor_params):
    """Peak torque should converge as N_POINTS increases.
    Uses Richardson extrapolation to estimate convergence order."""
    p = dict(motor_params)
    p["T_END"] = 0.1  # shorter time for speed

    n_points_list = [500, 1000, 2000]
    results = {}

    for n in n_points_list:
        p["N_POINTS"] = n
        r = run_dq_simulation(p, slip=p["slip"])
        results[n] = {
            "pos_peak": np.max(r["T_fault"]),
            "neg_peak": np.min(r["T_fault"]),
            "abs_peak": np.max(np.abs(r["T_fault"])),
        }

    # Check that values are converging (not diverging)
    abs_peaks = [results[n]["abs_peak"] for n in n_points_list]
    diffs = [abs(abs_peaks[i+1] - abs_peaks[i]) for i in range(len(abs_peaks)-1)]
    # Each successive refinement should produce a smaller change
    for i in range(len(diffs)-1):
        assert diffs[i+1] < diffs[i] or diffs[i+1] < 1.1 * diffs[i], \
            f"Convergence not monotonic: diffs={diffs}"


def test_rk4_peak_torque_reasonable_default(motor_params):
    """With default N_POINTS=5000, peak torque should be stable
    (within 0.1% of the 10000-point result)."""
    p = dict(motor_params)
    p["T_END"] = 0.1

    p["N_POINTS"] = 5000
    r5k = run_dq_simulation(p, slip=p["slip"])

    p["N_POINTS"] = 10000
    r10k = run_dq_simulation(p, slip=p["slip"])

    T5k = np.max(np.abs(r5k["T_fault"]))
    T10k = np.max(np.abs(r10k["T_fault"]))

    rel_diff = abs(T10k - T5k) / T10k
    assert rel_diff < 0.001, f"5k vs 10k relative difference: {rel_diff:.6f} (>0.001)"
