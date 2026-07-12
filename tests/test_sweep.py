"""
Test the angle-sweep module.

Verifies that the two-phase coarse/refine sweep produces sensible results:
nominal torque matches, asymmetric peaks are detected, and the peak absolute
torque is nearly invariant with inception angle.
"""

import numpy as np
import pytest

from scim_calc.sweep import angle_sweep


def test_sweep_tnom_matches(motor_params, expected_sweep):
    """Nominal torque from the sweep should match the saved baseline."""
    s = angle_sweep(
        motor_params, slip=motor_params["slip"],
        coarse_step=10.0, refine_width=2.0, refine_step=1.0,
    )
    assert s["T_nom"] == pytest.approx(expected_sweep["T_nom"], rel=1e-2)


def test_sweep_peaks_reasonable(motor_params):
    """Negative peak should be more severe than positive, and both should
    be within expected bounds for the test motor."""
    s = angle_sweep(
        motor_params, slip=motor_params["slip"],
        coarse_step=10.0, refine_width=2.0, refine_step=1.0,
    )
    assert abs(s["worst_neg_pct"]) > s["worst_pos_pct"]
    assert s["worst_neg_pct"] < -400.0  # >400% negative
    assert s["worst_pos_pct"] < 350.0   # <350% positive


def test_sweep_angle_range(motor_params):
    """The returned angle array should span [0, 180) deg with reasonable size."""
    s = angle_sweep(
        motor_params, slip=motor_params["slip"],
        coarse_step=10.0, refine_width=2.0, refine_step=1.0,
    )
    assert np.all(s["angles"] >= 0.0)
    assert np.all(s["angles"] <= 180.0)
    assert len(s["angles"]) > 10


def test_sweep_peak_abs_uniform(motor_params):
    """The peak absolute torque should be nearly constant across all inception
    angles (variation < 0.1 pu of T_nom)."""
    s = angle_sweep(
        motor_params, slip=motor_params["slip"],
        coarse_step=10.0, refine_width=2.0, refine_step=1.0,
    )
    variation = (np.max(s["abs_peaks"]) - np.min(s["abs_peaks"])) / s["T_nom"]
    assert variation < 0.1
