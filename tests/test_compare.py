"""
Test the comparison module that quantifies differences between quick and dq results.

Verifies that the metrics are physically sensible and the print function
produces the expected output format.
"""

import pytest

from scim_calc.quick import run_quick_calc
from scim_calc.dq import run_dq_simulation
from scim_calc.compare import compute_comparison, print_comparison


def test_comparison_metrics(motor_params):
    """Comparison metrics should reflect known method differences."""
    q = run_quick_calc(motor_params, slip=motor_params["slip"])
    d = run_dq_simulation(motor_params, slip=motor_params["slip"])
    cmp = compute_comparison(q, d)

    assert cmp["T_nom"] == pytest.approx(3989, rel=1e-2)
    # DQ negative peak should be more severe than quick
    assert abs(cmp["dq_neg_peak"]) > abs(cmp["quick_neg_peak"])
    # RMS difference should be positive (methods are not identical)
    assert cmp["diff_rms_pct"] > 0
    # Late-time disagreement should be smaller than early-time
    if cmp["diff_late_rms_pct"] is not None and cmp["diff_early_rms_pct"] is not None:
        assert cmp["diff_late_rms_pct"] < cmp["diff_early_rms_pct"]


def test_comparison_print(capsys, motor_params):
    """The print function should produce correctly formatted output."""
    q = run_quick_calc(motor_params, slip=motor_params["slip"])
    d = run_dq_simulation(motor_params, slip=motor_params["slip"])
    cmp = compute_comparison(q, d)
    print_comparison(cmp)
    captured = capsys.readouterr()
    assert "Quick Calc" in captured.out
    assert "DQ" in captured.out
    assert "416" in captured.out or "376" in captured.out
