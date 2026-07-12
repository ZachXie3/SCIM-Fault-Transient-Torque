"""
Test E: Eigenvalue/stability diagnostic.

For the constant-speed (speed dynamics disabled) linearized 4-state
electrical system, compute the system matrix eigenvalues and verify
that the RK4 stability condition is satisfied at the chosen step size.
"""

import numpy as np

from scim_calc.dq import run_dq_simulation


def _system_matrix(p, omega_m):
    """Build the 4x4 electrical system matrix for constant-speed operation.
    State order: [psi_ds, psi_qs, psi_dr, psi_qr]."""
    Ls, Lr, Lm, Delta = p["Ls"], p["Lr"], p["Lm"], p["Delta"]
    Rs, Rr = p["Rs"], p["Rr"]
    p_pole = p["pole_pairs"]
    ore = p_pole * omega_m  # rotor electrical speed

    return np.array([
        [-Rs*Lr/Delta,       0,          Rs*Lm/Delta,       0        ],
        [      0,        -Rs*Lr/Delta,       0,          Rs*Lm/Delta ],
        [ Rr*Lm/Delta,       0,         -Rr*Ls/Delta,    -ore      ],
        [      0,         Rr*Lm/Delta,     ore,        -Rr*Ls/Delta ],
    ])


def _rk4_stability_function(z):
    """RK4 stability function R(z) = 1 + z + z^2/2 + z^3/6 + z^4/24."""
    return 1.0 + z + z**2/2.0 + z**3/6.0 + z**4/24.0


def test_eigenvalues_computed(motor_params):
    """The system matrix eigenvalues should be computable and sensible."""
    p = motor_params
    d = run_dq_simulation(p, slip=p["slip"])
    omega_m = d["omega_m_trace"][0]

    A = _system_matrix(p, omega_m)
    evals = np.linalg.eigvals(A)

    # Should have two complex-conjugate pairs
    assert len(evals) == 4, f"Expected 4 eigenvalues, got {len(evals)}"
    # All eigenvalues should have negative real parts (stable system)
    assert np.all(evals.real < 0), "All eigenvalues should have negative real parts"
    # Should have at least one complex pair
    assert np.any(np.abs(evals.imag) > 1.0), "Expected at least one oscillatory mode"


def test_rk4_stability_condition(motor_params):
    """|R(lambda * h)| should be < 1 for all eigenvalues at the default step size."""
    p = motor_params
    d = run_dq_simulation(p, slip=p["slip"])
    omega_m = d["omega_m_trace"][0]

    A = _system_matrix(p, omega_m)
    evals = np.linalg.eigvals(A)

    h = p["T_END"] / (p["N_POINTS"] - 1)

    max_R = 0.0
    for lam in evals:
        z = lam * h
        Rz = _rk4_stability_function(z)
        max_R = max(max_R, abs(Rz))

    assert max_R < 1.0, f"Max |R(z)| = {max_R:.6f} (should be < 1 for stability)"
    # Stability margin: how close to the stability boundary
    margin = 1.0 - max_R
    assert margin > 0.0, f"No stability margin (|R| = {max_R:.6e})"


def test_eigenvalue_magnitude_reasonable(motor_params):
    """The largest eigenvalue magnitude should be on the order of
    the electrical frequency (~377 rad/s for 60 Hz)."""
    p = motor_params
    d = run_dq_simulation(p, slip=p["slip"])
    omega_m = d["omega_m_trace"][0]

    A = _system_matrix(p, omega_m)
    evals = np.linalg.eigvals(A)

    max_mag = np.max(np.abs(evals))
    # For a 60 Hz machine, the eigenvalues should have magnitudes
    # in the range of hundreds of rad/s (not ~12 as previously claimed)
    assert max_mag > 100, f"Max eigenvalue magnitude {max_mag:.1f} should be >100 rad/s"
    assert max_mag < 1000, f"Max eigenvalue magnitude {max_mag:.1f} should be <1000 rad/s"


def test_eigenvalue_stability_diagnostic_print(motor_params, capsys):
    """Print the stability diagnostic for documentation purposes."""
    p = motor_params
    d = run_dq_simulation(p, slip=p["slip"])
    omega_m = d["omega_m_trace"][0]

    A = _system_matrix(p, omega_m)
    evals = np.linalg.eigvals(A)

    h = p["T_END"] / (p["N_POINTS"] - 1)

    print(f"\nEigenvalue analysis (constant-speed, omega_m={omega_m:.2f} rad/s):")
    print(f"  Step size h = {h:.6e} s")
    for i, lam in enumerate(evals):
        z = lam * h
        Rz = _rk4_stability_function(z)
        print(f"  lambda[{i}] = {lam.real:.4f} +/- {abs(lam.imag):.4f}j, "
              f"|lambda| = {abs(lam):.4f}, |z| = {abs(z):.6e}, |R(z)| = {abs(Rz):.6e}")

    max_R = max(abs(_rk4_stability_function(lam * h)) for lam in evals)
    print(f"  Max |R(z)| = {max_R:.6e} (should be < 1)")
    print(f"  Max |lambda| = {np.max(np.abs(evals)):.2f} rad/s")
    print(f"  Max |z| = {np.max(np.abs(evals * h)):.6e}")

    captured = capsys.readouterr()
    assert "Eigenvalue analysis" in captured.out
    assert "Max |R(z)|" in captured.out
