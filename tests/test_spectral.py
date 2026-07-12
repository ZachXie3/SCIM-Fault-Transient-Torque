"""
Test D: FFT validation of torque waveform spectral content.

Computes the FFT of the simulated torque waveform and verifies that
the dominant spectral components are physically reasonable.
"""

import numpy as np
import pytest

from scim_calc.dq import run_dq_simulation
from scim_calc.spectral import compute_fft, find_peaks


def test_fft_dominant_peak_near_line_frequency(motor_params):
    """The dominant spectral peak should be near the electrical fundamental
    frequency (60 Hz for the test motor)."""
    d = run_dq_simulation(motor_params, slip=motor_params["slip"])
    T_pct = 100.0 * d["T_fault"] / d["T_nom"]
    freqs, mag = compute_fft(d["t"], T_pct)
    peaks = find_peaks(freqs, mag, min_rel_amplitude=0.01)

    assert len(peaks) >= 1, "No spectral peaks found above 1% threshold"
    f_dominant = peaks[0][0]
    # Dominant peak should be near 60 Hz
    assert 55 < f_dominant < 65, f"Dominant peak at {f_dominant:.2f} Hz (expected ~60 Hz)"


def test_fft_slip_peak_present(motor_params):
    """A secondary peak near the slip frequency should be present
    or the spectrum should show energy at low frequencies."""
    d = run_dq_simulation(motor_params, slip=motor_params["slip"])
    T_pct = 100.0 * d["T_fault"] / d["T_nom"]
    freqs, mag = compute_fft(d["t"], T_pct)
    peaks = find_peaks(freqs, mag, min_rel_amplitude=0.01)

    # At least one peak should be significantly below line frequency
    low_freq_peaks = [p for p in peaks if p[0] < 40]
    assert len(low_freq_peaks) >= 1, "No low-frequency spectral content detected"


def test_fft_amplitude_reasonable(motor_params):
    """The peak spectral amplitude should be on the order of the
    peak-to-peak torque variation."""
    d = run_dq_simulation(motor_params, slip=motor_params["slip"])
    T_pct = 100.0 * d["T_fault"] / d["T_nom"]
    freqs, mag = compute_fft(d["t"], T_pct)
    peaks = find_peaks(freqs, mag, min_rel_amplitude=0.01)

    if len(peaks) > 0:
        # The RMS of the AC component should be comparable to the
        # dominant FFT amplitude
        T_ac = T_pct - np.mean(T_pct)
        rms = np.sqrt(np.mean(T_ac ** 2))
        assert peaks[0][1] > rms * 0.5, f"Dominant peak {peaks[0][1]:.2f} too small vs RMS {rms:.2f}"


def test_fft_no_unexpected_high_frequency_content(motor_params):
    """Above 500 Hz, the spectral amplitude should be negligible
    (no high-frequency numerical noise)."""
    d = run_dq_simulation(motor_params, slip=motor_params["slip"])
    T_pct = 100.0 * d["T_fault"] / d["T_nom"]
    freqs, mag = compute_fft(d["t"], T_pct)

    high_freq = freqs > 500
    if np.any(high_freq):
        max_high = np.max(mag[high_freq])
        max_all = np.max(mag)
        ratio = max_high / max_all if max_all > 0 else 0
        assert ratio < 0.05, f"High-frequency content at {ratio*100:.2f}% of peak (>5%)"
