"""
Spectral analysis utilities for SCIM short-circuit torque waveforms.

Provides FFT-based frequency-domain analysis of the simulated torque
waveform to identify dominant spectral components and validate the
modal interpretation of the transient.
"""

import numpy as np


def compute_fft(t, signal):
    """Compute single-sided amplitude spectrum of a time-domain signal.

    Parameters
    ----------
    t : ndarray — time vector (s).
    signal : ndarray — signal values (same length as t).

    Returns
    -------
    freqs : ndarray — positive frequencies (Hz).
    mag : ndarray — single-sided amplitude magnitudes (same units as signal).
    """
    N = len(t)
    dt = t[1] - t[0]
    fs = 1.0 / dt

    fft_vals = np.fft.fft(signal - np.mean(signal))
    fft_freqs = np.fft.fftfreq(N, d=dt)

    pos = fft_freqs >= 0
    freqs = fft_freqs[pos]
    mag = np.abs(fft_vals[pos]) / N * 2.0  # single-sided amplitude
    return freqs, mag


def find_peaks(freqs, mag, min_rel_amplitude=0.01):
    """Find dominant spectral peaks above a relative amplitude threshold.

    Parameters
    ----------
    freqs : ndarray — frequencies (Hz).
    mag : ndarray — amplitudes.
    min_rel_amplitude : float — minimum peak amplitude as fraction of max.

    Returns
    -------
    list of (freq_Hz, amplitude, relative_percent) sorted by amplitude descending.
    """
    threshold = min_rel_amplitude * np.max(mag)
    peaks = []
    for i in range(1, len(mag) - 1):
        if mag[i] > mag[i-1] and mag[i] > mag[i+1] and mag[i] > threshold:
            peaks.append((freqs[i], mag[i], 100.0 * mag[i] / np.max(mag)))
    peaks.sort(key=lambda x: -x[1])
    return peaks
