# SCIM Three-Phase Terminal Short Circuit — Quick Calculation Workflow

## Purpose

Estimate the dominant electrical time constants and a first-pass transient electromagnetic torque envelope after a balanced three-phase bolted terminal short circuit on a squirrel-cage induction machine (SCIM).

Use this workflow for screening and order-of-magnitude checks. Use the dq/state-space workflow when point-on-wave, first-cycle shape, speed dynamics, or detailed torque waveform is important.

---

## 1. Implementation

The quick calculation is implemented in `scim_calc/quick.py`, function `run_quick_calc()`.
It is invoked automatically by the top-level entry point `run.py` alongside the dq model.

### Inputs

Motor parameters are read from `input.jsonc` at the project root.
See that file for the full list of per-phase, stator-referred equivalent-circuit quantities.

| Symbol | `input.jsonc` key | Meaning | Unit |
|---|---|---|---|
| $R_s$ | `Rs` | stator resistance per phase | $\Omega$ |
| $R_r$ | `Rr` | rotor resistance referred to stator | $\Omega$ |
| $X_s$ | `Xs` | stator leakage reactance | $\Omega$ at frequency $f$ |
| $X_r$ | `Xr` | rotor leakage reactance referred to stator | $\Omega$ at frequency $f$ |
| $X_m$ | `Xm` | magnetizing reactance | $\Omega$ at frequency $f$ |
| $s$ | `slip` | pre-fault slip (see below) | per unit |
| $V_{LL}$ | `V_LL` | pre-fault line-line RMS voltage | V |
| | `CONNECTION` | stator winding connection | `"wye"` or `"delta"` |
| $f$ | `f` | electrical frequency | Hz |
| $p$ | `poles` | total number of poles | 2, 4, 6… |
| $P_{HP}$ | `HP` | rated mechanical output (alt. `kW`) | horsepower |
| | `DC_OFFSET_FACTOR` | first-cycle asymmetry sensitivity | pu (0 = none) |
| | `FAULT_ANGLE_DEG` | fault inception angle for dc term | degrees |
| | `TORQUE_FREQUENCY_MODE` | torque oscillation frequency | `"line"` or `"rotor"` |

> **Note:** The code accepts **either** reactances ($X_s, X_r, X_m$) **or** inductances ($L_{\ell s}, L_{\ell r}, L_m$) in `input.jsonc`. When inductances are given, they are converted internally via $X = \omega_s L$.

### Slip derivation from nameplate power

When either `HP` or `kW` is provided, `circuit.py:derive_slip()` derives the full-load slip by solving

$$
T_{nom}(s) = \frac{P_{mech}}{(1-s)\,\omega_{syn,m}}
$$

for $s$ using bisection. If neither `HP` nor `kW` is given, the user-supplied `slip` is used as-is.

Base conversions:

$$
\omega_s = 2\pi f ,\qquad
V_{ph}=
\begin{cases}
\dfrac{V_{LL}}{\sqrt{3}}, & \text{wye connection}\\[8pt]
V_{LL}, & \text{delta connection}
\end{cases}
$$

---

## 2. Initial steady state from slip

Rotor branch at slip $s$:

$$
Z_r(s)=\frac{R_r}{s}+jX_r
$$

Magnetizing branch and parallel air-gap impedance:

$$
Z_m=jX_m,\qquad
Z_p=Z_m \parallel Z_r=\frac{Z_m Z_r}{Z_m+Z_r}
$$

Pre-fault stator current, air-gap voltage, and rotor current:

$$
I_{s0}=\frac{V_{\phi}}{R_s+jX_s+Z_p},\qquad
E_{ag0}=V_{\phi}-I_{s0}\left(R_s+jX_s\right),\qquad
I_{r0}=\frac{E_{ag0}}{R_r/s+jX_r}
$$

Nominal torque:

$$
T_{nom}=\frac{3\left|I_{r0}\right|^2\left(R_r/s\right)}{\omega_{syn,m}},\qquad
\omega_{syn,m}=\frac{\omega_s}{p}
$$

---

## 3. Short-circuit reactances and time constants

Stator- and rotor-side short-circuit reactances:

$$
X_{s,sc}=X_s+\frac{X_m X_r}{X_m+X_r},\qquad
X_{r,sc}=X_r+\frac{X_m X_s}{X_m+X_s}
$$

Time constants:

$$
\begin{aligned}
T_{r0} &= \frac{X_m+X_r}{\omega_s R_r} &\text{(rotor open-circuit)}\\[4pt]
T_{r,sc} &= \frac{X_{r,sc}}{\omega_s R_r} &\text{(rotor short-circuit)}\\[4pt]
T_{s,dc} &= \frac{X_{s,sc}}{\omega_s R_s} &\text{(stator dc)}\\[4pt]
T_T &\approx \frac{T_{r,sc}}{2} &\text{(rotor-dominated envelope)}\\[4pt]
T_T &\approx \frac{T_{s,dc}T_{r,sc}}{T_{s,dc}+T_{r,sc}} &\text{(mixed envelope)}
\end{aligned}
$$

---

## 4. Initial short-circuit current and torque envelope

Initial symmetrical RMS current and air-gap flux (RMS phasor basis):

$$
I_{sc,ac,0}\approx\frac{\left|E_{ag0}\right|}{\left|R_s+jX_{s,sc}\right|},\qquad
\left|\Psi_{m0}\right|\approx\frac{\left|E_{ag0}\right|}{\omega_s}
$$

Initial torque envelope:

$$
T_{env,0}\approx 3p\left|\Psi_{m0}\right|I_{sc,ac,0},\qquad
T_{env}(t)=T_{env,0}e^{-t/T_T}
$$

Signed engineering waveform:

$$
T_e(t)\approx T_{env,0}e^{-t/T_T}\cos\left(\omega_T t+\phi\right),\qquad
\phi=\cos^{-1}\!\left(\operatorname{clip}\!\left(\frac{T_{nom}}{T_{env,0}},-1,1\right)\right)
$$

where $\omega_T = \omega_s$ (line mode) or $\omega_T = (1-s)\omega_s$ (rotor mode).

---

## 5. DC component note

The quick script includes an optional `DC_OFFSET_FACTOR` to add a conservative first-cycle asymmetry sensitivity term. The dq model (`scim_calc/dq.py`) naturally captures the dc component through flux continuity at fault inception without requiring a separate equation.

---

## 6. Script output

Running `python run.py` at the project root:

1. reads motor parameters from `input.jsonc`;
2. runs both the quick calculation and the dq model;
3. prints a comparison of peak torques and RMS differences;
4. saves CSV files and an overlay PNG to the `data/` directory.

The quick calculation alone can also be called programmatically:

```python
from scim_calc.config import load_jsonc, normalize_params
from scim_calc.quick import run_quick_calc

cfg = load_jsonc("input.jsonc")
p = normalize_params(cfg)
q = run_quick_calc(p)
```
