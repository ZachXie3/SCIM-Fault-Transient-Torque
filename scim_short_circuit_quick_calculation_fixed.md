# SCIM Three-Phase Terminal Short Circuit — Quick Calculation Workflow

## Purpose

Estimate the dominant electrical time constants and a first-pass transient electromagnetic torque envelope after a balanced three-phase bolted terminal short circuit on a squirrel-cage induction machine (SCIM).

Use this workflow for screening and order-of-magnitude checks. Use the dq/state-space workflow when point-on-wave, first-cycle shape, speed dynamics, or detailed torque waveform is important.

---

## 1. Inputs and units

Use per-phase, stator-referred equivalent-circuit quantities.

| Symbol | Python variable | Meaning | Unit / note |
|---|---|---|---|
| $R_s$ | `Rs` | stator resistance per phase | $\Omega$ |
| $R_r$ | `Rr` | rotor resistance referred to stator | $\Omega$, not divided by slip |
| $X_s$ | `Xs` | stator leakage reactance | $\Omega$ at frequency $f$ |
| $X_r$ | `Xr` | rotor leakage reactance referred to stator | $\Omega$ at frequency $f$ |
| $X_m$ | `Xm` | magnetizing reactance | $\Omega$ at frequency $f$ |
| $s$ | `slip` | pre-fault slip | per unit |
| $J$ | `J` | total inertia referred to motor shaft | kg·m$^2$, optional in quick method |
| $V_{LL}$ | `V_LL` | pre-fault line-line RMS voltage | V |
| $f$ | `f` | electrical frequency | Hz |
| $p$ | `pole_pairs` | pole pairs | dimensionless |

Base conversions:

$$
\omega_s = 2\pi f
$$

$$
V_{\phi}=\frac{V_{LL}}{\sqrt{3}}
$$

---

## 2. Initial steady state from slip

Rotor branch at slip $s$:

$$
Z_r(s)=\frac{R_r}{s}+jX_r
$$

Magnetizing branch:

$$
Z_m=jX_m
$$

Parallel air-gap branch:

$$
Z_p=Z_m \parallel Z_r=\frac{Z_m Z_r}{Z_m+Z_r}
$$

Pre-fault stator current:

$$
I_{s0}=\frac{V_{\phi}}{R_s+jX_s+Z_p}
$$

Air-gap internal voltage:

$$
E_{ag0}=V_{\phi}-I_{s0}\left(R_s+jX_s\right)
$$

Rotor branch current:

$$
I_{r0}=\frac{E_{ag0}}{R_r/s+jX_r}
$$

Synchronous mechanical angular speed:

$$
\omega_{syn,m}=\frac{\omega_s}{p}
$$

Nominal/pre-fault electromagnetic torque:

$$
T_{nom}=\frac{3\left|I_{r0}\right|^2\left(R_r/s\right)}{\omega_{syn,m}}
$$

---

## 3. Short-circuit reactances and time constants

Stator-side short-circuit reactance:

$$
X_{s,sc}=X_s+\frac{X_m X_r}{X_m+X_r}
$$

Rotor-side short-circuit reactance:

$$
X_{r,sc}=X_r+\frac{X_m X_s}{X_m+X_s}
$$

Open-circuit rotor flux decay time constant:

$$
T_{r0}=\frac{X_m+X_r}{\omega_s R_r}
$$

Short-circuit rotor/ac component decay time constant:

$$
T_{r,sc}=\frac{X_{r,sc}}{\omega_s R_r}
$$

Equivalent NEMA-style form:

$$
T_{r,sc}=\frac{X_{s,sc}}{X_s+X_m}T_{r0}
$$

Stator dc/asymmetrical component decay time constant:

$$
T_{s,dc}=\frac{X_{s,sc}}{\omega_s R_s}
$$

Approximate torque envelope decay if dominated by the rotor-flux component:

$$
T_T \approx \frac{T_{r,sc}}{2}
$$

Approximate torque envelope decay for a product of stator dc and rotor ac quantities:

$$
T_T \approx \frac{T_{s,dc}T_{r,sc}}{T_{s,dc}+T_{r,sc}}
$$

---

## 4. Initial short-circuit current and torque envelope

Approximate initial symmetrical short-circuit RMS current:

$$
I_{sc,ac,0}\approx\frac{\left|E_{ag0}\right|}{\left|R_s+jX_{s,sc}\right|}
$$

Approximate air-gap flux linkage on RMS phasor basis:

$$
\left|\Psi_{m0}\right|\approx\frac{\left|E_{ag0}\right|}{\omega_s}
$$

Initial torque envelope:

$$
T_{env,0}\approx 3p\left|\Psi_{m0}\right|I_{sc,ac,0}
$$

Positive envelope:

$$
T_{env}(t)=T_{env,0}e^{-t/T_T}
$$

A signed engineering waveform for plotting is:

$$
T_e(t)\approx T_{env,0}e^{-t/T_T}\cos\left(\omega_T t+\phi\right)
$$

where the default torque oscillation frequency is:

$$
\omega_T=(1-s)\omega_s
$$

For low slip, $\omega_T\approx\omega_s$. Choose $\phi$ so that the plotted waveform starts near the pre-fault torque:

$$
\phi=\cos^{-1}\left(\operatorname{clip}\left(\frac{T_{nom}}{T_{env,0}},-1,1\right)\right)
$$

Then normalize the plotted curve as:

$$
T_{\%}(t)=100\frac{T_e(t)}{T_{nom}}
$$

---

## 5. DC component note

A balanced three-phase short circuit can still produce a dc/asymmetrical component in each phase current. This is not zero-sequence current; it is the natural response required by current and flux continuity at the instant the terminal voltage is forced to zero.

Include dc sensitivity when first-cycle peak torque matters, motor $X/R$ is high, or $T_{s,dc}$ is comparable to the first few electrical cycles. It can often be neglected for symmetrical RMS current estimates or after several $T_{s,dc}$.

The quick script includes an optional `DC_OFFSET_FACTOR` to add a conservative first-cycle asymmetry sensitivity term. The dq script naturally includes point-on-wave and dc/asymmetry through state continuity.

---

## 6. Script output

The companion Python script:

1. declares user-editable parameters at the top;
2. computes $T_{nom}$ from slip;
3. computes $T_{r0}$, $T_{r,sc}$, and $T_{s,dc}$;
4. estimates short-circuit current and torque envelope;
5. plots $T_{fault}/T_{nominal}$ in percent from $0$ to $0.2$ s;
6. saves a CSV and PNG file.
