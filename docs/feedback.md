# Feedback to Programming Agent — d–q Transient Model Documentation Fact Check

**Target document:** `d_q_transient_model.md`  
**Purpose:** Provide actionable technical corrections for the programming/documentation agent before this document is treated as a reliable engineering reference.  
**Review focus:** d–q / stationary-frame induction machine transient model, balanced three-phase short-circuit behavior, numerical method claims, and implementation-facing documentation consistency.

---

## Executive Summary

The document is a strong tutorial draft, but several statements should be corrected before release. The highest-priority issue is the treatment of **fault inception angle**. For the stated model — ideal balanced three-phase terminal short, linear magnetic circuit, symmetric machine, stationary αβ/dq equations — changing the initial voltage angle is only a common rotation of all space vectors. The scalar electromagnetic torque should therefore be **invariant to fault inception angle**, not merely “nearly” invariant in peak absolute value.

The document also incorrectly states that current can change discontinuously at fault inception. In the presented flux-linkage model, if flux linkages are continuous and the inductance matrix does not change, currents computed from the same algebraic flux-current relation are also continuous. The discontinuity is in the applied voltage boundary condition and derivatives, not in current.

---

## Required Corrections

### 1. Correct the fault-inception-angle discussion

**Affected sections:**

- Section 9.3 — Relationship Between Inception Angle and DC Offset Magnitude
- Section 10.3 — Effect on the Transient Waveform
- Section 10.4 — Periodicity: Why 180°
- Section 10.5 — Why Peak |T| Is Nearly Angle-Invariant for Balanced Faults
- Section 10.6 — Asymmetry Mechanism
- Section 10.7 — Angle Sweep Implementation
- Section 15.2 — Quantified Differences for the Test Motor
- Section 16.6 — Fault Inception Angle θ₀

**Issue:**  
The document says the fault inception angle changes the positive and negative torque peaks while peak absolute torque is nearly invariant. For the ideal balanced, linear, symmetric stationary-frame model, this is not correct. A change in θ₀ multiplies the initial flux and current space vectors by a common factor `exp(jθ₀)`. The ODEs are rotationally symmetric, so the full trajectory is rotated by the same factor. The torque expression is invariant to that common rotation:

```math
T_e = -\frac{3}{2}p\,\operatorname{Im}\{\psi_s i_s^*\}
```

If:

```math
\psi_s' = e^{j\theta_0}\psi_s, \qquad i_s' = e^{j\theta_0}i_s
```

then:

```math
\psi_s' i_s'^* = \left(e^{j\theta_0}\psi_s\right)\left(e^{j\theta_0}i_s\right)^* = \psi_s i_s^*
```

Therefore, the entire electromagnetic torque waveform should be invariant, not just the peak absolute value.

**Recommended replacement text:**

> For the ideal balanced, linear, symmetric three-phase terminal short-circuit model, the electromagnetic torque waveform is invariant to the global fault inception angle. Changing θ₀ rotates the initial stator and rotor flux/current space vectors together in the stationary αβ plane. Because the torque expression is a rotationally invariant cross product, the scalar torque waveform is unchanged. The phase currents and axis components redistribute with θ₀, but the electromagnetic torque does not.
>
> Angle dependence may appear if the model includes unbalanced faults, phase-asymmetric winding parameters, saturation with spatial saliency or axis-dependent effects, non-sinusoidal winding distribution, zero-sequence paths, supply/network asymmetry, or phase-specific protection quantities.

**Implementation action:**

- If `sweep.py` is intended only for balanced three-phase bolted terminal faults in the present model, either remove the angle sweep or convert it into a validation routine that confirms angle invariance.
- If angle sweep is retained for future extensions, clearly label it as useful for unbalanced/asymmetric/saturated extensions, not necessary for the current ideal balanced model.

---

### 2. Correct the statement that current can jump discontinuously

**Affected section:** Section 8.4 — Why Current Can Change Discontinuously but Flux Cannot

**Issue:**  
The document says current can change discontinuously while flux remains continuous. In the stated model, currents are algebraically recovered from fluxes using a constant inductance matrix:

```math
i_s = \frac{L_r\psi_s - L_m\psi_r}{\Delta}
```

```math
i_r = \frac{-L_m\psi_s + L_s\psi_r}{\Delta}
```

If `ψs` and `ψr` are continuous and `Ls`, `Lr`, `Lm`, and `Δ` do not change at the fault, then `is` and `ir` are also continuous.

**Recommended replacement title:**

> Why Flux and Current Are Continuous but Their Derivatives Change

**Recommended replacement text:**

> Flux linkages are continuous at fault inception because an instantaneous flux jump would require infinite voltage. In this linear model, the inductance matrix is unchanged at the instant of fault. Therefore, stator and rotor currents, which are algebraic functions of the continuous flux linkages, are also continuous. What changes abruptly is the voltage boundary condition. As a result, `dψs/dt`, `dψr/dt`, and the subsequent current derivatives change immediately after the fault.

---

### 3. Re-check the expanded mechanical equation sign

**Affected section:** Section 4.5 — The Complete System in Real Form

**Issue:**  
The complex torque equation is given as:

```math
T_e = -\frac{3}{2}p\,\operatorname{Im}\{\psi_s i_s^*\}
```

Using:

```math
i_s = \frac{L_r\psi_s - L_m\psi_r}{\Delta}
```

the `Lr ψs` self-term contributes no imaginary cross-product torque because `ψs ψs*` is real. The torque should reduce to a mutual-flux term proportional to:

```math
\psi_{qs}\psi_{dr} - \psi_{ds}\psi_{qr}
```

up to the sign convention used by `torque_from_flux_current()`.

The row-5 expression in the current document appears to have the opposite sign relative to the stated complex torque expression, and the inclusion of `Lr` terms in the displayed matrix row is confusing because those terms cancel in the torque expression.

**Recommended correction:**

Avoid presenting the fifth row as a linear matrix row. The system is nonlinear when speed dynamics are enabled, and the safest documentation is:

```math
\dot{\omega}_m = \frac{1}{J}\left(T_e - T_L\right)
```

with:

```math
T_e = -\frac{3}{2}p\,\operatorname{Im}\{\psi_s i_s^*\}
```

If an expanded real-form torque is included, derive it directly from the code’s `torque_from_flux_current()` function and add a unit test confirming sign consistency.

**Implementation action:**

Add a regression test:

```python
# pseudo-test intent
Te_complex = torque_from_flux_current(psi_s, i_s, pole_pairs)
Te_expanded = expanded_real_form_torque(psi_ds, psi_qs, psi_dr, psi_qr, params)
assert np.isclose(Te_complex, Te_expanded)
```

---

### 4. Fix RK4 stability discussion

**Affected sections:**

- Section 11.5 — Stability Region
- Section 11.6 — Step Size Selection for This Application

**Issue:**  
The document uses approximately `12 s^-1` as the largest eigenvalue magnitude. This appears to be a damping rate, not the full complex eigenvalue magnitude. For modes near 60 Hz, the imaginary part is approximately:

```math
\omega_s = 2\pi \cdot 60 \approx 377\ \text{rad/s}
```

Therefore the complex eigenvalue magnitude is dominated by hundreds of rad/s, not only the damping rate. The stability conclusion is still likely correct because:

```math
377 \times 40\ \mu\text{s} \approx 0.015
```

which is still very small for RK4. But the current calculation is technically incorrect.

**Recommended replacement text:**

> The RK4 stability check should use the full complex eigenvalues of the linearized 4-state electrical system, not only the damping rates. For each eigenvalue λ, compute `z = λh` and verify that the RK4 stability function satisfies `|R(z)| < 1`, where:
>
> ```math
> R(z)=1+z+\frac{z^2}{2}+\frac{z^3}{6}+\frac{z^4}{24}
> ```
>
> With `h ≈ 40 μs`, the expected 60 Hz electrical components give `|z| ≈ 0.015`, which is well inside the RK4 stability region.

**Implementation action:**

Add a diagnostic function that computes and prints:

- eigenvalues of the constant-speed 4-state electrical system,
- `max(abs(lambda))`,
- `max(abs(R(lambda*h)))`,
- recommended stability margin.

---

### 5. Remove or soften the “below machine precision” RK4 error claim

**Affected section:** Section 11.6 — Step Size Selection for This Application

**Issue:**  
The document states that the local truncation error is far below double-precision machine precision and that actual error is dominated by round-off. This is too strong unless proven by convergence testing.

**Recommended replacement text:**

> With `h ≈ 40 μs`, RK4 is expected to provide high accuracy for the 60 Hz electrical transient and the dominant exponential decays. Accuracy should be verified using a step-halving convergence test, comparing peak torque, time of peak torque, and RMS waveform difference between `N_POINTS`, `2*N_POINTS`, and `4*N_POINTS`.

**Implementation action:**

Add an automated convergence check:

```python
# pseudo-test intent
results_5k = run_dq_simulation(..., n_points=5000)
results_10k = run_dq_simulation(..., n_points=10000)
results_20k = run_dq_simulation(..., n_points=20000)

# interpolate to common time grid and compare
# report peak torque convergence and RMS waveform convergence
```

---

### 6. Revise skin-effect limitation

**Affected section:** Section 17.5 — No Skin Effect

**Issue:**  
The document says skin effect is negligible because rotor frequency is close to slip frequency. That is valid for normal steady-state operation near rated slip, but it is not generally valid during a short-circuit transient. A stator-frame DC-offset component is stationary relative to the stator and therefore is seen by the rotor at approximately rotor electrical speed, near line frequency for a near-synchronous 2-pole motor.

**Recommended replacement text:**

> Skin effect may be negligible for low-slip steady-state rotor currents, but short-circuit transients can include rotor-current components at higher effective frequencies, including components associated with stator-frame DC offset. For deep-bar or double-cage rotors, frequency-dependent rotor resistance and leakage may affect early transient torque. The single-cage constant-parameter model should therefore be treated as an approximation.

---

### 7. Revise stationary-frame justification

**Affected section:** Section 1.6 — Why the Stationary Frame for This Application

**Issue:**  
The document implies the short-circuit condition `vs = 0` is uniquely simple in the stationary frame and more complicated in the synchronous frame. A zero vector remains zero under any orthogonal frame transformation. The stationary frame is still convenient, but the stated reason should be refined.

**Recommended replacement text:**

> The stationary frame is convenient because no stator frame-rotation voltage term appears, and the stator terminal voltage boundary condition can be written directly in stator coordinates. The zero-voltage condition itself remains zero in any orthogonal reference frame; the main difference is the location of the speed-voltage terms and the transformation of nonzero pre-fault steady-state quantities.

---

### 8. Clarify “bolted three-phase short” wording

**Affected sections:**

- Section 8.1 — What a Bolted Three-Phase Short Circuit Means
- Section 17.2 — Balanced Fault Only

**Issue:**  
The document describes a bolted three-phase short as all three phases shorted to each other and to ground. For the positive-sequence balanced dq model, the essential condition is zero terminal voltage space vector / zero line-to-line terminal voltage. Grounding only matters if zero-sequence/common-mode paths are included.

**Recommended replacement text:**

> A balanced bolted three-phase terminal short forces the stator terminal voltage space vector to zero. In a three-wire positive-sequence dq model, this is equivalent to zero line-to-line voltage at the motor terminals. Whether the short point is grounded is not material unless zero-sequence paths are explicitly modeled.

---

### 9. Revise leakage coefficient interpretation

**Affected section:** Section 2.3 — Leakage vs Magnetizing Inductance

**Issue:**  
The document says `σ ≈ 0.051` indicates a loose-coupled machine. That wording is misleading. A leakage coefficient near 0.05 is typical of an induction motor with strong magnetizing coupling and finite leakage.

**Recommended replacement text:**

> The leakage coefficient is approximately `σ ≈ 0.051`, indicating typical induction-machine leakage with strong magnetizing coupling and finite leakage flux.

---

### 10. Verify torque spectral labels using FFT instead of verbal assumptions

**Affected section:** Section 15.3 — Component-by-Component Spectral Analysis

**Issue:**  
The document uses labels such as rotor-frequency, line-frequency, and slip-frequency components, but the frequency definitions are not always consistent. Earlier sections describe a rotor mode near slip frequency, while later sections refer to rotor-frequency components near `(1-s)ωs`. This can confuse users.

**Recommended correction:**

Replace the purely verbal spectral description with an FFT-based validation result from the actual simulated torque waveform.

**Implementation action:**

Add a post-processing utility that reports:

- dominant frequency peaks in `T_fault`,
- corresponding amplitudes,
- whether the detected frequencies match the expected modal/electrical components,
- warning if labels are ambiguous.

---

## Recommended New Validation Tests

### A. Angle-invariance test

Run the same balanced three-phase terminal short at multiple values of `INITIAL_VOLTAGE_ANGLE_DEG`, for example:

```text
0°, 15°, 30°, 60°, 90°, 135°, 180°
```

For the current ideal balanced model, verify:

```python
max_abs_difference(T_fault_angle_a, T_fault_angle_b) < tolerance
```

Expected result: torque traces should be identical within numerical tolerance.

---

### B. Current-continuity test

Verify that current immediately before and immediately after the fault is continuous when computed from the same flux states and same inductance matrix.

Expected result:

```python
is_0_minus == is_0_plus
ir_0_minus == ir_0_plus
```

within floating-point tolerance. The derivatives should change due to the voltage boundary condition.

---

### C. Torque sign-consistency test

Compare:

1. complex function implementation,
2. expanded real-form expression,
3. expected motoring sign at pre-fault condition.

Expected result: all three agree.

---

### D. RK4 convergence test

Run the simulation with at least three time-step resolutions:

```text
N_POINTS = 5,000
N_POINTS = 10,000
N_POINTS = 20,000
```

Compare:

- peak positive torque,
- peak negative torque,
- peak absolute torque,
- time of peak absolute torque,
- RMS waveform difference after interpolation to a common grid.

Expected result: differences should decrease approximately with 4th-order behavior until round-off or interpolation error dominates.

---

### E. Eigenvalue/stability diagnostic

For constant-speed mode, compute the 4-state electrical system matrix and evaluate:

```math
R(z)=1+z+\frac{z^2}{2}+\frac{z^3}{6}+\frac{z^4}{24}
```

where:

```math
z = \lambda h
```

Expected result:

```python
max(abs(R(lambda_i * h))) < 1
```

with a comfortable stability margin.

---

## Suggested Documentation Restructure

To avoid future confusion, restructure the controversial sections as follows:

```text
9. DC Offset
   9.1 Origin of DC Offset in Phase Currents
   9.2 Flux Continuity and Current Continuity
   9.3 DC Offset Distribution Among Phases
   9.4 Why Balanced Three-Phase Torque Is Angle-Invariant

10. Fault Inception Angle
   10.1 Physical Meaning of θ₀
   10.2 Effect on Phase Quantities
   10.3 No Effect on Electromagnetic Torque in the Ideal Balanced Model
   10.4 When Angle Dependence Can Appear

15. Validation and Analysis
   15.1 Quick vs d–q Method Differences
   15.2 Angle-Invariance Validation
   15.3 RK4 Step-Convergence Validation
   15.4 FFT-Based Spectrum Check
```

---

## Priority List for Programming Agent

### Must fix before release

1. Correct the fault-angle invariance issue.
2. Correct the current-discontinuity statement.
3. Verify and correct the expanded mechanical torque sign.
4. Fix the RK4 stability calculation.
5. Remove unsupported “machine precision” accuracy claims.

### Should fix soon

6. Revise skin-effect limitation.
7. Clarify stationary-frame justification.
8. Clarify bolted-fault wording.
9. Improve leakage coefficient interpretation.
10. Replace verbal spectral claims with FFT-backed validation.

---

## Final Note to Agent

The code may still be structurally correct even if the explanation is not. The most important task is to ensure the documentation matches the mathematical invariances of the implemented model. In particular, if the current simulation produces torque traces that vary with `INITIAL_VOLTAGE_ANGLE_DEG` for an ideal balanced short, that is a red flag and should be investigated as a possible implementation inconsistency, sign convention problem, or non-rotationally-invariant initialization.
