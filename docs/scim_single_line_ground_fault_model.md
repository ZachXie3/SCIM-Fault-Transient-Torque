# SCIM Single-Line-to-Ground Fault Torque - Implementation Instructions

## 0. Purpose

Implement a **SCIM-focused motor terminal single-line-to-ground (SLG) fault torque simulation** for the source-connected interval before breaker clearing.

This is not a ground-grid, relay coordination, or power-system short-circuit-current study. The purpose is to estimate transient electromagnetic torque during a motor-side or motor-lead ground fault while the upstream grounded source remains connected.

Typical real-life failures represented by this model:

- motor lead insulation failure to grounded frame,
- terminal lug or lead touching the grounded enclosure,
- contamination or moisture creating a phase-to-ground path,
- foreign object or animal intrusion creating a grounded fault path.

SLG is lower priority than L-L for torque, but this document is prepared so that the same unbalanced-fault framework can later support SLG after L-L is working.

---

## 1. User Inputs

For Phase 1, user-facing fault inputs shall be limited to:

```json
{
  "FAULT_TYPE": "SLG",
  "BREAKER_CLEARING_TIME_S": 0.1
}
```

All other inputs are inherited from the existing SCIM model:

```text
Rs, Rr, Xs, Xr, Xm, V_LL, CONNECTION, f, poles, slip,
ROTOR_GD2 or ROTOR_WK2, USE_SPEED_DYNAMICS, T_END, N_POINTS
```

Do not expose these as Phase-1 user inputs:

```text
FAULTED_PHASE
FAULT_INCEPTION_ANGLE_DEG
GROUNDING_MODEL
GROUND_IMPEDANCE
CABLE_IMPEDANCE
SOURCE_NETWORK_MODEL
```

---

## 2. Main Phase-1 Assumptions

### 2.1 Motor-side study, not system study

Assume the fault is at the motor terminal or motor lead close enough that cable impedance is neglected.

```text
FAULT_LOCATION = motor_terminal_or_motor_lead
CABLE_IMPEDANCE = 0
```

### 2.2 Canonical phase selection

Use a canonical A-G fault internally:

```text
canonical SLG fault = A-G fault
```

Do not ask the user to select A-G, B-G, or C-G. For a symmetric motor and balanced source, the three choices are equivalent after a 120-degree phase shift. Since the solver sweeps inception angle over 0-360 degrees, explicit faulted-phase selection does not add useful worst-case torque information in Phase 1.

Metadata shall include:

```text
SLG_FAULT_PHASE_ASSUMPTION = canonical A-G fault; phase not user-selected
```

### 2.3 Grounding assumption

For Phase 1, assume grounding exists upstream only.

```text
SLG_GROUNDING_ASSUMPTION = upstream grounded only; no motor-side ground return after breaker clearing
```

Meaning:

```text
0 <= t < t_clear
    upstream source remains connected
    upstream grounding provides ground return path

t >= t_clear
    breaker opens
    upstream ground return path is removed
    SLG simulation terminates
```

Do not ask the user to define grounding.

### 2.4 Source-connected-only SLG simulation

Unlike L-L, Phase-1 SLG shall be simulated only while the source is connected.

At breaker clearing:

```text
terminate SLG simulation
```

Reason:

```text
Under the upstream-grounded-only assumption, the conductive ground return path disappears when the breaker opens. A meaningful motor-only SLG continuation would require a motor-side ground/capacitance/leakage model, which is outside Phase 1 scope.
```

---

## 3. Reused SCIM State Model

Reuse the existing stationary alpha-beta / dq SCIM state equations for torque-producing dynamics.

Complex stationary-frame variables:

```text
psi_s = psi_s_alpha + j*psi_s_beta
psi_r = psi_r_alpha + j*psi_r_beta
i_s   = i_s_alpha   + j*i_s_beta
i_r   = i_r_alpha   + j*i_r_beta
```

Derived quantities:

```text
omega_s = 2*pi*f
pole_pairs = poles / 2
Lls = Xs / omega_s
Llr = Xr / omega_s
Lm  = Xm / omega_s
Ls  = Lls + Lm
Lr  = Llr + Lm
D_L = Ls*Lr - Lm^2
```

Flux-current inversion:

```text
i_s = (Lr*psi_s - Lm*psi_r) / D_L
i_r = (-Lm*psi_s + Ls*psi_r) / D_L
```

State equations:

```text
dpsi_s/dt = v_s - Rs*i_s
dpsi_r/dt = j*omega_re*psi_r - Rr*i_r
omega_re = pole_pairs*omega_m
```

Mechanical equation, if enabled:

```text
J*domega_m/dt = T_e - T_load
```

Default:

```text
USE_SPEED_DYNAMICS = false
```

---

## 4. Pre-Fault Initial Condition

Use the existing balanced pre-fault solution.

Steps:

```text
1. Convert V_LL to phase voltage using CONNECTION.
2. Solve pre-fault steady-state induction-machine phasor equations.
3. Construct initial psi_s0, psi_r0, and omega_m0.
4. Rotate the initial state according to the internal swept fault inception angle.
```

At fault inception:

```text
psi_s(0+) = psi_s(0-)
psi_r(0+) = psi_r(0-)
omega_m(0+) = omega_m(0-)
```

Flux states are continuous. Terminal boundary condition changes at t = 0.

---

## 5. Terminal Voltage Boundary Conditions

The SLG fault is represented through terminal phase-voltage constraints and then converted to alpha-beta voltage for the existing machine model.

### 5.1 Source-connected A-G fault

For canonical A-G fault while the upstream source is connected:

```text
0 <= t < t_clear
```

Use the following reduced terminal voltage boundary condition:

```text
v_a(t) = 0
v_b(t) = e_b(t)
v_c(t) = e_c(t)
```

where the balanced source phase references are:

```text
e_a(t) = Vpk*cos(omega_s*t + theta0)
e_b(t) = Vpk*cos(omega_s*t + theta0 - 2*pi/3)
e_c(t) = Vpk*cos(omega_s*t + theta0 + 2*pi/3)
```

This represents a motor terminal / motor lead A-G fault with upstream grounding available while the breaker is closed.

Then compute:

```text
v_alpha, v_beta = Clarke(v_a, v_b, v_c)
v_s = v_alpha + j*v_beta
```

### 5.2 Breaker clearing

At:

```text
t = t_clear
```

terminate the SLG simulation unless a future motor-side ground-return model is explicitly added.

Do not extrapolate post-clearing SLG torque by assuming the ground path remains available.

---

## 6. Torque Calculation

Compute electromagnetic torque from instantaneous alpha-beta quantities:

```text
T_e = (3/2)*pole_pairs*(psi_s_alpha*i_s_beta - psi_s_beta*i_s_alpha)
```

Do not add a separate zero-sequence torque term in Phase 1.

Reason:

```text
The existing SCIM air-gap torque model is based on alpha-beta torque-producing components. Zero-sequence ground current is primarily a protection/current quantity and is outside this Phase-1 torque model.
```

Normalize using the same convention as the existing dq model:

```text
T_percent = 100*T_e/T_nom
```

---

## 7. Angle Sweep Workflow

For each `BREAKER_CLEARING_TIME_S`, perform an internal angle sweep.

Recommended sequence:

```text
1. Coarse sweep: 0 to 360 degrees using coarse_step_deg.
2. Identify candidate angles for:
       maximum positive torque
       maximum negative torque
       maximum absolute torque
3. Refine around each candidate.
4. Re-run final high-resolution simulations for selected worst angles.
5. Save summary and worst-absolute-torque waveform.
```

Default sweep settings:

```json
{
  "coarse_step_deg": 10,
  "refine_width_deg": 10,
  "refine_step_deg": 1
}
```

---

## 8. Numerical Integration

Use the same RK4 structure as the existing dq implementation.

Important implementation requirement:

```text
Boundary condition evaluation must occur inside each RK4 substep.
```

For each `rhs(state, t, theta0)` call:

```text
1. Verify t <= t_clear for SLG simulation.
2. Compute terminal phase voltages from A-G boundary condition.
3. Convert terminal phase voltages to alpha-beta voltage.
4. Recover stator and rotor currents from flux states.
5. Compute derivatives.
```

If `BREAKER_CLEARING_TIME_S = 0`, the source-connected SLG interval has zero duration. The solver may return a no-transient result or a clear message that SLG torque is not evaluated because the breaker opens immediately and no motor-side ground path is modeled.

---

## 9. Outputs

For each final worst-angle run, write:

```text
time_s
Va_V
Vb_V
Vc_V
Ialpha_motor_A
Ibeta_motor_A
Ia_motor_A
Ib_motor_A
Ic_motor_A
I1_A
I2_A
I0_A_diagnostic
torque_Nm
torque_percent
omega_m_rad_per_s
source_connected_bool
```

Write summary output:

```text
fault_type
fault_phase_assumption
grounding_assumption
breaker_clearing_time_s
worst_positive_torque_Nm
worst_positive_torque_percent
worst_positive_angle_deg
worst_negative_torque_Nm
worst_negative_torque_percent
worst_negative_angle_deg
worst_abs_torque_Nm
worst_abs_torque_percent
worst_abs_angle_deg
post_clearing_status
```

Sequence currents are diagnostic outputs only; they are not the primary state representation.

---

## 10. Validation Checks

Required checks:

```text
1. During source-connected SLG interval, verify |Va| is near zero.
2. Verify negative-sequence current is present during the unbalanced fault.
3. Verify torque waveform changes with fault inception angle.
4. Confirm time-step convergence of peak torque.
5. Confirm angle sweep repeatability.
6. Confirm BREAKER_CLEARING_TIME_S = "inf" simulates through T_END.
7. Confirm BREAKER_CLEARING_TIME_S = 0 produces no source-connected SLG interval.
8. Confirm simulation terminates at t_clear for finite clearing time.
```

Engineering sanity checks:

```text
- SLG torque is expected to be less central than LL torque for SCIM mechanical studies.
- SLG current outputs from this model are diagnostic only and shall not be used for protection duty.
- Do not assume post-clearing SLG torque without a defined motor-side ground return path.
```

---

## 11. Hard Stop Conditions

Stop with a clear error if:

```text
1. FAULT_TYPE is not "SLG".
2. User attempts to specify FAULTED_PHASE in Phase 1.
3. User attempts to specify GROUNDING_MODEL or GROUND_IMPEDANCE in Phase 1.
4. D_L = Ls*Lr - Lm^2 is non-positive or near zero.
5. Code attempts to continue SLG after t_clear without a future motor-side ground-return model.
```

---

## 12. Recommended Implementation Order

```text
1. Complete and validate LL first.
2. Reuse existing dq parameter normalization.
3. Reuse existing pre-fault steady-state solver.
4. Reuse abc <-> alpha-beta transforms from LL implementation.
5. Implement source-connected canonical A-G boundary condition.
6. Add finite-clearing-time termination logic.
7. Add angle sweep and SLG summary outputs.
8. Add validation tests.
9. Compare SLG torque severity against LL and three-phase short-circuit results.
```
