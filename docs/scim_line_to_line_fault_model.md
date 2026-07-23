# SCIM Line-to-Line Fault Torque - Implementation Instructions

## 0. Purpose

Implement a **SCIM-focused motor terminal line-to-line (L-L) fault torque simulation**.

This is not a power-system short-circuit-current study. Power-system current duty, breaker interrupting current, relay pickup, and cable duty should remain outside this model and can be handled by ETAP/SKM/EasyPower or similar tools.

The purpose here is to calculate:

- transient electromagnetic torque,
- worst positive torque,
- worst negative torque,
- worst absolute torque,
- torque waveform after a motor-side phase-to-phase insulation failure.

Typical real-life failures represented by this model:

- motor lead phase-to-phase insulation failure,
- terminal box phase-to-phase flashover,
- winding insulation failure causing two phase circuits to touch,
- contamination or tracking between phase terminals.

---

## 1. User Inputs

For Phase 1, user-facing fault inputs shall be limited to:

```json
{
  "FAULT_TYPE": "LL",
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
FAULTED_PHASES
FAULT_INCEPTION_ANGLE_DEG
CABLE_IMPEDANCE
GROUNDING_MODEL
SOURCE_NETWORK_MODEL
```

### 1.1 Breaker clearing time convention

```text
BREAKER_CLEARING_TIME_S = 0
    source is disconnected immediately at fault inception

0 < BREAKER_CLEARING_TIME_S < T_END
    source remains connected from t = 0 to t = t_clear
    source is disconnected after t_clear

BREAKER_CLEARING_TIME_S = "inf"
    source remains connected for the entire simulation
```

### 1.2 Fault inception angle convention

The fault inception angle is **not** a normal user input. The solver shall sweep it internally.

For each breaker clearing time, the solver shall report:

```text
worst positive torque and associated angle
worst negative torque and associated angle
worst absolute torque and associated angle
```

---

## 2. Main Phase-1 Assumptions

### 2.1 Motor-side study, not system study

Assume the fault is at the motor terminals or in the motor leads close enough that cable impedance is neglected.

```text
FAULT_LOCATION = motor_terminal
CABLE_IMPEDANCE = 0
```

This is intentional. Adding cable impedance usually reduces fault severity and is not needed for a conservative motor torque study.

### 2.2 Canonical phase selection

Use a canonical B-C fault internally:

```text
canonical LL fault = B-C terminal short
```

Do not ask the user to select A-B, B-C, or C-A. For a symmetric motor and balanced source, the three choices are equivalent after a 120-degree phase shift. Since the solver sweeps inception angle over 0-360 degrees, explicit phase-pair selection does not add useful worst-case torque information in Phase 1.

Metadata shall include:

```text
LL_FAULT_PHASE_ASSUMPTION = canonical B-C fault; phase pair not user-selected
```

### 2.3 Source representation for torque study

During the source-connected interval, use a balanced ideal voltage reference to define terminal voltage boundary conditions under the fault.

Do not calculate source current duty from this model. If source current duty is needed, use a power-system fault-current program.

The torque model only needs terminal voltage excitation applied to the motor model.

### 2.4 Motor connection

The existing input `CONNECTION = "wye" | "delta"` shall continue to define the conversion from line-line voltage to phase voltage for the pre-fault equivalent machine model.

For Phase 1 unbalanced fault torque simulation, the model uses the same stator-referred equivalent phase model as the existing dq implementation. The implementation shall clearly state that detailed delta internal circulating branch currents are not separately resolved unless a future full phase-winding model is added.

---

## 3. Reused SCIM State Model

Reuse the existing stationary alpha-beta / dq SCIM state equations.

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

Use `D_L` or `det_L`, not only `Delta`, to avoid confusion with delta connection.

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

Default for first-cycle torque studies:

```text
USE_SPEED_DYNAMICS = false
```

---

## 4. Pre-Fault Initial Condition

Use the existing balanced pre-fault solution.

Steps:

1. Convert `V_LL` to phase voltage using `CONNECTION`.
2. Solve pre-fault steady-state induction-machine phasor equations at the requested slip.
3. Construct initial `psi_s0`, `psi_r0`, and `omega_m0`.
4. Rotate the initial steady-state voltage/flux/current according to the internal fault inception angle being swept.

At fault inception:

```text
psi_s(0+) = psi_s(0-)
psi_r(0+) = psi_r(0-)
omega_m(0+) = omega_m(0-)
```

Flux states are continuous. Terminal boundary conditions change at t = 0.

---

## 5. Terminal Voltage Boundary Conditions

The line-to-line fault is represented through terminal phase-voltage constraints and then converted to alpha-beta voltage for the machine model.

Use amplitude-invariant Clarke transform consistent with the existing dq model.

### 5.1 Source-connected interval

For canonical B-C fault while the upstream source is still connected:

```text
0 <= t < t_clear
```

Use the following reduced terminal voltage boundary condition:

```text
v_a(t) = e_a(t)
v_b(t) = v_c(t) = 0.5*(e_b(t) + e_c(t))
```

For a balanced source:

```text
e_a(t) = Vpk*cos(omega_s*t + theta0)
e_b(t) = Vpk*cos(omega_s*t + theta0 - 2*pi/3)
e_c(t) = Vpk*cos(omega_s*t + theta0 + 2*pi/3)
```

This enforces:

```text
v_b(t) - v_c(t) = 0
```

and represents a symmetric upstream source with the B and C terminals shorted together. It is a motor-terminal torque boundary condition, not a source-current-duty calculation.

Then compute:

```text
v_alpha, v_beta = Clarke(v_a, v_b, v_c)
v_s = v_alpha + j*v_beta
```

### 5.2 Source-disconnected interval

For canonical B-C fault after breaker clearing:

```text
t >= t_clear
```

The motor-side B-C short remains, but the source is open.

For Phase 1 implementation, use a constrained terminal-voltage model with the following physical constraints:

```text
v_b(t) = v_c(t)
i_a_terminal(t) = 0
```

Interpretation:

```text
phase A lead is open after the breaker opens
phases B and C remain shorted together on the motor side
```

The coding agent shall solve for the terminal voltage vector that satisfies these constraints using the instantaneous machine state.

Recommended implementation path:

```text
1. At an RK4 substep, compute the unconstrained stator current response as a function of unknown terminal voltages.
2. Express phase currents from alpha-beta current using inverse Clarke transform.
3. Solve the two scalar constraints:
       v_b - v_c = 0
       i_a = 0
4. Convert solved v_a, v_b, v_c to v_alpha, v_beta.
5. Use v_s in dpsi_s/dt.
```

Do not silently replace this stage with `v_s = 0`; that would become a three-phase short circuit and is not an L-L fault.

If the post-clearing constrained solve is not implemented in the first coding pass, the code shall stop at `t_clear` and label the result as source-connected-only LL torque.

---

## 6. Torque Calculation

Compute electromagnetic torque from instantaneous alpha-beta quantities:

```text
T_e = (3/2)*pole_pairs*(psi_s_alpha*i_s_beta - psi_s_beta*i_s_alpha)
```

Do not compute final torque as separate positive- and negative-sequence scalar torques unless all cross terms are explicitly included and validated.

The alpha-beta torque expression naturally includes positive/negative sequence interaction and twice-line-frequency torque pulsation.

Normalize using the same sign convention as the existing dq model:

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
4. Re-run final high-resolution simulations for the selected worst angles.
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
1. Determine whether source is connected based on t and t_clear.
2. Compute terminal phase voltages from LL boundary condition.
3. Convert terminal phase voltages to alpha-beta voltage.
4. Recover stator and rotor currents from flux states.
5. Compute derivatives.
```

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
I0_A
torque_Nm
torque_percent
omega_m_rad_per_s
source_connected_bool
```

Write summary output:

```text
fault_type
fault_phase_assumption
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
post_clearing_model_status
```

Sequence currents are diagnostic outputs only; they are not the primary state representation.

---

## 10. Validation Checks

Required checks:

```text
1. During source-connected LL interval, verify |Vb - Vc| is near zero.
2. During source-disconnected LL interval, verify |Vb - Vc| is near zero and Ia is near zero.
3. Confirm negative-sequence current is present during the LL fault.
4. Confirm torque contains unbalanced-fault pulsation components.
5. Confirm time-step convergence of peak torque.
6. Confirm angle sweep repeatability.
7. Confirm BREAKER_CLEARING_TIME_S = "inf" skips the disconnected stage.
8. Confirm BREAKER_CLEARING_TIME_S = 0 skips the source-connected stage.
```

Engineering sanity checks:

```text
- LL torque shall not be assumed lower than three-phase short-circuit torque.
- Worst torque shall be determined by sweep, not by fixed inception angle.
- Fault current values from this model are diagnostic only and shall not be used for breaker duty.
```

---

## 11. Hard Stop Conditions

Stop with a clear error if:

```text
1. FAULT_TYPE is not "LL".
2. User attempts to specify FAULTED_PHASES in Phase 1.
3. D_L = Ls*Lr - Lm^2 is non-positive or near zero.
4. Post-clearing LL run is requested but the constrained source-disconnected solver is not implemented.
5. RK4 boundary-condition solve fails to satisfy voltage/current constraints.
```

---

## 12. Recommended Implementation Order

```text
1. Reuse existing dq parameter normalization.
2. Reuse existing pre-fault steady-state solver.
3. Add abc <-> alpha-beta transforms.
4. Implement source-connected canonical B-C boundary condition.
5. Implement angle sweep and LL summary outputs.
6. Add validation tests for source-connected LL.
7. Implement source-disconnected constrained B-C model.
8. Add validation tests for post-clearing LL.
9. Compare LL worst torque against existing three-phase-short model.
```
