# LL Fault Terminal Constraints for SCIM Motor-Side Torque Study

## 0. Purpose

This document defines the **motor-terminal electrical boundary conditions** for a squirrel-cage induction motor (SCIM) line-to-line (L-L) fault torque simulation.

This document does **not** redefine the SCIM magnetic model, torque equation, pre-fault steady-state solver, or RK4 integration method. Those are inherited from the existing dq / alpha-beta transient model documents.

The purpose of this document is to remove ambiguity for the coding agent:

```text
Given a motor terminal B-C line-to-line fault,
what terminal voltages or current constraints must be applied
before and after breaker clearing?
```

---

## 1. Phase-1 Scope

### 1.1 Study type

This is a **SCIM motor-side transient torque study**, not a power-system fault-current study.

The model is intended to estimate electromagnetic torque after a motor-side failure such as:

```text
motor lead B touching motor lead C
terminal box B-C flashover
phase-to-phase insulation failure near the motor terminals
phase-to-phase stator winding fault represented at the terminals
```

Do not use this model for:

```text
breaker interrupting duty
relay pickup setting
cable thermal duty
ground grid duty
utility fault-current contribution
```

Those belong in a power-system study tool.

### 1.2 User inputs

The user-facing fault inputs remain:

```json
{
  "FAULT_TYPE": "LL",
  "BREAKER_CLEARING_TIME_S": 0.1
}
```

The user shall not define:

```text
faulted phase pair
fault inception angle
cable impedance
power-system topology
grounding model
```

### 1.3 Canonical fault selection

Use canonical B-C fault internally:

```text
canonical LL fault = B-C short at the motor terminals
```

Reason:

```text
For a symmetric motor and balanced source, A-B, B-C, and C-A faults are equivalent after a 120-degree phase shift. Since the solver sweeps fault inception angle over 0-360 degrees, explicit faulted phase selection does not add useful Phase-1 worst-case torque information.
```

---

## 2. Electrical Variables and Conventions

Use instantaneous phase voltages at the motor terminals:

```text
v_a(t), v_b(t), v_c(t)
```

Use instantaneous motor phase currents entering the motor terminals:

```text
i_a(t), i_b(t), i_c(t)
```

Use amplitude-invariant Clarke transform consistent with the existing dq model:

```text
v_alpha = (2/3)*(v_a - 0.5*v_b - 0.5*v_c)
v_beta  = (2/3)*(sqrt(3)/2)*(v_b - v_c)
v_s     = v_alpha + j*v_beta
```

For balanced three-wire quantities, inverse Clarke is:

```text
v_a = v_alpha
v_b = -0.5*v_alpha + (sqrt(3)/2)*v_beta
v_c = -0.5*v_alpha - (sqrt(3)/2)*v_beta
```

For unbalanced faults, do not assume the zero-sequence component is automatically zero unless the boundary condition explicitly enforces a three-wire balanced condition.

---

## 3. Pre-Fault State

Before the fault:

```text
t < 0
```

The motor is in balanced steady-state operation.

The existing dq pre-fault solver provides:

```text
psi_s0
i_s0
psi_r0
i_r0
omega_m0
```

For each internally swept fault inception angle `theta0`, rotate the pre-fault phasor/state consistently before applying the fault.

Flux continuity at fault inception:

```text
psi_s(0+) = psi_s(0-)
psi_r(0+) = psi_r(0-)
omega_m(0+) = omega_m(0-)
```

Terminal voltages change according to the fault boundary condition at `t = 0`.

---

## 4. Source-Connected B-C Fault Boundary

### 4.1 Time interval

Applicable for:

```text
0 <= t < t_clear
```

where:

```text
t_clear = BREAKER_CLEARING_TIME_S
```

If:

```text
BREAKER_CLEARING_TIME_S = 0
```

this interval has zero duration.

If:

```text
BREAKER_CLEARING_TIME_S = "inf"
```

this interval lasts until the end of the simulation.

### 4.2 Defining fault condition

Canonical B-C short means:

```text
v_b(t) = v_c(t)
```

or:

```text
v_bc(t) = v_b(t) - v_c(t) = 0
```

### 4.3 Reduced terminal-voltage boundary for Phase 1

For a SCIM torque study, use a reduced stiff balanced source boundary:

```text
v_a(t) = e_a(t)
v_b(t) = v_c(t) = 0.5*(e_b(t) + e_c(t))
```

where:

```text
e_a(t) = Vpk*cos(omega_s*t + theta0)
e_b(t) = Vpk*cos(omega_s*t + theta0 - 2*pi/3)
e_c(t) = Vpk*cos(omega_s*t + theta0 + 2*pi/3)
```

This boundary enforces the B-C short:

```text
v_b(t) - v_c(t) = 0
```

For a balanced source, the B-C common-node voltage simplifies to:

```text
v_b(t) = v_c(t) = -0.5*e_a(t)
```

because:

```text
e_a(t) + e_b(t) + e_c(t) = 0
```

Therefore the source-connected canonical B-C fault terminal voltages can be written as:

```text
v_a(t) = e_a(t)
v_b(t) = -0.5*e_a(t)
v_c(t) = -0.5*e_a(t)
```

### 4.4 Resulting alpha-beta voltage

Using the amplitude-invariant Clarke transform:

```text
v_alpha = e_a(t)
v_beta  = 0
v_s     = e_a(t) + j*0
```

This is a useful coding simplification for the canonical B-C source-connected fault.

Important:

```text
This does not mean the motor is balanced.
It only means the imposed terminal voltage vector lies on the alpha axis for the canonical B-C fault.
The resulting current and torque can still be unbalanced/transient because the terminal boundary changed abruptly and the motor states are inherited from pre-fault operation at theta0.
```

### 4.5 Engineering limitation of this source-connected boundary

An ideal voltage source with zero source impedance directly shorted line-to-line is physically singular for source current.

The reduced boundary above should be interpreted as:

```text
balanced stiff source with symmetric source impedance,
used only to impose motor terminal voltage for torque calculation,
not to calculate source fault current.
```

The model shall label current outputs during this interval as diagnostic motor-current quantities only.

Do not report source current duty from this simulation.

---

## 5. Source-Disconnected B-C Fault Boundary

### 5.1 Time interval

Applicable for:

```text
t >= t_clear
```

only when:

```text
0 <= t_clear < T_END
```

After breaker clearing:

```text
source is open
motor remains connected to the B-C short on the motor side
```

### 5.2 Physical terminal constraints

Breaker opens all three phases upstream. The motor-side B and C terminals remain shorted together.

Physical constraints:

```text
v_b(t) = v_c(t)
i_a_terminal(t) = 0
```

For the B-C shorted pair, the external current into the shorted node must also satisfy:

```text
i_b_terminal(t) + i_c_terminal(t) = 0
```

if the B-C short is the only external connection remaining on the motor side.

Thus, after source disconnection, the physically complete motor-terminal constraints are:

```text
v_b = v_c
i_a = 0
i_b + i_c = 0
```

For a three-wire motor terminal set:

```text
i_a + i_b + i_c = 0
```

so `i_a = 0` and `i_b + i_c = 0` are not independent. In practice, enforce:

```text
v_b = v_c
i_a = 0
```

and monitor:

```text
i_a + i_b + i_c
```

as a numerical consistency check.

### 5.3 Critical implementation warning

The existing positive-sequence dq/alpha-beta voltage-driven model does not automatically enforce open-phase current constraints after the breaker opens.

In the existing flux-state model:

```text
i_s = f(psi_s, psi_r)
```

current is recovered algebraically from flux states. At an ideal topology change, enforcing:

```text
i_a = 0
```

may require a voltage impulse or a constrained differential-algebraic formulation.

Therefore, the coding agent must not simply impose:

```text
v_b = v_c
```

inside the existing dq model and assume the post-clearing open-phase condition is satisfied.

### 5.4 Required post-clearing implementation approach

For rigorous post-clearing L-L torque, implement one of the following approaches.

#### Preferred approach: phase-domain constrained stator model

Use a phase-domain stator representation with terminal constraints. The state equations must allow phase terminal constraints such as:

```text
v_b = v_c
i_a = 0
```

to be solved consistently.

The coding agent may still transform results to alpha-beta for torque calculation, but the boundary solve must be phase-domain.

Required elements:

```text
phase-domain stator flux states or equivalent constrained alpha-beta-zero formulation
terminal constraint equations
consistent current recovery under open-phase condition
constrained derivative solve at every RK4 substep or use of a DAE/integration method suitable for constraints
```

#### Acceptable temporary approach: stop at breaker clearing

If the constrained post-clearing solver is not yet implemented, then for finite breaker clearing time:

```text
simulate only 0 <= t <= t_clear
terminate at t_clear
report post-clearing LL motor-only continuation as not implemented
```

The output metadata shall include:

```text
POST_CLEARING_LL_MODEL = not implemented; simulation terminated at breaker clearing
```

Do not substitute the existing three-phase short-circuit boundary `v_s = 0` for the post-clearing L-L condition.

---

## 6. Relationship to Existing dq Model

### 6.1 What can be reused

The following can be reused directly:

```text
parameter normalization
pre-fault steady-state solution
flux-current inversion for balanced alpha-beta machine model
RK4 framework
torque calculation from alpha-beta flux/current
angle sweep infrastructure
CSV and summary output structure
```

### 6.2 What cannot be assumed without review

The following cannot be assumed for post-clearing L-L:

```text
v_s = 0
positive-sequence-only behavior
current continuity without constraint handling
zero-sequence-free behavior for all connection types
source current from ideal source boundary
```

---

## 7. Recommended Phase-1 Coding Plan

### Step 1 - Source-connected LL only

Implement source-connected canonical B-C fault using:

```text
v_s(t) = e_a(t) + j*0
```

for:

```text
0 <= t < min(t_clear, T_END)
```

Then perform angle sweep and report torque up to `t_clear`.

This is the lowest-risk first implementation.

### Step 2 - Validate source-connected LL

Validation checks:

```text
1. v_b - v_c = 0 within numerical tolerance.
2. alpha-beta voltage equals e_a + j*0 for canonical B-C fault.
3. torque changes with fault inception angle.
4. current and torque waveforms differ from three-phase short result.
5. time-step convergence is acceptable.
```

### Step 3 - Decide whether post-clearing LL is required

If source-connected LL torque already captures the peak torque before typical breaker clearing, post-clearing LL may be lower priority.

If post-clearing LL torque is required, proceed to Step 4.

### Step 4 - Implement constrained post-clearing LL

Implement the phase-domain constrained terminal model satisfying:

```text
v_b = v_c
i_a = 0
```

This should be treated as a separate implementation task and test suite.

---

## 8. Outputs Required for Constraint Validation

For debugging and validation, output the following for every saved time step:

```text
time_s
Va_V
Vb_V
Vc_V
Vbc_V
Ialpha_A
Ibeta_A
Ia_A
Ib_A
Ic_A
Ia_plus_Ib_plus_Ic_A
torque_Nm
torque_percent
source_connected_bool
constraint_residual_Vbc
constraint_residual_Ia
```

For source-connected-only implementation, `constraint_residual_Ia` may be blank or marked not applicable.

---

## 9. Hard Stop Rules

Stop with a clear error if:

```text
1. FAULT_TYPE is not LL.
2. A user attempts to select faulted phase pair in Phase 1.
3. BREAKER_CLEARING_TIME_S is finite and post-clearing simulation is requested, but constrained post-clearing solver is not implemented.
4. The source-connected voltage boundary does not satisfy |Vb - Vc| tolerance.
5. The post-clearing constrained solver does not satisfy |Vb - Vc| and |Ia| tolerances.
```

---

## 10. Engineering Conclusion

For a SCIM-focused motor-side L-L torque study:

```text
source-connected LL can be implemented with a reduced canonical voltage boundary:
    v_s = e_a + j*0

post-clearing LL is a constrained open-phase / B-C short problem:
    v_b = v_c
    i_a = 0
```

The source-connected part can be implemented using the existing dq/alpha-beta model.

The post-clearing part should not be approximated as a three-phase short. It requires a constrained phase-domain terminal model or should be explicitly deferred.
