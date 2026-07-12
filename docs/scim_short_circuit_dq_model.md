# SCIM Three-Phase Terminal Short Circuit - dq / State-Space Workflow

## Purpose

Compute transient electromagnetic torque after a balanced three-phase bolted terminal short circuit using a stationary-reference-frame dq/space-vector induction-machine model.

This workflow is preferred when first-cycle torque shape, point-on-wave, stator dc offset, initial flux angle, or speed dynamics matter.

---

## 1. Implementation

The dq model is implemented in `scim_calc/dq.py`, function `run_dq_simulation()`.
It is invoked automatically by the top-level entry point `run.py` alongside the quick calculation.

### Inputs

Motor parameters are read from `input.jsonc` at the project root.
See that file for the full list of per-phase, stator-referred parameters.

| Symbol | `input.jsonc` key | Meaning | Unit |
|---|---|---|---|
| $R_s$ | `Rs` | stator resistance per phase | $\Omega$ |
| $R_r$ | `Rr` | rotor resistance referred to stator | $\Omega$ |
| $X_s$ | `Xs` | stator leakage reactance | $\Omega$ at $f$ |
| $X_r$ | `Xr` | rotor leakage reactance referred to stator | $\Omega$ at $f$ |
| $X_m$ | `Xm` | magnetizing reactance | $\Omega$ at $f$ |
| $GD^2_{\text{rotor}}$ | `ROTOR_GD2` | rotor flywheel effect | kg·m$^2$ |
| $WK^2_{\text{rotor}}$ | `ROTOR_WK2` | rotor flywheel effect (US) | lb·ft$^2$ |
| $s$ | `slip` | pre-fault slip | per unit |
| $V_{LL}$ | `V_LL` | line-line RMS voltage | V |
| | `CONNECTION` | stator connection | `"wye"` or `"delta"` |
| $f$ | `f` | electrical frequency | Hz |
| $p$ | `poles` | total number of poles | 2, 4, 6… |
| | `INITIAL_VOLTAGE_ANGLE_DEG` | fault inception angle | degrees |
| | `USE_SPEED_DYNAMICS` | enable mechanical equation | boolean |

> **Note:** The code accepts **either** reactances ($X_s, X_r, X_m$) **or** inductances ($L_{\ell s}, L_{\ell r}, L_m$). When inductances are given they are used directly; otherwise reactances are converted via $L = X / \omega_s$.

### Derived quantities

$$
\omega_s=2\pi f,\qquad
L_{\ell s}=\frac{X_s}{\omega_s},\qquad
L_{\ell r}=\frac{X_r}{\omega_s},\qquad
L_m=\frac{X_m}{\omega_s}
$$

$$
L_s=L_{\ell s}+L_m,\qquad
L_r=L_{\ell r}+L_m,\qquad
\Delta=L_sL_r-L_m^2
$$

---

## 2. Stationary-frame space-vector model

Complex stationary-frame space vectors:

$$
\psi_s=\psi_{ds}+j\psi_{qs},\qquad
\psi_r=\psi_{dr}+j\psi_{qr}
$$

Flux-current relations and their inversion (evaluated at every time step):

$$
\begin{aligned}
\psi_s &= L_s i_s+L_m i_r \\
\psi_r &= L_m i_s+L_r i_r
\end{aligned}
\quad\Longrightarrow\quad
\begin{aligned}
i_s &= \frac{L_r\psi_s-L_m\psi_r}{\Delta} \\
i_r &= \frac{-L_m\psi_s+L_s\psi_r}{\Delta}
\end{aligned}
$$

Stator and rotor voltage equations (squirrel-cage rotor, $v_r=0$):

$$
\frac{d\psi_s}{dt}=v_s-R_s i_s,\qquad
\frac{d\psi_r}{dt}=j\omega_r^e\psi_r-R_r i_r,\qquad
\omega_r^e=p\omega_m
$$

For a bolted three-phase terminal short circuit at $t=0$:

$$
v_s(t\ge 0)=0
$$

---

## 3. Torque equation

Using stationary dq quantities:

$$
T_e=\frac{3}{2}p\left(\psi_{ds}i_{qs}-\psi_{qs}i_{ds}\right)
     =-\frac{3}{2}p\operatorname{Im}\left\{\psi_s i_s^*\right\}
$$

The normalized output is:

$$
T_{\%}(t)=100\frac{T_e(t)}{T_{nom}}
$$

The sign convention is chosen so that pre-fault torque is positive.

---

## 4. Initial steady state from slip

Phase voltage depends on winding connection:

$$
V_{ph}=
\begin{cases}
\dfrac{V_{LL}}{\sqrt{3}}, & \text{wye}\\[8pt]
V_{LL}, & \text{delta}
\end{cases}
$$

The peak phase voltage $\sqrt{2}\,V_{ph}$ is used as the initial space-vector magnitude.

Solve the synchronous steady-state phasor equations:

$$
\begin{aligned}
V_s &= \left(R_s+j\omega_sL_s\right)I_{s0}+j\omega_sL_m I_{r0} \\
0   &= js\omega_sL_m I_{s0}+\left(R_r+js\omega_sL_r\right)I_{r0}
\end{aligned}
$$

Initial fluxes and speed:

$$
\psi_{s0}=L_s I_{s0}+L_m I_{r0},\qquad
\psi_{r0}=L_m I_{s0}+L_r I_{r0},\qquad
\omega_{m0}=\frac{(1-s)\omega_s}{p}
$$

---

## 5. Numerical integration

1. Convert reactances to inductances (handled by `config.py`).
2. Solve the pre-fault phasor equations from slip (`circuit.py:solve_prefault()`).
3. Set initial flux states $\psi_{s0}$, $\psi_{r0}$.
4. At $t\ge 0$, set $v_s=0$.
5. Integrate the flux ODE system with classical RK4 (`scim_calc/dq.py`).
6. At each step, compute $i_s$, $i_r$, and $T_e$.
7. Save and plot results.

Optional speed dynamics:

$$
J\frac{d\omega_m}{dt}=T_e-T_L
$$

For first-cycle electrical torque, constant speed is normally used (`USE_SPEED_DYNAMICS = false`).

---

## 6. DC component handling

No separate dc-offset equation is required. The stator dc/asymmetrical component appears naturally because flux linkages are continuous at fault inception while terminal voltage is forced to zero.

Use `INITIAL_VOLTAGE_ANGLE_DEG` in `input.jsonc` to vary the fault inception angle. The `scim_calc/sweep.py` module provides a two-phase coarse/refine sweep to find the worst-case angle automatically.

---

## 7. Script output

Running `python run.py` at the project root:

1. reads motor parameters from `input.jsonc`;
2. runs both the dq model and the quick calculation;
3. prints a comparison summary;
4. saves CSV files and an overlay PNG to the `data/` directory.

The dq model alone can also be called programmatically:

```python
from scim_calc.config import load_jsonc, normalize_params
from scim_calc.dq import run_dq_simulation

cfg = load_jsonc("input.jsonc")
p = normalize_params(cfg)
d = run_dq_simulation(p)
```
