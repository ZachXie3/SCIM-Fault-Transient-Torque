# SCIM Three-Phase Terminal Short Circuit — dq / State-Space Workflow

## Purpose

Compute transient electromagnetic torque after a balanced three-phase bolted terminal short circuit using a stationary-reference-frame dq/space-vector induction-machine model.

This workflow is preferred when first-cycle torque shape, point-on-wave, stator dc offset, initial flux angle, or speed dynamics matter.

---

## 1. Inputs and units

Use per-phase, stator-referred parameters.

| Symbol | Python variable | Meaning | Unit / note |
|---|---|---|---|
| $R_s$ | `Rs` | stator resistance per phase | $\Omega$ |
| $R_r$ | `Rr` | rotor resistance referred to stator | $\Omega$ |
| $X_s$ | `Xs` | stator leakage reactance | $\Omega$ at frequency $f$ |
| $X_r$ | `Xr` | rotor leakage reactance referred to stator | $\Omega$ at frequency $f$ |
| $X_m$ | `Xm` | magnetizing reactance | $\Omega$ at frequency $f$ |
| $J$ | `J` | total inertia referred to motor shaft | kg·m$^2$ |
| $s$ | `slip` | pre-fault slip | per unit |
| $V_{LL}$ | `V_LL` | pre-fault line-line RMS voltage | V |
| $f$ | `f` | electrical frequency | Hz |
| $p$ | `pole_pairs` | pole pairs | dimensionless |

Electrical angular frequency:

$$
\omega_s=2\pi f
$$

Reactance-to-inductance conversion:

$$
L_{\ell s}=\frac{X_s}{\omega_s},\qquad
L_{\ell r}=\frac{X_r}{\omega_s},\qquad
L_m=\frac{X_m}{\omega_s}
$$

Total stator and rotor inductances:

$$
L_s=L_{\ell s}+L_m
$$

$$
L_r=L_{\ell r}+L_m
$$

Leakage determinant:

$$
\Delta=L_sL_r-L_m^2
$$

---

## 2. Stationary-frame space-vector model

Use complex stationary-frame space vectors:

$$
\psi_s=\psi_{ds}+j\psi_{qs}
$$

$$
\psi_r=\psi_{dr}+j\psi_{qr}
$$

Flux-current relations:

$$
\psi_s=L_s i_s+L_m i_r
$$

$$
\psi_r=L_m i_s+L_r i_r
$$

Inverted current equations used at every time step:

$$
i_s=\frac{L_r\psi_s-L_m\psi_r}{\Delta}
$$

$$
i_r=\frac{-L_m\psi_s+L_s\psi_r}{\Delta}
$$

Stator voltage equation:

$$
\frac{d\psi_s}{dt}=v_s-R_s i_s
$$

Rotor voltage equation for a squirrel-cage rotor:

$$
\frac{d\psi_r}{dt}=j\omega_r^e\psi_r-R_r i_r
$$

where the rotor electrical speed is:

$$
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
$$

Equivalent complex form used in the script:

$$
T_e=-\frac{3}{2}p\operatorname{Im}\left\{\psi_s i_s^*\right\}
$$

The normalized output is:

$$
T_{\%}(t)=100\frac{T_e(t)}{T_{nom}}
$$

The script automatically chooses the sign convention so that the pre-fault torque is positive on the plot.

---

## 4. Initial steady state from slip

The dynamic model uses peak phase voltage:

$$
V_{\phi,pk}=\sqrt{2}\frac{V_{LL}}{\sqrt{3}}
$$

With phasors in the synchronous steady-state frame, solve:

$$
V_s=\left(R_s+j\omega_sL_s\right)I_{s0}+j\omega_sL_m I_{r0}
$$

$$
0=js\omega_sL_m I_{s0}+\left(R_r+js\omega_sL_r\right)I_{r0}
$$

Then compute initial fluxes:

$$
\psi_{s0}=L_s I_{s0}+L_m I_{r0}
$$

$$
\psi_{r0}=L_m I_{s0}+L_r I_{r0}
$$

Initial rotor mechanical speed:

$$
\omega_{m0}=\frac{(1-s)\omega_s}{p}
$$

Nominal/pre-fault torque:

$$
T_{nom}=T_e(0^-)
$$

---

## 5. Numerical integration workflow

1. Convert reactances to inductances.
2. Solve the pre-fault phasor equations from slip.
3. Set the initial flux states $\psi_{s0}$ and $\psi_{r0}$.
4. At $t\ge 0$, set $v_s=0$.
5. Integrate the flux equations with RK4.
6. At each step, compute $i_s$, $i_r$, and $T_e$.
7. Plot $T_{fault}/T_{nominal}$ in percent over $0$ to $0.2$ s.

Optional speed dynamics:

$$
J\frac{d\omega_m}{dt}=T_e-T_L
$$

For first-cycle electrical torque, constant speed is normally used. Enable speed dynamics only when the 0.2 s window and inertia/load interaction are important.

---

## 6. DC component handling

No separate dc-offset equation is required in the dq method. The stator dc/asymmetrical component appears naturally because flux linkages are continuous at fault inception while terminal voltage is forced to zero.

The actual first-cycle waveform depends on the point-on-wave angle. Use `INITIAL_VOLTAGE_ANGLE_DEG` in the script to sweep the fault inception angle.

---

## 7. Script output

The companion Python script:

1. declares user-editable parameters at the top;
2. computes initial steady-state currents and fluxes from slip;
3. integrates the short-circuit dynamic model over $0$ to $0.2$ s;
4. plots $T_{fault}/T_{nominal}$ in percent;
5. saves a CSV and PNG file.
