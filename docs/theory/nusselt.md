# Nusselt and Reynolds numbers

Notation: $\langle\cdot\rangle$ is the volume average, $\langle\cdot\rangle_h$ the
horizontal average, and an overline the time average over the averaging window. All
quantities are in {ref}`free-fall units <free-fall-units>`, where $\kappa = 1/\sqrt{Ra\,Pr}$
and $\nu = \sqrt{Pr/Ra}$.

## Definition

The vertical heat flux is advection plus conduction, $J = wT - \kappa\,\partial_z T$. The
Nusselt number is its horizontal and time average divided by the conductive flux of the
same plates with the fluid at rest, which is $\kappa$:

$$
Nu = \frac{\overline{\langle J\rangle_h}}{\kappa}
   = \sqrt{Ra\,Pr}\;\overline{\langle wT\rangle_h} - \partial_z\overline{\langle T\rangle_h}.
$$

In a statistically steady state this is the same at every height: the horizontally and
time-averaged temperature equation reads $\partial_z\overline{\langle J\rangle_h} = 0$.
That makes the profile $Nu(z)$ written by {mod}`rbc.analysis.profiles` a check on
stationarity.

## Three identities

The code computes Nu three ways, from volume averages saved every `A_FREQ`:

$$
Nu \;=\; 1 + \sqrt{Ra\,Pr}\,\overline{\langle wT\rangle}
   \;=\; \overline{\langle|\nabla T|^2\rangle}
   \;=\; 1 + Pr\,\overline{\langle|\boldsymbol\omega|^2\rangle}.
$$

They are equal for a long-time average of an exact solution, so their disagreement measures
two things at once: the window being too short, and the resolution being too coarse.

```{dropdown} Derivation: Nu = 1 + √(Ra Pr) ⟨wT⟩
:class-container: derivation

Average the definition over $0 \le z \le 1$. The conductive part integrates to
$-\int_0^1 \partial_z\overline{\langle T\rangle_h}\,dz = T(0) - T(1) = 1$, and the
advective part is the volume average $\sqrt{Ra\,Pr}\,\overline{\langle wT\rangle}$.
```

```{dropdown} Derivation: Nu = ⟨|∇T|²⟩
:class-container: derivation

Multiply the temperature equation by $T$ and average over the volume:

$$
\frac{d}{dt}\Big\langle\frac{T^2}{2}\Big\rangle
 = -\Big\langle \boldsymbol u\cdot\nabla\frac{T^2}{2}\Big\rangle + \kappa\langle T\nabla^2 T\rangle.
$$

The advective term is a divergence, $\nabla\cdot(\boldsymbol u T^2/2)$, since
$\nabla\cdot\boldsymbol u = 0$; it vanishes because $w = 0$ on the plates and the sides
are periodic. Integrating the last term by parts,
$\langle T\nabla^2T\rangle = \frac{1}{V}\oint T\,\partial_n T\,dS - \langle|\nabla T|^2\rangle$.
On the top plate $T = 0$; on the bottom plate $T = 1$ and the outward normal is
$-\hat{\boldsymbol z}$, so the boundary term is the area average of $-\partial_z T$ at
$z = 0$, the bottom-plate Nusselt number $Nu_\text{bot}(t)$ (the volume has unit height).
Hence

$$
\frac{d}{dt}\Big\langle\frac{T^2}{2}\Big\rangle = \kappa\big(Nu_\text{bot}(t) - \langle|\nabla T|^2\rangle\big).
$$

The time average of the left side vanishes in a steady state, and that of $Nu_\text{bot}$
is $Nu$.
```

(nu-vorticity)=
```{dropdown} Derivation: Nu = 1 + Pr ⟨|ω|²⟩
:class-container: derivation

Dot the momentum equation with $\boldsymbol u$ and average: advection and pressure are
divergences that vanish with no-slip plates and periodic sides, leaving

$$
\frac{d}{dt}\Big\langle\frac{|\boldsymbol u|^2}{2}\Big\rangle = \langle wT\rangle - \nu\langle|\nabla\boldsymbol u|^2\rangle .
$$

For incompressible flow, $|\nabla\boldsymbol u|^2 = |\boldsymbol\omega|^2$ plus a
divergence whose boundary integral vanishes when $\boldsymbol u = 0$ on the plates, so
$\langle|\nabla\boldsymbol u|^2\rangle = \langle|\boldsymbol\omega|^2\rangle$. In a steady
state, $\nu\overline{\langle|\boldsymbol\omega|^2\rangle} = \overline{\langle wT\rangle} = (Nu-1)/\sqrt{Ra\,Pr}$,
and dividing by $\nu = \sqrt{Pr/Ra}$:

$$
\overline{\langle|\boldsymbol\omega|^2\rangle} = \frac{Nu - 1}{\sqrt{Ra\,Pr}}\sqrt{\frac{Ra}{Pr}} = \frac{Nu - 1}{Pr}.
$$
```

### Why they differ over a finite window

The derivations above contain time derivatives that only vanish on average. Over a window
of length $T_w$ they leave exact remainders:

$$
Nu_{|\nabla T|^2} = \overline{Nu_\text{bot}} - \sqrt{Ra\,Pr}\;\frac{\Delta\langle T^2/2\rangle}{T_w},
\qquad
Nu_{|\omega|^2} = Nu_{wT} - \sqrt{Ra\,Pr}\;\frac{\Delta\langle|\boldsymbol u|^2/2\rangle}{T_w},
$$

where $\Delta$ is the change between the ends of the window. The prefactor
$\sqrt{Ra\,Pr}$ can be large, so small fluctuations of the stored thermal or kinetic
energy at the window's ends shift these two estimates noticeably, and the shift shrinks
only as $1/T_w$. The $\langle wT\rangle$ form has no such term, which is why it is the
one reported as `Nu`.

Resolution enters differently: $\langle|\nabla T|^2\rangle$ and above all
$\langle|\boldsymbol\omega|^2\rangle$ weight the smallest scales, so an under-resolved
run overshoots in the vorticity form well before $\langle wT\rangle$ changes. That
difference is reported as `vort_diff_pct` and is a resolution diagnostic; it is not part of
the convergence checks, which use only $\langle wT\rangle$ and $\langle|\nabla T|^2\rangle$
({ref}`which identities are checked <convergence-identities>`). Their spread `qoi_spread_pct` is likewise a check of
resolution once the window is long.

## Plate fluxes

From the horizontally averaged temperature saved in `horizontal_analysis`,

$$
Nu_\text{bot}(t) = -\partial_z\langle T\rangle_h(0,t),\qquad Nu_\text{top}(t) = -\partial_z\langle T\rangle_h(1,t),
$$

evaluated exactly through the Chebyshev interpolant
({func}`rbc.analysis.nusselt.wall_gradient_weights`). Integrating the temperature equation
over the volume gives

$$
\frac{d\langle T\rangle}{dt} = \kappa\,\big(Nu_\text{bot} - Nu_\text{top}\big),
$$

so in a steady state the two plates carry the same mean flux, and a persistent difference
means the mean temperature is still drifting. This is one of the transient checks in
{doc}`convergence`.

## Reynolds number and rms vorticity

$$
Re = \frac{U_\text{rms}H}{\nu} = \sqrt{\frac{Ra}{Pr}}\,\overline{\langle|\boldsymbol u|^2\rangle}^{1/2},
\qquad
\omega_\text{rms} = \overline{\langle|\boldsymbol\omega|^2\rangle}^{1/2} = \sqrt{\frac{Nu-1}{Pr}}.
$$

`Re` comes from the volume average `avg_K`; `Re_rms` from the time-averaged profiles. The
exact $\omega_\text{rms}$ sets the fixed colour scale of the vorticity visualizations.

## Why Nu − 1

$Nu - 1$ is the convective part of the heat flux, $\sqrt{Ra\,Pr}\,\langle wT\rangle$; it is
what the identities above and the scaling theories refer to. At moderate Nu the
difference matters for exponents:

$$
\frac{d\ln Nu}{d\ln Ra} = \frac{Nu-1}{Nu}\,\frac{d\ln(Nu-1)}{d\ln Ra},
$$

so at $Nu = 3$ the local exponent of Nu is only two thirds of that of $Nu - 1$. The
exponent analysis and the fixed-parameter figure use $Nu - 1$ by default.
