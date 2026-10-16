# Equations and numerics

## The physical problem

A layer of fluid of depth $H$ between two horizontal plates, the bottom one at
temperature $T_\text{bot}$ and the top one cooler by $\Delta T$, in gravity $g$. In the
Boussinesq approximation the density varies only through the buoyancy term,
$\rho = \rho_0[1 - \alpha_T(T - T_0)]$, and the flow is incompressible:

$$
\partial_t \boldsymbol u + \boldsymbol u\cdot\nabla\boldsymbol u
  = -\nabla p + \nu\nabla^2\boldsymbol u + g\alpha_T T\,\hat{\boldsymbol z},\qquad
\partial_t T + \boldsymbol u\cdot\nabla T = \kappa\nabla^2 T,\qquad
\nabla\cdot\boldsymbol u = 0,
$$

with kinematic viscosity $\nu$ and thermal diffusivity $\kappa$ (the constant part of the
buoyancy is absorbed into $p$). Two numbers control the flow:

$$
Ra = \frac{g\alpha_T\Delta T H^3}{\nu\kappa},\qquad Pr = \frac{\nu}{\kappa}.
$$

The domain is periodic in $x$ and $y$ with $L_x = L_y = \Gamma H$, $\Gamma = 2$. The plates
are no-slip and isothermal.

(free-fall-units)=
## Free-fall units

Lengths are measured in $H$, temperature (relative to the top plate) in $\Delta T$, and
velocity in the free-fall velocity $U_f = \sqrt{g\alpha_T\Delta T H}$, so time is in
$H/U_f$. The equations become

$$
\partial_t \boldsymbol u + \boldsymbol u\cdot\nabla\boldsymbol u
  = -\nabla p + \sqrt{\frac{Pr}{Ra}}\,\nabla^2\boldsymbol u + T\,\hat{\boldsymbol z},\qquad
\partial_t T + \boldsymbol u\cdot\nabla T = \frac{1}{\sqrt{Ra\,Pr}}\nabla^2 T,
$$

with $T = 1$ at $z = 0$, $T = 0$ at $z = 1$, and $\boldsymbol u = 0$ on both plates. In the
code, `nu = (Ra/Pr)**(-1/2)` and `kappa = (Ra*Pr)**(-1/2)`.

```{dropdown} Derivation
:class-container: derivation

Write $\boldsymbol x = H\boldsymbol x'$, $\boldsymbol u = U_f\boldsymbol u'$,
$t = (H/U_f)t'$, $T - T_\text{top} = \Delta T\,T'$, $p = U_f^2 p'$. The momentum equation,
divided by $U_f^2/H$, becomes

$$
\partial_{t'}\boldsymbol u' + \boldsymbol u'\cdot\nabla'\boldsymbol u'
 = -\nabla'p' + \frac{\nu}{U_f H}\nabla'^2\boldsymbol u'
 + \frac{g\alpha_T\Delta T H}{U_f^2}T'\hat{\boldsymbol z}.
$$

The buoyancy coefficient is 1 by the choice of $U_f$. For the viscous one,
$U_f H = \sqrt{g\alpha_T\Delta T H^3}$, and $Ra\,\nu\kappa = g\alpha_T\Delta T H^3$, so
$U_f H = \sqrt{Ra\,\nu\kappa}$ and

$$
\frac{\nu}{U_f H} = \frac{\nu}{\sqrt{Ra\,\nu\kappa}} = \sqrt{\frac{\nu}{\kappa Ra}} = \sqrt{\frac{Pr}{Ra}}.
$$

The same steps in the temperature equation give $\kappa/(U_f H) = \sqrt{\kappa/(\nu Ra)}
= 1/\sqrt{Ra\,Pr}$.
```

Why these units: the velocity scale is set by the buoyancy, not by a diffusivity, so
velocities, timesteps and run lengths stay of order one across the whole (Ra, Pr) range.
In thermal-diffusion units, by contrast, the unit of time would be $\sqrt{Ra\,Pr}$
free-fall times, which ranges from about 10 ($Ra = 10^5$, $Pr = 10^{-3}$) to $10^4$
($Ra = 10^8$, $Pr = 1$) over the runs here. A free-fall time is roughly
the time a plume needs to cross the layer.

## Conduction and the perturbation

Without motion the temperature is linear, $T_c = 1 - z$. The code evolves the departure
$b = T - (1 - z)$, which vanishes on both plates:

$$
\partial_t b + \boldsymbol u\cdot\nabla(1 - z) - \kappa\nabla^2 b = -\boldsymbol u\cdot\nabla b,
\qquad b(0) = b(1) = 0,
$$

and the buoyancy is written $b\,\hat{\boldsymbol z}$ (the part $(1-z)\hat{\boldsymbol z}$ is
a gradient, absorbed into the pressure). The momentum equation is used in rotational form,
$\boldsymbol u\cdot\nabla\boldsymbol u = \nabla(|\boldsymbol u|^2/2) - \boldsymbol u\times\boldsymbol\omega$,
with $\boldsymbol\omega = \nabla\times\boldsymbol u$ and the gradient absorbed into $p$;
that is the `cross(u, omega)` on the right-hand side of the code.

## Spatial discretization

| Direction | Basis | Modes | Grid |
|---|---|---|---|
| $x$, $y$ | real Fourier (`RealFourier`) | $N_x = N_y = 2N_z$ | uniform |
| $z$ | Chebyshev polynomials (`ChebyshevT`) | $N_z$ | Gauss–Chebyshev points, clustered at the plates |

Nonlinear products are evaluated on a grid 3/2 times finer in each direction
({term}`dealiasing`), which removes aliasing errors of quadratic products.
{mod}`rbc.spectra.core` documents exactly how Dedalus stores the coefficients.

### Boundary conditions: the tau method

A Chebyshev series of degree $N_z - 1$ cannot satisfy both a second-order differential
equation and two boundary conditions exactly. The tau method makes room: each second-order
equation gets two extra unknown functions of $(x, y)$, the *tau variables*, multiplying
fixed high-degree polynomials (`lift`), and the two boundary conditions determine them. The
error this introduces is confined to the highest modes, and is small when the solution is
resolved. In the code, `tau_u1`, `tau_u2`, `tau_b1`, `tau_b2` play this role for velocity
and temperature, and `tau_p` for the pressure gauge.

```{dropdown} The horizontally uniform mode
:class-container: derivation

For $k_x = k_y = 0$, continuity reduces to $\partial_z \bar w = 0$, and $\bar w(1) = 0$
then forces $\bar w \equiv 0$, so the condition $\bar w(0) = 0$ is redundant. Meanwhile
the mean pressure is defined only up to a constant. The code therefore replaces
$\bar w(0) = 0$ by the gauge $\bar p(0) = 0$ and sets the vertical component of the first
velocity tau to zero for this mode; the horizontal components keep their no-slip
conditions. These are the four equations with `condition="nx == 0 and ny == 0"`.
```

(timestepping)=
## Timestepping

The diffusion, pressure and buoyancy terms are treated implicitly and the advection
explicitly, with Dedalus' second-order, two-stage implicit–explicit Runge–Kutta scheme
**RK222**. The two-step multistep scheme SBDF2 was tried and gave oscillatory spectra and
needed smaller timesteps.

The timestep follows a {term}`CFL condition`: $\Delta t = C\min(\Delta x/|u|)$ over the
grid, with `CFL_SAFETY` $= C = 1/3$, recomputed every `CFL_CADENCE` iterations and changed
only when it moves by more than `CFL_THRESHOLD` (10%), since every change refactorizes the
implicit matrices. The largest timestep allowed is 0.1.

## Initial condition

A fresh run starts with $\boldsymbol u = 0$ and
$b = 10^{-5}\,\xi(\boldsymbol x)\,z(1-z)$, with $\xi$ Gaussian noise (fixed seed) low-pass
filtered to the lowest quarter of the modes in each direction. The noise seeds the convective instability;
the factor $z(1-z)$ makes it satisfy the boundary conditions. Convection grows
exponentially, overshoots, and settles into a statistically steady state; that transient
is removed from all averages ({ref}`transient`).

## Choices

**Horizontally periodic.** No side walls means no side-wall boundary layers, which would
need resolving and would change Nu, and Fourier bases, which are the most efficient. The
cost is that the box size can constrain the largest flow structures.

**Γ = 2.** A box twice as wide as it is tall, as in many periodic-domain studies.

```{admonition} To fill in
:class: todo-note

The reasons for choosing Γ = 2 (and whether larger boxes were tested).
```
