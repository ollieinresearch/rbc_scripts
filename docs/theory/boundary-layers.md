# Boundary layers

Near each plate the mean temperature changes from the plate value to the bulk value of
about 1/2 over a thin **thermal boundary layer**, and the horizontal velocity rises from
zero to its bulk magnitude over a **viscous boundary layer**. They are the hardest features
to resolve in $z$, and at small Pr the thermal layer is much thicker than the viscous one.

## Definitions

{mod}`rbc.analysis.profiles` measures, from the time-averaged horizontal-mean profiles
(averaged over the bottom and top plates):

| Symbol | Column | Definition |
|---|---|---|
| $\delta_T$ | `delta_T` | where the wall tangent of $\overline T$ reaches the bulk value 1/2: $1/(2Nu)$ |
| $\delta_{T,\text{rms}}$ | `delta_T_rms` | distance from the wall to the peak of $\theta_\text{rms}$ |
| $\delta_u$ | `delta_u` | distance from the wall to the peak of the rms horizontal velocity |
| $\delta_u^\text{slope}$ | `delta_u_slope` | where the wall tangent of the rms horizontal velocity reaches its peak value |

```{dropdown} Derivation: δ_T = 1/(2 Nu)
:class-container: derivation

At the wall $w = 0$, so the flux is purely conductive and $|\partial_z\overline T| = Nu$
there ({doc}`nusselt`). The tangent at the wall falls from the plate temperature by 1/2,
to the bulk value, over a distance $\delta_T$ with $Nu\,\delta_T = 1/2$.
```

The slope definitions are robust; the peak definitions are closer to what one sees in a
profile but the peak can be broad. The ratio figure uses the slope definitions for both
layers, as Prandtl–Blasius theory does.

## Prandtl–Blasius ratios

For a laminar boundary layer over a plate, [Shishkina et al.
(2010)](https://arxiv.org/abs/1109.6870) give the ratio of the layers' thicknesses as

$$
\frac{\delta_T}{\delta_u} \approx
\begin{cases}
0.588\,Pr^{-1/2}, & Pr < 3\times10^{-4},\\
Pr^{-0.357 + 0.022\log Pr}, & 3\times10^{-4} \le Pr \le 3 \text{ (a fit)},\\
0.982\,Pr^{-1/3}, & Pr > 3,
\end{cases}
$$

the two power laws crossing at $Pr^* \approx 0.046$. {mod}`rbc.plots.boundary_layers`
draws the two asymptotes for reference. Note that most runs here lie in the transition
range, where neither asymptote applies exactly; the fit for that range (and the base of its
logarithm) should be checked against the paper before it is drawn.

(chebyshev-nodes)=
## Nodes inside a boundary layer

Dedalus' Chebyshev grid on $[0, 1]$ with $N_z$ modes is

$$
z_j = \tfrac12\Big(1 - \cos\frac{\pi(j + \frac12)}{N_z}\Big),\qquad j = 0,\dots,N_z - 1,
$$

so the number of nodes within a distance $\delta$ of a plate is

$$
n(\delta) \approx \frac{N_z}{\pi}\arccos(1 - 2\delta) \approx \frac{2N_z}{\pi}\sqrt{\delta}
\quad(\delta \ll 1).
$$

```{dropdown} Derivation
:class-container: derivation

$z_j < \delta$ exactly when $\cos(\pi(j+\frac12)/N_z) > 1 - 2\delta$, that is
$j + \frac12 < (N_z/\pi)\arccos(1 - 2\delta)$, which counts the nodes. For small
$\delta$, $\arccos(1 - 2\delta) \approx 2\sqrt\delta$, since $\cos x \approx 1 - x^2/2$.
```

Because of the square root, halving a layer's thickness costs only $\sqrt2$ more nodes
for the same coverage: the Chebyshev grid's clustering is what makes thin layers
affordable. The `n_delta_*` columns of `results.csv` are these counts for the simulation
grid, and the `grid` figure of {mod}`rbc.plots.boundary_layers` draws every node against the
profiles.

For orientation, Shishkina et al. recommend lower bounds on the number of nodes inside each
layer; for $Pr \approx 0.7$ their fits are $N_{\text{th}} \gtrsim 0.35\,Ra^{0.15}$ and
$N_{\text{v}} \gtrsim 0.31\,Ra^{0.15}$, with other prefactors at other Pr. Those were
derived for finite-volume grids; with spectral accuracy the decisive test here is the
resolution study ({doc}`resolution`).
