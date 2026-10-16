# Pictures and movies

Visualizations of the simulations, made with {doc}`../user-guide/visualization`.

```{admonition} To fill in
:class: todo-note

Add pictures and movies here as they are chosen.
```

## Adding one

Put the file in `docs/_static/media/` (keep movies small: a few MB, H.264 `.mp4`), then add
a card to the grid below. A picture:

````markdown
:::{grid-item-card} Ra = 10^8, Pr = 10^-2.5: temperature isosurfaces
:img-top: ../_static/media/ra8_pr-2.5_iso.png

Warm (red) and cool (blue) plumes at $\theta = \pm1.5\,\theta_\text{rms}$.
`python -m rbc.vis.render3d RUN --mode iso`
:::
````

A movie:

````markdown
:::{grid-item-card} Ra = 10^7, Pr = 1 and 10^-2: midplane slices
```{raw} html
<video controls loop muted width="100%" src="../_static/media/slices_ra7.mp4"></video>
```
:::
````

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} Example card
Replace with the first picture.
:::
::::
