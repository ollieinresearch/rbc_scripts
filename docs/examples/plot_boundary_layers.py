"""
Boundary layers
===============

Left: the ratio of viscous to thermal boundary-layer thickness against Pr, with the
Prandtl–Blasius asymptotes. Right: Chebyshev nodes inside each layer against Ra. See
:doc:`/theory/boundary-layers`.

Command::

    python -m rbc.plots.boundary_layers ratio
"""
import tempfile

from rbc import style
from rbc.plots import boundary_layers

style.KEEP_OPEN = True
boundary_layers.main(["ratio", "--no-tex", "--out", tempfile.mkdtemp()])
