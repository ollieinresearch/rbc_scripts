"""
Contours over the parameter plane
=================================

Filled contours of :math:`\\log_{10} Nu` over (Ra, Pr), interpolated on a triangulation of
the runs (white dots) and drawn only inside their convex hull, with the slices overlaid.

Commands::

    python -m rbc.plots.contour
    python -m rbc.plots.contour --comp 0.25
"""
import tempfile

from rbc import style
from rbc.plots import contour

style.KEEP_OPEN = True
out = tempfile.mkdtemp()
contour.main(["--no-tex", "--out", out])

# %%
# Compensated by :math:`Ra^{1/4}`: along a slice the Shishkina–Lohse prediction would
# follow a contour line.
contour.main(["--comp", "0.25", "--no-tex", "--out", out])
