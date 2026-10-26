"""
Local exponents α and β
=======================

Local exponents of :math:`Nu - 1 \\sim Pr^{\\alpha}Ra^{\\beta}`: two-point estimates at the
midpoints between neighbouring runs (circles, joined along each row of constant Pr or Ra)
and both exponents at cell centres (diamonds), against Ra and against Pr, coloured by the
other parameter. The black line is a weighted fit of the exponent against the parameter;
its slope and the mean are in each panel's title. The printed tables compare surface fits
and list the slope along the slices at the cell centres. ``--draw-surface`` overlays the
surface's exponents, and ``--min-q 1`` leaves out estimates that use a run with
:math:`Nu - 1 < 1`. Method: :doc:`/theory/exponents`.

Command::

    python -m rbc.plots.exponents --surface 2
"""
import tempfile

from rbc import style
from rbc.plots import exponents

style.KEEP_OPEN = True
exponents.main(["--surface", "2", "--no-tex", "--out", tempfile.mkdtemp()])
