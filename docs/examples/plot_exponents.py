"""
Local exponents α and β
=======================

Local exponents of :math:`Nu - 1 \\sim Pr^{\\alpha}Ra^{\\beta}`: two-point estimates at the
midpoints between neighbouring runs (circles) and both exponents at cell centres
(diamonds), against Ra and against Pr, coloured by the other parameter. Each panel has a
weighted constant (dashed) and line (solid) fit, to test whether the exponent changes.
Thin curves are the exponents of the quadratic surface fit. The printed tables compare
surface fits and list the slope along the slices at the cell centres. Method:
:doc:`/theory/exponents`.

Command::

    python -m rbc.plots.exponents --surface 2
"""
import tempfile

from rbc import style
from rbc.plots import exponents

style.KEEP_OPEN = True
exponents.main(["--surface", "2", "--no-tex", "--out", tempfile.mkdtemp()])
