"""
Nu − 1 at fixed Ra and at fixed Pr
==================================

Left: :math:`Nu - 1` against Pr, one curve per Ra. Right: against Ra, one curve per Pr.
With ``--slopes`` each legend entry gives the curve's fitted exponent.

Command::

    python -m rbc.plots.fixed --slopes
"""
import tempfile

from rbc import style
from rbc.plots import fixed

style.KEEP_OPEN = True
fixed.main(["--slopes", "--no-tex", "--out", tempfile.mkdtemp()])
