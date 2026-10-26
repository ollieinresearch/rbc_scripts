"""
Nu − 1 at fixed Ra and at fixed Pr
==================================

Top: :math:`Nu - 1` against Pr, one curve per Ra. Bottom: against Ra, one curve per Pr.
With ``--slopes`` each legend entry gives the curve's fitted exponent.

The second figure is compensated by the reference law fitted to all points, with a
quadratic fit through each curve, so drifting exponents show as curvature; the third
divides by the fixed parameter, so the curves collapse onto the reference law if
:math:`Nu - 1` separates into a function of Pr times a function of Ra.

Commands::

    python -m rbc.plots.fixed --slopes
    python -m rbc.plots.fixed --slopes --norm comp --fit quadratic
    python -m rbc.plots.fixed --slopes --norm collapse --fit linear
"""
import tempfile

from rbc import style
from rbc.plots import fixed

style.KEEP_OPEN = True
out = tempfile.mkdtemp()
fixed.main(["--slopes", "--no-tex", "--out", out])

# %%
# Compensated, with a quadratic fit per curve
fixed.main(["--slopes", "--norm", "comp", "--fit", "quadratic", "--no-tex", "--out", out])

# %%
# Collapsed onto the reference law
fixed.main(["--slopes", "--norm", "collapse", "--fit", "linear", "--no-tex", "--out", out])
