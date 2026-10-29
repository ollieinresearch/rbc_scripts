"""
Nu − 1 against Ra Pr^a along slices, and its local exponent
===========================================================

Left: slices :math:`Pr \\propto Ra^{-1/2}` against :math:`x = Ra\\,Pr^{-1/2}`. Right: slices
:math:`Pr \\propto Ra^{-1}` against :math:`x = Ra/Pr = Re_{\\rm ff}^2`. Bottom: the local
exponent :math:`\\gamma = d\\log(Nu-1)/d\\log x` of each slice, with the value the
Shishkina–Lohse low-Pr ultimate law gives and the one a single power law fitted to all
points gives. Along a slice :math:`\\gamma` is the slice exponent divided by
:math:`1 + a^2`; what the shared axis adds is whether the slices' :math:`\\gamma(x)` fall on
one curve. See :ref:`reduced-exponent`.

Command::

    python -m rbc.plots.reduced
"""
import tempfile

from rbc import style
from rbc.plots import reduced

style.KEEP_OPEN = True
reduced.main(["--no-tex", "--out", tempfile.mkdtemp()])
