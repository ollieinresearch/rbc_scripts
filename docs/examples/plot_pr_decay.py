"""
How Nu − 1 decays with Pr
=========================

Top: :math:`(Nu-1)/(Nu-1)_{Pr=1}` against Pr, one curve per Ra, so every curve ends at
(1, 1) and only the shape of the Pr dependence is compared. Bottom: the local exponent
:math:`d\\log(Nu-1)/d\\log Pr` from neighbouring Pr, with the Pr exponents of the low-Pr
Grossmann–Lohse regimes as guides. A curve whose local exponent climbs past all of them as
Pr falls is leaving those regimes towards weak convection. The synthetic data here follow
one power law, so the exponent is flat.

Command::

    python -m rbc.plots.pr_decay
"""
import tempfile

from rbc import style
from rbc.plots import pr_decay

style.KEEP_OPEN = True
pr_decay.main(["--no-tex", "--out", tempfile.mkdtemp()])
