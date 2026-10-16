"""
Parameter map
=============

Every point of ``cases.csv``: filled circles are finished runs, open circles running ones,
stars the points with a resolution study, and coloured lines the slices
:math:`Pr \\propto Ra^{-1/2}`, labelled by their Prandtl number at :math:`Ra = 10^5`.

Command::

    python -m rbc.plots.param_grid --dim 3 --legend
"""
import tempfile

from rbc import style
from rbc.plots import param_grid

style.KEEP_OPEN = True
param_grid.main(["--dim", "3", "--legend", "--no-tex", "--out", tempfile.mkdtemp()])

# %%
# With ``--ra-inv`` the slices :math:`Pr \propto Ra^{-1}` are added, dashed.
param_grid.main(["--dim", "3", "--ra-inv", "--no-tex", "--out", tempfile.mkdtemp()])
