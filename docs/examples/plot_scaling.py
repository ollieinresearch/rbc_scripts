"""
Compensated Nu along slices
===========================

Nu divided by :math:`Ra^{1/4}`, the Shishkina–Lohse prediction along
:math:`Pr \\propto Ra^{-1/2}`, one curve per slice, each with its fitted slope :math:`m` in
the legend. The four normalizations place the slices differently but give identical
slopes; see :doc:`/theory/normalization`.

Commands::

    python -m rbc.plots.scaling --q Nu --norm none
    python -m rbc.plots.scaling --q Nu --norm law
    python -m rbc.plots.scaling --q Nu --norm common
    python -m rbc.plots.scaling --q Nu --norm ref --ref-ra 7
"""
import tempfile

from rbc import style
from rbc.plots import scaling

style.KEEP_OPEN = True
out = tempfile.mkdtemp()

# %%
# No normalization: offsets between slices carry the Pr dependence.
scaling.main(["--q", "Nu", "--norm", "none", "--no-tex", "--out", out])

# %%
# ``law``: :math:`Nu/(Pr^{1/2}Ra^{1/2})` for every point, no free constant. Collapse onto
# one horizontal line would mean the reference law holds.
scaling.main(["--q", "Nu", "--norm", "law", "--no-tex", "--out", out])

# %%
# ``common``: each slice divided by its geometric mean over the Ra values all slices share
# (shaded).
scaling.main(["--q", "Nu", "--norm", "common", "--no-tex", "--out", out])

# %%
# ``ref``: every slice equal to 1 at :math:`Ra = 10^7`.
scaling.main(["--q", "Nu", "--norm", "ref", "--ref-ra", "7", "--no-tex", "--out", out])

# %%
# The Reynolds number, compensated by :math:`Ra^{3/4}`, its predicted scaling along the
# same slices.
scaling.main(["--q", "Re", "--no-tex", "--out", out])
