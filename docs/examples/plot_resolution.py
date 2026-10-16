"""
Resolution study
================

For each point with several resolutions: Nu from the three identities relative to the
finest run (top) and the spectral decades of the velocity and temperature spectra (bottom),
against :math:`N_z`. Grey lines mark the production run. See :doc:`/theory/resolution`.

Command::

    python -m rbc.plots.resolution study
"""
import tempfile

from rbc import style
from rbc.plots import resolution

style.KEEP_OPEN = True
resolution.main(["study", "--no-tex", "--out", tempfile.mkdtemp()])
