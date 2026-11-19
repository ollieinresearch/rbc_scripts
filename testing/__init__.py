"""Resolution testing kept apart from the main ``rbc`` package while it is being validated.

Run from the repository root (``python -m testing.<module>``); see ``testing/README.md``.

* :mod:`testing.convergence`: the runs of one (Ra, Pr) at several resolutions, the plateau
  they converge to and each run's distance from it.
* :mod:`testing.calibrate`: every multi-resolution point at once, and which resolution
  measure (k_max η, decades, the vorticity identity, ...) separates converged runs from the
  rest with one threshold.
* :mod:`testing.check`: other runs against those thresholds.
"""
