"""Simulation helpers, analysis and figures for 3D Rayleigh–Bénard convection DNS.

Run from the command line as ``python -m rbc.<module>``; see the documentation for the
workflow. Subpackages:

``rbc.analysis``
    per-run analysis (Nusselt and Reynolds numbers, profiles) and scaling exponents.
``rbc.spectra``
    spectral resolution diagnostics.
``rbc.plots``
    paper figures, built from ``data/cases.csv`` and ``data/results.csv``.
``rbc.vis``
    field visualizations from snapshots.
"""

__version__ = "0.2.0"
