"""Regenerate every paper figure from ``data/cases.csv`` and ``data/results.csv``.

Usage::

    python -m rbc.plots.all [--out figures] [--png] [--no-tex]

A figure that fails (for example for lack of data) is reported and skipped; the others
are still made.
"""
from __future__ import annotations

import argparse
import traceback

from rbc.plots import boundary_layers, contour, exponents, fixed, param_grid, resolution, scaling


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.plots.all", description=__doc__.split("\n")[0])
    p.add_argument("--out", default="figures", help="output folder")
    p.add_argument("--png", action="store_true", help="also write PNG previews")
    p.add_argument("--no-tex", action="store_true", help="never use LaTeX for text (default: use it where it is installed)")
    return p


def main(argv=None) -> None:
    """Command-line entry point.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments (default: ``sys.argv[1:]``).
    """
    a = get_parser().parse_args(argv)
    common = ["--out", a.out] + (["--png"] if a.png else []) + (["--no-tex"] if a.no_tex else [])
    jobs = [
        (param_grid.main, ["--dim", "3"]),
        (scaling.main, ["--q", "Nu"]),
        (scaling.main, ["--q", "Nu", "--norm", "law"]),
        (scaling.main, ["--q", "Re"]),
        (fixed.main, ["--slopes"]),
        (exponents.main, ["--surface", "2"]),
        (boundary_layers.main, ["ratio"]),
        (resolution.main, ["study"]),
        (contour.main, []),
    ]
    for fn, args in jobs:
        try:
            fn(args + common)
        except SystemExit as e:
            print(f"{fn.__module__}: {e}")
        except Exception:
            print(f"{fn.__module__} failed:\n{traceback.format_exc(limit=2)}")


if __name__ == "__main__":
    main()
