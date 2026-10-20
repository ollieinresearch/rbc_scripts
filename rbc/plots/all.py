"""Regenerate every paper figure from ``data/cases.csv`` and ``data/results.csv``.

Usage::

    python -m rbc.plots.all [--out figures] [--png] [--no-tex]
    python -m rbc.plots.all --results nusselt.csv     # a simple Ra,Pr,Nusselt file

With a simple three-column ``Ra,Pr,Nusselt`` file (see
:func:`rbc.cases.read_simple_results`) every figure that needs only Nu is made, with the
parameter map drawn from the file's points; the Re, boundary-layer and resolution
figures need the full results table and are skipped.

A figure that fails (for example for lack of data) is reported and skipped; the others
are still made.
"""
from __future__ import annotations

import argparse
import traceback

from rbc import cases as rc
from rbc.plots import boundary_layers, contour, exponents, fixed, param_grid, resolution, scaling


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.plots.all", description=__doc__.split("\n")[0])
    p.add_argument("--out", default="figures", help="output folder")
    p.add_argument("--results", help="results table: results.csv, or a simple Ra,Pr,Nusselt CSV "
                                     "(default: data/results.csv)")
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
    simple = bool(a.results) and rc.is_simple_results(a.results)
    res = ["--results", a.results] if a.results else []
    if simple:
        jobs = [
            (param_grid.main, ["--dim", "3", "--results", a.results]),
            (scaling.main, ["--q", "Nu"] + res),
            (scaling.main, ["--q", "Nu", "--norm", "law"] + res),
            (scaling.main, ["--q", "Nu-1"] + res),
            (fixed.main, ["--slopes"] + res),
            (exponents.main, ["--surface", "2"] + res),
            (contour.main, [] + res),
        ]
        print("simple Ra,Pr,Nusselt file: skipping the Re, boundary-layer and resolution figures")
    else:
        jobs = [
            (param_grid.main, ["--dim", "3"]),
            (scaling.main, ["--q", "Nu"] + res),
            (scaling.main, ["--q", "Nu", "--norm", "law"] + res),
            (scaling.main, ["--q", "Re"] + res),
            (fixed.main, ["--slopes"] + res),
            (exponents.main, ["--surface", "2"] + res),
            (boundary_layers.main, ["ratio"] + res),
            (resolution.main, ["study"] + res),
            (contour.main, [] + res),
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