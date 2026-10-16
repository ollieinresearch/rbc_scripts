"""Overlay time-averaged spectra of every resolution under one (Ra, Pr) folder.
Curves stack at low k and peel off at each run's cutoff; divergence at low k
means those modes are genuinely under-resolved. Usage::

    python -m rbc.spectra.compare PARAM_DIR [--panels nz h] [--t0 T]

Each run's averaging start is taken from its ``outputs/nusselt.json`` unless ``--t0`` is
given. Writes ``PARAM_DIR/outputs/spectra_compare.pdf``.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


from rbc import style
from rbc.cases import parse_run_path
from rbc.spectra import core as cs

PANEL_LABEL = {"nz": r"Chebyshev degree $n_z$", "h": r"horizontal mode $k$",
               "perp": r"$k_\perp$", "3d": r"3D index shell $s$"}
FIELD_LABEL = {"vel": "velocity", "temp": r"temperature fluctuation $b$"}


def mean_spectra(run: Path, t0: float | None = None) -> dict | None:
    """Time-averaged spectra of one run.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.
    t0 : float, optional
        Averaging start (default: the run's ``nusselt.json`` window start, else 0).

    Returns
    -------
    dict or None
        ``attrs``, ``t0``, ``t1`` and, per ``(field, panel)``, ``(n, E_mean, D)``; None
        if the run has no spectra or fewer than two writes after ``t0``.
    """
    files = cs.spectra_files(run)
    if not files:
        return None
    if t0 is None:
        nj = run / "outputs" / "nusselt.json"
        t0 = json.loads(nj.read_text())["t0"] if nj.exists() else 0.0
    time, modes, data, attrs = cs.load_spectra(files)
    keep = time >= t0
    if keep.sum() < 2:
        return None
    time = time[keep]
    out = dict(attrs=attrs, t0=float(time[0]), t1=float(time[-1]))
    for fld in cs.FIELDS:
        series = cs.panel_series(modes, {k: v[keep] for k, v in data[fld].items()}, attrs)
        for key, (n, E) in series.items():
            Em = cs.time_average(time, E)
            out[(fld, key)] = (n, Em, float(cs.decades(n, Em)["D"]))
    return out


def runs_in(param_dir: Path) -> list[Path]:
    """Run directories with spectra under one (Ra, Pr) folder, by resolution.

    Parameters
    ----------
    param_dir : pathlib.Path
        e.g. ``.../ra1e8/pr1e0``.

    Returns
    -------
    list of pathlib.Path
    """
    runs = [p for p in param_dir.iterdir() if p.is_dir() and parse_run_path(p) and (p / "spectra").exists()]
    return sorted(runs, key=lambda p: (parse_run_path(p)["nz"] or 0, p.name))


def plot_overlay(specs: dict[str, dict], panels, path: Path, title=""):
    """Overlay several runs' mean spectra, one panel per (field, panel).

    Parameters
    ----------
    specs : dict
        Run label -> output of :func:`mean_spectra`.
    panels : list of str
        Panels to draw (keys of ``PANEL_LABEL``).
    path : pathlib.Path
        File to write.
    title : str, optional
        Figure title.
    """
    import matplotlib.pyplot as plt
    from matplotlib import colormaps

    labels = list(specs)
    cmap = colormaps["viridis"]
    colors = {lab: cmap(i / max(len(labels) - 1, 1) * 0.9) for i, lab in enumerate(labels)}
    fig, axes = plt.subplots(len(cs.FIELDS), len(panels), figsize=(5.2 * len(panels), 4.2 * len(cs.FIELDS)),
                             squeeze=False)
    for i, fld in enumerate(cs.FIELDS):
        for j, key in enumerate(panels):
            ax = axes[i, j]
            arrays = []
            for lab in labels:
                n, E, D = specs[lab][(fld, key)]
                m = (n >= 1) & (E > 0)
                ax.loglog(n[m], E[m], color=colors[lab], label=f"{lab}  ($D={D:.1f}$)")
                arrays.append(E[m])
            cs.setup_yaxis(ax, *cs.auto_ylim(*arrays))
            ax.set_xlabel(PANEL_LABEL[key])
            ax.set_ylabel(f"$E$, {FIELD_LABEL[fld]}")
            ax.legend(fontsize=8, loc="lower left")
    fig.suptitle(title)
    fig.savefig(path)
    plt.close(fig)


def get_parser() -> argparse.ArgumentParser:
    """Command-line options.

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.spectra.compare", description=__doc__.split("\n")[0])
    p.add_argument("param_dir", type=Path, help="folder holding every resolution of one (Ra, Pr)")
    p.add_argument("--panels", nargs="+", default=["nz", "h"], choices=list(PANEL_LABEL),
                   help="spectra to show: nz (Chebyshev), h (horizontal), perp, 3d")
    p.add_argument("--t0", type=float, help="averaging start for every run")
    return p


def main(argv=None) -> None:
    """Overlay the spectra of every resolution under one parameter folder.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments (default: ``sys.argv[1:]``).
    """
    a = get_parser().parse_args(argv)
    style.use("diag")
    specs = {}
    for run in runs_in(a.param_dir):
        s = mean_spectra(run, a.t0)
        if s is not None:
            specs[run.name] = s
    if not specs:
        raise SystemExit(f"no spectra under {a.param_dir}")
    out = a.param_dir / "outputs"
    out.mkdir(exist_ok=True)
    plot_overlay(specs, a.panels, out / "spectra_compare.pdf", a.param_dir.as_posix().split("3d/")[-1])
    print(out / "spectra_compare.pdf")


if __name__ == "__main__":
    main()
