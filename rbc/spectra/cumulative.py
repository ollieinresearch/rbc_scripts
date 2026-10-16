"""Time-averaged power spectra and "decades spanned" diagnostics for one run.

Reads ``RUN/spectra/spectra_s*.h5`` (in-simulation, preferred) or
``RUN/snapshots/snapshots_s*_spectra.h5`` (from :mod:`rbc.spectra.snapshots`) and
writes to ``RUN/outputs/``: ``cumulative_spectra_velocity.pdf``,
``cumulative_spectra_temperature.pdf``, ``decades_timeseries.pdf``, ``spec_range.txt``
and ``spectra.json`` (read by :mod:`rbc.collect`). Usage::

    python -m rbc.spectra.cumulative RUN_DIR [--source auto|sim|snapshots] [--time T0] [--end_time T1]

Two summaries of "decades spanned" are reported side by side:

``D[<E>]``
    decades of the time-averaged spectrum (the primary number).
``<D(t)>``
    time average of the decades of each instantaneous spectrum, with its standard
    deviation, minimum and 5th percentile.

They differ because averaging and taking logs do not commute. The tail of <E> is an
arithmetic mean, dominated by the bursts when the grid-scale energy is largest, so
D[<E>] <= typical D(t). The truncation error that biases a time-averaged Nusselt
number is set by that same arithmetic mean, so D[<E>] is the conservative number to
quote. The D(t) distribution shows whether a shortfall is steady or intermittent: a
low minimum or 5th percentile with a healthy mean means the grid is adequate most of
the time but under-resolves the strongest events (plume impacts, boundary-layer
eruptions).

Decades are log10(peak) minus the mean of log10 E over the last 3 modes, not max - min;
see :func:`rbc.spectra.core.decades` for why this matters when energy piles up at the
grid scale.
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from rbc import style
from rbc.spectra import core as cs

style.use("diag")

#: Panels of the spectra figures (see :data:`rbc.spectra.core.PANEL_META`).
PANELS = ("3d", "nz", "h", "perp")
FIELD_NAMES = {"vel": "Velocity", "temp": "Temperature fluctuation"}


def select_window(time, data, t0, t1):
    """Restrict spectra to a time window.

    Parameters
    ----------
    time : numpy.ndarray
        Write times.
    data : dict
        ``data[field][key]`` arrays with time first.
    t0, t1 : float
        Window; ``t1 = 0`` means the end of the data.

    Returns
    -------
    tuple
        ``(time, data)`` restricted to the window.

    Raises
    ------
    ValueError
        With fewer than two spectra in the window.
    """
    keep = time >= t0
    if t1 > 0:
        keep &= time <= t1
    if keep.sum() < 2:
        raise ValueError(f"Fewer than 2 spectra in the window [{t0}, {t1 or 'end'}].")
    return time[keep], {f: {k: v[keep] for k, v in d.items()} for f, d in data.items()}


def summarize(time, modes, data, attrs):
    """Time-mean spectra and decades for every field and panel.

    Parameters
    ----------
    time, modes, data, attrs
        As returned by :func:`rbc.spectra.core.load_spectra` (windowed).

    Returns
    -------
    dict
        ``out[field][panel]`` with ``n``, ``E_mean``, ``E_lo``/``E_hi`` (5th/95th
        percentiles over time), ``d_mean`` (decades of the mean spectrum) and ``d_t``
        (decades of each write); ``out[field]["mean_profile"]`` for the Chebyshev
        spectrum of the horizontal mean.
    """
    out = {}
    for fld in cs.FIELDS:
        series = cs.panel_series(modes, data[fld], attrs)
        res = {}
        for key in PANELS:
            n, E_t = series[key]                    # E_t: (nt, n)
            E_mean = cs.time_average(time, E_t)
            d_mean = cs.decades(n, E_mean)
            d_t = cs.decades(n, E_t)                # arrays over time
            res[key] = dict(n=n, E_mean=E_mean,
                            E_lo=np.percentile(E_t, 5, axis=0),
                            E_hi=np.percentile(E_t, 95, axis=0),
                            d_mean=d_mean, d_t=d_t)
        mean_prof = cs.time_average(time, data[fld]["mean_E_nz"])
        res["mean_profile"] = dict(n=modes["nz"], E_mean=mean_prof,
                                   d_mean=cs.decades(modes["nz"], mean_prof))
        out[fld] = res
    return out


def plot_field(fld, res, t_range, savepath):
    """Four-panel figure of one field's time-averaged spectra.

    Parameters
    ----------
    fld : {"vel", "temp"}
        Field.
    res : dict
        ``summarize(...)[fld]``.
    t_range : str
        Subtitle describing the window.
    savepath : pathlib.Path
        File to write.
    """
    fig, axes = plt.subplots(2, 2, figsize=(16, 14), layout="constrained")
    for key in PANELS:
        (r, c), ptitle, xlabel = cs.PANEL_META[key]
        ax = axes[r, c]
        p = res[key]
        n = p["n"]
        m = (n >= 1) & (p["E_mean"] > 0)
        ax.fill_between(n[m], p["E_lo"][m], p["E_hi"][m], alpha=0.25, linewidth=0,
                        label="5-95% of snapshots")
        ax.loglog(n[m], p["E_mean"][m], ".-", linewidth=1.5, markersize=4, label="time average")
        extra = [p["E_mean"][m], p["E_hi"][m]]
        if key == "nz":
            q = res["mean_profile"]
            mm = (q["n"] >= 1) & (q["E_mean"] > 0)
            ax.loglog(q["n"][mm], q["E_mean"][mm], "--", linewidth=1.2, label="horizontal mean")
            extra.append(q["E_mean"][mm])
        dm, dt = p["d_mean"], p["d_t"]
        txt = (f"D[<E>] = {dm['D']:.2f}\n"
               f"<D(t)> = {np.nanmean(dt['D']):.2f} ± {np.nanstd(dt['D']):.2f}\n"
               f"min D(t) = {np.nanmin(dt['D']):.2f}")
        if dm["pileup"] > 0.1:
            txt += f"\npile-up {dm['pileup']:.2f} dec"
        ax.text(0.97, 0.97, txt, transform=ax.transAxes, ha="right", va="top", fontsize=11,
                bbox=dict(boxstyle="round", fc="white", alpha=0.8))
        ax.set_title(ptitle)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("E")
        cs.setup_yaxis(ax, *cs.auto_ylim(*extra))
        ax.legend(fontsize=10, loc="lower left")
    fig.suptitle(f"Time-averaged {FIELD_NAMES[fld]} spectra\n{t_range}")
    fig.savefig(savepath, transparent=True, dpi=300)
    plt.close(fig)


def plot_decades_timeseries(time, summ, savepath):
    """Decades of each write against time, for every field and panel.

    Parameters
    ----------
    time : numpy.ndarray
        Write times.
    summ : dict
        Output of :func:`summarize`.
    savepath : pathlib.Path
        File to write.
    """
    fig, axes = plt.subplots(len(PANELS), 2, figsize=(16, 3.2 * len(PANELS)),
                             layout="constrained", sharex=True)
    for j, fld in enumerate(cs.FIELDS):
        for i, key in enumerate(PANELS):
            ax = axes[i, j]
            p = summ[fld][key]
            ax.plot(time, p["d_t"]["D"], linewidth=0.8, label="D(t)")
            ax.axhline(p["d_mean"]["D"], color="C1", label="D[<E>]")
            ax.axhline(np.nanmean(p["d_t"]["D"]), color="C2", linestyle="--", label="<D(t)>")
            ax.set_ylabel("decades")
            ax.set_title(f"{FIELD_NAMES[fld]} - {cs.REPORT_NAMES[key]}", fontsize=12)
            if i == 0 and j == 0:
                ax.legend(fontsize=10)
        axes[-1, j].set_xlabel("t")
    fig.savefig(savepath, dpi=200)
    plt.close(fig)


def write_report(path, time, summ, attrs, files):
    """Write the text table ``spec_range.txt``.

    Parameters
    ----------
    path : pathlib.Path
        File to write.
    time : numpy.ndarray
        Write times.
    summ : dict
        Output of :func:`summarize`.
    attrs : dict
        Spectra file attributes.
    files : list of pathlib.Path
        Files read (for the header).
    """
    with open(path, "w") as r:
        r.write(f"Spectra source: {attrs.get('source', '?')}  ({len(files)} file(s))\n")
        r.write(f"Resolution: Nx={attrs.get('Nx')}, Ny={attrs.get('Ny')}, Nz={attrs.get('Nz')}\n")
        r.write(f"Averaging window: t in [{time[0]:.4f}, {time[-1]:.4f}], {len(time)} spectra\n")
        r.write("D = log10(peak) - mean(log10 E over last 3 modes); peak over n >= 1.\n")
        r.write("=" * 96 + "\n")
        hdr = (f"{'field':<12}{'panel':<34}{'D[<E>]':>8}{'<D(t)>':>8}{'std':>7}"
               f"{'min':>7}{'p05':>7}{'maxmin':>8}{'pileup':>8}\n")
        r.write(hdr)
        r.write("-" * 96 + "\n")
        for fld in cs.FIELDS:
            for key in PANELS:
                p = summ[fld][key]
                D_t = p["d_t"]["D"]
                r.write(f"{fld:<12}{cs.REPORT_NAMES[key]:<34}"
                        f"{p['d_mean']['D']:8.3f}{np.nanmean(D_t):8.3f}{np.nanstd(D_t):7.3f}"
                        f"{np.nanmin(D_t):7.3f}{np.nanpercentile(D_t, 5):7.3f}"
                        f"{p['d_mean']['D_maxmin']:8.3f}{p['d_mean']['pileup']:8.3f}\n")
            q = summ[fld]["mean_profile"]["d_mean"]
            r.write(f"{fld:<12}{'horizontal-mean profile (z)':<34}{q['D']:8.3f}"
                    f"{'':>8}{'':>7}{'':>7}{'':>7}{q['D_maxmin']:8.3f}{q['pileup']:8.3f}\n")
        r.write("-" * 96 + "\n")
        r.write("maxmin = old max-min definition; pileup = tail - smallest 3-mode running mean "
                "of log10 E; > 0.1 means energy accumulating at the grid scale.\n")


def write_json(path, time, summ, attrs):
    """Write ``spectra.json`` for :mod:`rbc.collect`.

    Parameters
    ----------
    path : pathlib.Path
        File to write.
    time : numpy.ndarray
        Write times.
    summ : dict
        Output of :func:`summarize`.
    attrs : dict
        Spectra file attributes (Nx, Ny, Nz).

    Notes
    -----
    Keys: ``D_<field>_<panel>`` (decades of the mean spectrum), ``Dmin_<field>_<panel>``
    (smallest decades of any write), ``pileup_<field>_<panel>`` and
    ``D_<field>_meanz`` (horizontal-mean profile), with field in {vel, temp} and panel
    in :data:`PANELS`.
    """
    import json
    rec = dict(t0=float(time[0]), t1=float(time[-1]), n=int(len(time)),
               Nx=int(attrs.get("Nx", 0)), Ny=int(attrs.get("Ny", 0)), Nz=int(attrs.get("Nz", 0)))
    for fld in cs.FIELDS:
        for key in PANELS:
            p = summ[fld][key]
            rec[f"D_{fld}_{key}"] = float(p["d_mean"]["D"])
            rec[f"Dmin_{fld}_{key}"] = float(np.nanmin(p["d_t"]["D"]))
            rec[f"pileup_{fld}_{key}"] = float(p["d_mean"]["pileup"])
        rec[f"D_{fld}_meanz"] = float(summ[fld]["mean_profile"]["d_mean"]["D"])
    path.write_text(json.dumps(rec, indent=1))


def run(basepath, source="auto", t0=0.0, t1=0.0):
    """Summarize, plot and report the spectra of one run.

    Parameters
    ----------
    basepath : str or pathlib.Path
        Run directory.
    source : {"auto", "sim", "snapshots"}, optional
        Which spectra files to read (:func:`rbc.spectra.core.spectra_files`).
    t0, t1 : float, optional
        Averaging window (``t1 = 0``: end of the data).

    Returns
    -------
    dict
        Output of :func:`summarize`.
    """
    basepath = Path(basepath)
    out = basepath / "outputs"
    out.mkdir(exist_ok=True)
    files = cs.spectra_files(basepath, source)
    time, modes, data, attrs = cs.load_spectra(files)
    time, data = select_window(time, data, t0, t1)
    print(f"{len(files)} file(s), {len(time)} spectra, t in [{time[0]:.4f}, {time[-1]:.4f}]")

    summ = summarize(time, modes, data, attrs)
    t_range = f"t in [{time[0]:.3f}, {time[-1]:.3f}], {len(time)} spectra"
    plot_field("vel", summ["vel"], t_range, out / "cumulative_spectra_velocity.pdf")
    plot_field("temp", summ["temp"], t_range, out / "cumulative_spectra_temperature.pdf")
    plot_decades_timeseries(time, summ, out / "decades_timeseries.pdf")
    write_report(out / "spec_range.txt", time, summ, attrs, files)
    write_json(out / "spectra.json", time, summ, attrs)
    print(f"Wrote figures and spec_range.txt to {out}")
    return summ


def get_parser() -> argparse.ArgumentParser:
    """Command-line options.

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.spectra.cumulative",
                                description="Time-averaged spectra and decades for one run.")
    p.add_argument("basepath", type=Path, help="run directory")
    p.add_argument("--source", choices=("auto", "sim", "snapshots"), default="auto",
                   help="in-simulation spectra, snapshot-derived spectra, or whichever exists")
    p.add_argument("--time", "--t0", type=float, default=0.0, dest="t0",
                   help="start of the averaging window (exclude the transient)")
    p.add_argument("--end_time", "--t1", type=float, default=0.0, dest="t1",
                   help="end of the averaging window (0 = end of the data)")
    return p


def main(argv=None) -> None:
    """Command-line entry point for :func:`run`.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments (default: ``sys.argv[1:]``).
    """
    a = get_parser().parse_args(argv)
    run(a.basepath, a.source, a.t0, a.t1)


if __name__ == "__main__":
    main()
