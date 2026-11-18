"""Fits and small helpers shared by the paper figures.

Error bars (``--errorbars`` on every figure command)
    ``all`` (default): a capped bar of ±1 standard error (statistical) and a thin uncapped
    whisker out to ±(standard error + systematic bound), the bound ``Nu_sys`` from
    resolution and timestepping (:mod:`rbc.analysis.variants`); fitted slopes are quoted
    as ``m ± stat ± sys``. ``stat``: the statistical bar only. ``none``: no bars and no ±
    in the legends. The systematic part is drawn only where it is known.
"""
from __future__ import annotations

import numpy as np

LN10 = np.log(10)
#: Choices of ``--errorbars``.
ERRORBARS = ("all", "stat", "none")
_BAR_KEYS = ("ecolor", "elinewidth", "capsize", "capthick")


def add_errorbar_arg(p) -> None:
    """Add ``--errorbars`` to a figure's argument parser.

    Parameters
    ----------
    p : argparse.ArgumentParser
    """
    p.add_argument("--errorbars", choices=ERRORBARS, default="all",
                   help="all: statistical bars and statistical + systematic whiskers; stat: statistical only; "
                        "none: no error bars")


def sys_values(df, q):
    """Systematic bounds of the plotted quantity (absolute, NaN where unknown).

    Parameters
    ----------
    df : pandas.DataFrame
        Results rows.
    q : {"Nu", "Nu-1", "Re"}
        Quantity. Only Nu has a bound (the same for Nu and Nu - 1).

    Returns
    -------
    numpy.ndarray
    """
    if q in ("Nu", "Nu-1") and "Nu_sys" in df:
        return df["Nu_sys"].to_numpy(float)
    return np.full(len(df), np.nan)


def _arr(v, n):
    """``v`` as a float array of length ``n`` (NaN for None)."""
    if v is None:
        return np.full(n, np.nan)
    return np.broadcast_to(np.asarray(v, float), (n,)).copy()


def errorbar(ax, x, y, err=None, sys=None, mode="all", **kw):
    """Markers (and line) with statistical bars and systematic whiskers.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes.
    x, y : float or array_like
        Points.
    err : float or array_like, optional
        Statistical errors (1 standard error); NaN or None where unknown.
    sys : float or array_like, optional
        Systematic bounds; drawn with ``mode="all"`` as thin whiskers to ±(err + sys).
    mode : {"all", "stat", "none"}, optional
        Which bars to draw (see the module docstring).
    **kw
        Passed to :meth:`matplotlib.axes.Axes.errorbar` for the markers; ``ecolor``,
        ``elinewidth``, ``capsize`` and ``capthick`` style the bars.

    Returns
    -------
    matplotlib.container.ErrorbarContainer
        The markers (for legends).
    """
    bar = {k: kw.pop(k) for k in _BAR_KEYS if k in kw}
    kw.pop("yerr", None)
    h = ax.errorbar(x, y, **kw)
    x, y = np.atleast_1d(np.asarray(x, float)), np.atleast_1d(np.asarray(y, float))
    if mode == "none" or not len(x):
        return h
    colour = bar.get("ecolor") or (h.lines[0].get_color() if h.lines[0] is not None else "0.3")
    elw = bar.get("elinewidth", kw.get("lw", 0.8) or 0.8)
    z = kw.get("zorder", 2)
    e = _arr(err, len(x))
    f = np.isfinite(e) & (e > 0)
    if f.any():
        ax.errorbar(x[f], y[f], yerr=e[f], fmt="none", ecolor=colour, elinewidth=elw,
                    capsize=bar.get("capsize", 1.5), capthick=bar.get("capthick", elw), zorder=z)
    if mode == "all":
        b = _arr(sys, len(x))
        g = np.isfinite(b) & (b > 0)
        if g.any():
            tot = np.nan_to_num(e[g]) + b[g]
            ax.vlines(x[g], y[g] - tot, y[g] + tot, colors=colour, lw=0.6 * elw, alpha=0.55, zorder=z - 0.1)
    return h


def fit_map(X, w, deg):
    """The linear map from values to weighted least-squares polynomial coefficients.

    Parameters
    ----------
    X : array_like
        Abscissa.
    w : array_like
        Weights (1 / error, or ones).
    deg : int
        Degree.

    Returns
    -------
    numpy.ndarray
        Shape (deg + 1, n), highest power first.
    """
    from rbc.analysis.exponents import fit_map as fm
    return fm(X, w, deg)


def fit_loglog(x, y, sy=None, kind="linear", bsy=None):
    """Weighted least squares of log10 y on log10 x.

    Parameters
    ----------
    x, y : array_like
        Data (positive).
    sy : array_like, optional
        Absolute 1-sigma errors of y, converted to log space as sy / (y ln 10).
        Without them the fit is unweighted and the errors are scaled by the residuals.
    kind : {"linear", "quadratic"}, optional
        ``linear``: log y = b + m log x (a power law, the default);
        ``quadratic``: log y = b + m log x + q (log x)² (a drifting exponent).
    bsy : array_like, optional
        Absolute systematic bounds of y; propagated to the coefficients as worst-case
        bounds (:func:`rbc.analysis.exponents.propagate_sys`).

    Returns
    -------
    dict
        ``kind``; ``coef`` (highest power first), ``err`` (1-sigma) and ``sys`` (NaN
        where unknown); ``dof``;
        ``chi2_dof``; ``aic`` = chi² + 2k with k parameters; ``predict``, a function
        evaluating the fit at new x.
    """
    X, Y = np.log10(x), np.log10(y)
    deg = {"linear": 1, "quadratic": 2}[kind]
    k, n = deg + 1, len(X)
    w = None if sy is None else 1.0 / (np.asarray(sy) / (np.asarray(y) * LN10))
    if n <= k:
        coef = np.polyfit(X, Y, deg, w=w)
        cov = np.full((k, k), np.nan)
    else:
        coef, cov = np.polyfit(X, Y, deg, w=w, cov="unscaled" if w is not None else True)
    r = Y - np.polyval(coef, X)
    chi2 = float(np.sum((r * (w if w is not None else 1)) ** 2))
    dof = n - k
    from rbc.analysis.exponents import propagate_sys
    bL = _arr(bsy, n) / (np.asarray(y, float) * LN10)
    sys = propagate_sys(fit_map(X, np.ones(n) if w is None else w, deg), bL)
    return dict(kind=kind, coef=coef, err=np.sqrt(np.diag(cov)), sys=sys, dof=dof,
                chi2_dof=chi2 / dof if dof > 0 else np.nan, aic=chi2 + 2 * k,
                predict=lambda xx: 10 ** np.polyval(coef, np.log10(xx)))


def slope_text(fit, mode="all") -> str:
    """Legend text for a fit: the slope m with its errors (and q for quadratics).

    Parameters
    ----------
    fit : dict
        Output of :func:`fit_loglog`.
    mode : {"all", "stat", "none"}, optional
        Errors quoted: statistical and systematic (``m ± stat ± sys``), statistical, or none.

    Returns
    -------
    str
        e.g. ``$m=-0.047\\pm0.003\\pm0.004$``; "(exact)" marks fits with no degrees of freedom.
    """
    m, dm = fit["coef"][-2], fit["err"][-2]
    bm = fit.get("sys", [np.nan] * 2)[-2]
    s = rf"$m={m:+.3f}"
    if mode != "none" and np.isfinite(dm):
        s += rf"\pm{dm:.3f}"
    if mode == "all" and np.isfinite(bm):
        s += rf"\pm{bm:.3f}"
    s += "$"
    if fit["kind"] == "quadratic":
        s += rf", $q={fit['coef'][0]:+.3f}$"
    if fit["dof"] <= 0:
        s += " (exact)"
    return s


def log_axes(ax, x=True, y=True):
    """Set log scales.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes.
    x, y : bool, optional
        Which axes to make logarithmic.
    """
    if x:
        ax.set_xscale("log")
    if y:
        ax.set_yscale("log")


def tidy_log_y(ax, max_decades=1.2):
    """Plain-number ticks on a log y axis that spans less than about a decade.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes.
    max_decades : float, optional
        Apply only if the axis spans fewer decades than this.
    """
    from matplotlib import ticker
    lo, hi = ax.get_ylim()
    if np.log10(hi / lo) < max_decades:
        ax.yaxis.set_major_locator(ticker.MaxNLocator(5))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
        ax.yaxis.set_minor_formatter(ticker.NullFormatter())


def plain_log_y(ax, max_decades=2.5):
    """Plain-number ticks (0.02, 0.05, 3) on a log y axis spanning a few decades at most.

    Below a decade the ticks are round mantissas (1-9, or 1, 2, 3, 5); up to
    ``max_decades`` the decades are labelled as plain numbers and the 2s and 5s (or 3s)
    between them too, instead of the ``4 x 10^0`` labels Matplotlib gives by default.
    Very short axes fall back to :func:`tidy_log_y`.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes.
    max_decades : float, optional
        Leave the axis alone if it spans more decades than this.
    """
    from matplotlib import ticker
    lo, hi = ax.get_ylim()
    span = np.log10(hi / lo)
    if span < 0.15:
        tidy_log_y(ax, 1.0)
        return
    if span > max_decades:
        return
    if span < 1.0:          # round mantissas: 0.04, 0.05, 0.06 rather than 0.048, 0.056
        subs = np.arange(1, 10) if span < 0.7 else (1, 2, 3, 5)
        plain = ticker.FuncFormatter(lambda v, _: f"{v:g}")
        ax.yaxis.set_major_locator(ticker.LogLocator(base=10, subs=subs))
        ax.yaxis.set_major_formatter(plain)
        ax.yaxis.set_minor_locator(ticker.LogLocator(base=10, subs=np.arange(1, 10, 0.5 if span < 0.5 else 1)))
        ax.yaxis.set_minor_formatter(ticker.NullFormatter())
        return
    plain = ticker.FuncFormatter(lambda v, _: f"{v:g}")
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10, subs=(1.0,)))
    ax.yaxis.set_major_formatter(plain)
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10, subs=np.arange(2, 10)))
    keep = (2, 5) if span < 1.8 else (3,) if span < 2.2 else ()
    ax.yaxis.set_minor_formatter(ticker.FuncFormatter(
        lambda v, _: f"{v:g}" if round(v / 10 ** np.floor(np.log10(v) + 1e-9)) in keep else ""))


def proxy(kw, ls="-"):
    """Legend handle showing a marker and a line.

    Parameters
    ----------
    kw : dict
        Marker keywords (``color``, ``marker``, ``mfc``, ``mec``, ``path_effects``).
    ls : str, optional
        Line style.

    Returns
    -------
    matplotlib.lines.Line2D
    """
    from matplotlib.lines import Line2D
    return Line2D([], [], color=kw["color"], marker=kw.get("marker"), mfc=kw.get("mfc"),
                  mec=kw.get("mec"), ls=ls, lw=1, path_effects=kw.get("path_effects", []))


def decade_ticks(ax, axis="x"):
    """A labelled major tick at every decade (Matplotlib skips some on narrow axes).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes.
    axis : {"x", "y"}, optional
        Which axis.
    """
    from matplotlib import ticker
    a = ax.xaxis if axis == "x" else ax.yaxis
    a.set_major_locator(ticker.LogLocator(base=10, numticks=30))
    a.set_minor_locator(ticker.LogLocator(base=10, subs=np.arange(2, 10) * 0.1, numticks=30))
    a.set_minor_formatter(ticker.NullFormatter())
