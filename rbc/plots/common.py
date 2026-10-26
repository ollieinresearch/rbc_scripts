"""Fits and small helpers shared by the paper figures."""
from __future__ import annotations

import numpy as np

LN10 = np.log(10)


def fit_loglog(x, y, sy=None, kind="linear"):
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

    Returns
    -------
    dict
        ``kind``; ``coef`` (highest power first) and ``err`` (1-sigma); ``dof``;
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
    return dict(kind=kind, coef=coef, err=np.sqrt(np.diag(cov)), dof=dof,
                chi2_dof=chi2 / dof if dof > 0 else np.nan, aic=chi2 + 2 * k,
                predict=lambda xx: 10 ** np.polyval(coef, np.log10(xx)))


def slope_text(fit) -> str:
    """Legend text for a fit: the slope m with its error (and q for quadratics).

    Parameters
    ----------
    fit : dict
        Output of :func:`fit_loglog`.

    Returns
    -------
    str
        e.g. ``$m=-0.047$`` with its error; "(exact)" marks fits with no degrees of freedom.
    """
    m, dm = fit["coef"][-2], fit["err"][-2]
    s = rf"$m={m:+.3f}" + (rf"\pm{dm:.3f}$" if np.isfinite(dm) else "$")
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
        Marker keywords (``color``, ``marker``, ``mfc``, ``mec``).
    ls : str, optional
        Line style.

    Returns
    -------
    matplotlib.lines.Line2D
    """
    from matplotlib.lines import Line2D
    return Line2D([], [], color=kw["color"], marker=kw.get("marker"), mfc=kw.get("mfc"),
                  mec=kw.get("mec"), ls=ls, lw=1)


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
