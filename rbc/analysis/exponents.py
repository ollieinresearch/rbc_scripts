"""Local scaling exponents of Q = C Pr^alpha Ra^beta on the (Ra, Pr) grid.

Work in x = log10 Ra, y = log10 Pr, L = log10 Q, so alpha = ∂L/∂y and beta = ∂L/∂x.
Q is Nu - 1 by default (also Nu or Re).

edges
    Two-point secants along a row (fixed Pr) give beta, along a column (fixed Ra)
    give alpha, located at the midpoint. A secant equals the derivative at the
    midpoint up to O(h²) (a central difference), and exactly if L is quadratic.
    Error: sqrt(sL1² + sL2²) / h.
cells
    For a grid cell with corners Q00 = (x, y), Q10 = (x+h, y), Q01 = (x, y+k),
    Q11 = (x+h, y+k)::

        beta  = [(L10 - L00) + (L11 - L01)] / (2h)
        alpha = [(L01 - L00) + (L11 - L10)] / (2k)

    the average of the two edge secants on each side, located at the cell centre.
    With h = 1 decade in Ra and k = 1/2 decade in Pr, the diagonal Q01 -> Q10 is one
    step along a Pr ∝ Ra^-1/2 slice and Q00 -> Q11 one step along Pr ∝ Ra^+1/2; their
    secants are s- = (L10 - L01)/h = beta - alpha/2 and s+ = (L11 - L00)/h =
    beta + alpha/2, so beta = (s+ + s-)/2 and alpha = s+ - s-. The cell centre lies on
    the slice, and the slice exponent (what the compensated slice plots measure) is
    s- = beta - alpha/2.
lines
    Weighted log-log regression along each full row (beta) or column (alpha).
surface
    Weighted least-squares polynomial L(x, y) of total degree n; alpha and beta are
    its analytic derivatives, with errors from the coefficient covariance. n = 1 is
    the single global power law; n = 2 lets alpha and beta drift linearly (a
    crossover). Degrees are compared with AICc (AIC with the small-sample correction
    2k(k+1)/(n-k-1)).

Errors in L: sL = sQ / (Q ln 10), with sQ the standard error of the time mean of Nu
(or Re).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

LN10 = np.log(10)


def prepare(df: pd.DataFrame, q: str = "Nu-1") -> pd.DataFrame:
    """Log-space table of one quantity.

    Parameters
    ----------
    df : pandas.DataFrame
        Results rows (``log_ra``, ``log_pr``, ``Nu``, ``Nu_se``, ...).
    q : {"Nu-1", "Nu", "Re"}, optional
        Quantity Q.

    Returns
    -------
    pandas.DataFrame
        Columns ``x`` (log10 Ra), ``y`` (log10 Pr), ``L`` (log10 Q) and ``sL`` (its
        standard error, NaN if unknown); rows with Q <= 0 dropped.
    """
    if q == "Nu-1":
        Q, sQ = df["Nu"] - 1, df.get("Nu_se")
    else:
        Q, sQ = df[q], df.get(f"{q}_se")
    sQ = sQ if sQ is not None else pd.Series(np.nan, index=df.index)
    d = pd.DataFrame(dict(x=df["log_ra"].astype(float), y=df["log_pr"].astype(float),
                          L=np.log10(Q.where(Q > 0)), sL=sQ / (Q * LN10)))
    return d.dropna(subset=["L"]).reset_index(drop=True)


def _secants(g, along, other, name):
    """Two-point secants between neighbours of one row or column."""
    g = g.sort_values(along)
    a, L, s = g[along].values, g["L"].values, g["sL"].values
    out = []
    for i in range(len(a) - 1):
        h = a[i + 1] - a[i]
        out.append({"kind": name, along: 0.5 * (a[i] + a[i + 1]), other: g[other].iloc[0],
                    "value": (L[i + 1] - L[i]) / h, "err": np.hypot(s[i], s[i + 1]) / h,
                    "span": h, "source": "edge"})
    return out


def edges(d: pd.DataFrame) -> pd.DataFrame:
    """Two-point exponent estimates at the midpoints of neighbouring runs.

    Parameters
    ----------
    d : pandas.DataFrame
        Output of :func:`prepare`.

    Returns
    -------
    pandas.DataFrame
        Columns ``kind`` (alpha or beta), ``x``, ``y`` (midpoint), ``value``, ``err``,
        ``span`` (step in decades) and ``source`` (``"edge"``).
    """
    rows = []
    for _, g in d.groupby("y"):
        rows += _secants(g, "x", "y", "beta")
    for _, g in d.groupby("x"):
        rows += _secants(g, "y", "x", "alpha")
    return pd.DataFrame(rows, columns=["kind", "x", "y", "value", "err", "span", "source"])


def cells(d: pd.DataFrame) -> pd.DataFrame:
    """Both exponents at the centre of every complete grid cell.

    Parameters
    ----------
    d : pandas.DataFrame
        Output of :func:`prepare`.

    Returns
    -------
    pandas.DataFrame
        Columns ``x``, ``y`` (cell centre), ``alpha``, ``alpha_err``, ``beta``,
        ``beta_err``, ``cov`` (their covariance), ``s_slice`` = beta - alpha/2 (slope
        along the Pr ∝ Ra^-1/2 slice through the centre), ``s_slice_err`` and ``c``
        (that slice's label log10 Pr_0).
    """
    xs, ys = np.sort(d["x"].unique()), np.sort(d["y"].unique())
    look = {(r.x, r.y): (r.L, r.sL) for r in d.itertuples()}
    rows = []
    for x0, x1 in zip(xs[:-1], xs[1:]):
        for y0, y1 in zip(ys[:-1], ys[1:]):
            c = [look.get(k) for k in ((x0, y0), (x1, y0), (x0, y1), (x1, y1))]
            if any(v is None for v in c):
                continue
            (L00, s00), (L10, s10), (L01, s01), (L11, s11) = c
            h, k = x1 - x0, y1 - y0
            beta = (L10 - L00 + L11 - L01) / (2 * h)
            alpha = (L01 - L00 + L11 - L10) / (2 * k)
            v = np.array([s00, s10, s01, s11]) ** 2
            sb = np.sqrt(v.sum()) / (2 * h)
            sa = np.sqrt(v.sum()) / (2 * k)
            cov = (v[0] - v[1] - v[2] + v[3]) / (4 * h * k)
            s_minus = beta - 0.5 * alpha          # slope along Pr ~ Ra^-1/2 (dy/dx = -1/2)
            rows.append(dict(x=0.5 * (x0 + x1), y=0.5 * (y0 + y1), alpha=alpha, alpha_err=sa,
                             beta=beta, beta_err=sb, cov=cov, s_slice=s_minus,
                             s_slice_err=np.sqrt(sb ** 2 + 0.25 * sa ** 2 - cov),
                             c=0.5 * (y0 + y1) + 0.5 * (0.5 * (x0 + x1) - 5)))   # slice id log10 Pr_0
    return pd.DataFrame(rows)


def _wfit(t, v, e, deg=1):
    """Weighted polynomial fit with parameter covariance.

    Unweighted, with the covariance scaled by the residuals, when the errors are
    unknown.

    Parameters
    ----------
    t, v, e : array_like
        Abscissa, values and their standard errors.
    deg : int, optional
        Polynomial degree.

    Returns
    -------
    dict or None
        ``coef`` (highest power first), ``err``, ``chi2_dof`` and ``n``; None if there
        are too few points.
    """
    ok = np.isfinite(v)
    t, v, e = np.asarray(t)[ok], np.asarray(v)[ok], np.asarray(e)[ok]
    if len(t) <= deg:
        return None
    if np.all(np.isfinite(e)) and np.all(e > 0):
        coef, cov = np.polyfit(t, v, deg, w=1 / e, cov="unscaled")
        chi2 = float(np.sum(((v - np.polyval(coef, t)) / e) ** 2))
    else:
        if len(t) <= deg + 1:
            coef = np.polyfit(t, v, deg)
            return dict(coef=coef, err=np.full(deg + 1, np.nan), chi2_dof=np.nan, n=len(t))
        coef, cov = np.polyfit(t, v, deg, cov=True)
        chi2 = np.nan
    dof = len(t) - deg - 1
    return dict(coef=coef, err=np.sqrt(np.diag(cov)), chi2_dof=chi2 / dof if dof > 0 else np.nan, n=len(t))


def lines(d: pd.DataFrame, min_points: int = 3) -> pd.DataFrame:
    """One weighted log-log regression per full row (beta) and column (alpha).

    Parameters
    ----------
    d : pandas.DataFrame
        Output of :func:`prepare`.
    min_points : int, optional
        Rows and columns with fewer points are skipped.

    Returns
    -------
    pandas.DataFrame
        Columns ``kind``, ``x``, ``y`` (mean position), ``lo``, ``hi`` (range along the
        line), ``value``, ``err``, ``chi2_dof``, ``n`` and ``source`` (``"line"``).
    """
    rows = []
    for name, along, other in (("beta", "x", "y"), ("alpha", "y", "x")):
        for val, g in d.groupby(other):
            if len(g) < min_points:
                continue
            f = _wfit(g[along], g["L"], g["sL"], 1)
            rows.append({"kind": name, other: val, along: g[along].mean(), "lo": g[along].min(),
                         "hi": g[along].max(), "value": f["coef"][0], "err": f["err"][0],
                         "chi2_dof": f["chi2_dof"], "n": len(g), "source": "line"})
    return pd.DataFrame(rows)


def flatness(t, v, e) -> dict:
    """Is v(t) constant? Fits of a constant and of a line.

    Parameters
    ----------
    t, v, e : array_like
        Abscissa (log10 of the parameter), exponent estimates and their errors.

    Returns
    -------
    dict
        ``mean``, ``mean_err``, ``chi2_flat`` (chi² per degree of freedom of the
        constant), ``slope``, ``slope_err``, ``intercept``, ``chi2_lin`` and ``n``;
        empty if there are too few points. chi2_flat near 1 with a slope consistent
        with 0 means the exponent is constant.
    """
    c = _wfit(t, v, e, 0)
    l = _wfit(t, v, e, 1)
    if c is None or l is None:
        return {}
    return dict(mean=c["coef"][0], mean_err=c["err"][0], chi2_flat=c["chi2_dof"],
                slope=l["coef"][0], slope_err=l["err"][0], intercept=l["coef"][1],
                chi2_lin=l["chi2_dof"], n=c["n"])


class Surface:
    """Weighted least-squares polynomial L(x, y) of total degree ``deg``.

    Parameters
    ----------
    d : pandas.DataFrame
        Output of :func:`prepare`.
    deg : int
        Total degree.

    Attributes
    ----------
    coef : numpy.ndarray
        Coefficients of the monomials x^i y^j listed in ``terms``.
    cov : numpy.ndarray
        Their covariance (scaled by the residuals when the errors are unknown).
    weighted : bool
        Whether standard errors were available for every point.
    chi2, dof, rms, aic : float
        Weighted residual sum of squares, degrees of freedom, rms residual of
        log10 Q, and AICc.
    """

    def __init__(self, d: pd.DataFrame, deg: int):
        self.deg = deg
        self.terms = [(i, j) for i in range(deg + 1) for j in range(deg + 1 - i)]
        A = self._design(d["x"].values, d["y"].values)
        self.weighted = bool(np.all(np.isfinite(d["sL"])) and np.all(d["sL"] > 0))
        w = 1 / d["sL"].values if self.weighted else np.ones(len(d))
        Aw, Lw = A * w[:, None], d["L"].values * w
        self.coef, *_ = np.linalg.lstsq(Aw, Lw, rcond=None)
        r = Lw - Aw @ self.coef
        self.chi2 = float(r @ r)
        self.dof = len(d) - len(self.terms)
        self.cov = np.linalg.pinv(Aw.T @ Aw)
        self.n, k = len(d), len(self.terms)
        self.rms = float(np.sqrt(np.mean((d["L"].values - A @ self.coef) ** 2)))
        if self.weighted:
            self.aic = self.chi2 + 2 * k
        else:   # unknown errors: scale covariance by the residuals, Gaussian AIC from the RSS
            if self.dof > 0:
                self.cov *= self.chi2 / self.dof
            self.aic = self.n * np.log(max(self.chi2, 1e-300) / self.n) + 2 * k
        # small-sample correction (AICc): matters when n is not >> k, as here
        self.aic += 2 * k * (k + 1) / (self.n - k - 1) if self.n - k - 1 > 0 else np.inf

    def _design(self, x, y):
        """Design matrix: one column per monomial."""
        return np.stack([x ** i * y ** j for i, j in self.terms], axis=-1)

    def _grad_rows(self, x, y, wrt):
        """Rows mapping the coefficients to ∂L/∂x or ∂L/∂y at (x, y)."""
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        cols = []
        for i, j in self.terms:
            if wrt == "x":
                cols.append(i * x ** max(i - 1, 0) * y ** j if i else np.zeros_like(x))
            else:
                cols.append(j * x ** i * y ** max(j - 1, 0) if j else np.zeros_like(x))
        return np.stack(cols, axis=-1)

    def exponent(self, x, y, kind):
        """Local exponent of the fitted surface.

        Parameters
        ----------
        x, y : float or array_like
            log10 Ra and log10 Pr.
        kind : {"alpha", "beta"}
            ∂L/∂y or ∂L/∂x.

        Returns
        -------
        value, err : numpy.ndarray
            The exponent and its standard error from the coefficient covariance.
        """
        G = self._grad_rows(x, y, "y" if kind == "alpha" else "x")
        val = G @ self.coef
        err = np.sqrt(np.einsum("...i,ij,...j->...", G, self.cov, G))
        return val, err


def surfaces(d: pd.DataFrame, degs=(1, 2, 3)) -> list[Surface]:
    """Fit :class:`Surface` for several degrees.

    Parameters
    ----------
    d : pandas.DataFrame
        Output of :func:`prepare`.
    degs : iterable of int, optional
        Degrees; those with at least as many coefficients as points are skipped.

    Returns
    -------
    list of Surface
    """
    return [Surface(d, n) for n in degs if len(d) > (n + 1) * (n + 2) // 2]
