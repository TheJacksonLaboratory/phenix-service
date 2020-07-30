"""
This file is a slightly modified version of the `gr_metrics` code available here:
https://github.com/datarail/gr_metrics/blob/master/SRC/python/gr50/__init__.py
from commit e583df61a888cd3e8ba74bdfcf8ea13f817e67a3
"""

import pandas as pd
import numpy as np
import scipy.optimize, scipy.stats


def logistic(x, params):
    einf, log10_mid, slope = params
    emin = 1.0
    mid = 10 ** log10_mid
    return ((emin - einf) / (1 + ((x / mid) ** slope))) + einf


def _logistic_inv(y, params):
    einf, log10_mid, slope = params
    emin = 1.0
    mid = 10 ** log10_mid
    if y >= min(emin, einf) and y <= max(emin, einf):
        return mid * ((y - emin) / (einf - y)) ** (1 / slope)
    else:
        return np.inf


def _flat(x, params):
    (y,) = params
    return y


def _rss(params, fn, xdata, ydata):
    rss = 0.0
    for x, y in zip(xdata, ydata):
        rss += (y - fn(x, params)) ** 2
    return rss


def _tss(ydata):
    tss = 0.0
    y_mean = ydata.mean()
    for y in ydata:
        tss += (y - y_mean) ** 2
    return tss


def _rsquare(params, fn, xdata, ydata):
    ss_res = _rss(params, fn, xdata, ydata)
    ss_tot = _tss(ydata)
    return 1 - ss_res / ss_tot


def _fit(fn, xdata, ydata, prior, bounds):
    res = scipy.optimize.minimize(
        _rss, args=(fn, xdata, ydata), x0=prior, bounds=bounds
    )
    return res


def _calculate_pval(logistic_result, flat_result, n):
    rss2 = logistic_result.fun
    rss1 = flat_result.fun
    df1 = len(logistic_result.x) - len(flat_result.x)
    df2 = n - len(logistic_result.x) + 1
    f = ((rss1 - rss2) / df1) / (rss2 / df2)
    pval = 1 - scipy.stats.f.cdf(f, df1, df2)
    return pval


def _summarize(df, pval, alpha, logistic_result, flat_result, col):
    if pval > alpha or not logistic_result.success:
        if flat_result.x[0] > 0.5:
            Z_50 = np.inf
        else:
            Z_50 = -np.inf
        inf = flat_result.x[0]
        EZ_50 = 0.0
        slope = 0.01
        r2 = _rsquare(flat_result.x, _flat, df.concentration, df[col])
    else:
        Z_50 = _logistic_inv(0.5, logistic_result.x)
        inf = logistic_result.x[0]
        EZ_50 = 10 ** logistic_result.x[1]
        slope = logistic_result.x[2]
        r2 = _rsquare(logistic_result.x, logistic, df.concentration, df[col])
    max_ = min(df[col][-2:])
    log_conc = np.log10(df.concentration)
    aoc_width = log_conc.max() - log_conc.min()
    aoc = np.trapz(1 - df[col], log_conc) / aoc_width
    return [Z_50, max_, aoc, EZ_50, inf, slope, r2, pval]


def _metrics(df, gr_alpha=0.05, gr_value="GRvalue", ic_alpha=0.05, ic_value="t1"):
    df = df.sort_values(by="concentration")
    conc_min = df.concentration.min() / 100
    conc_max = df.concentration.max() * 100

    gr_bounds = np.array([[-1, 1], np.log10([conc_min, conc_max]), [0.1, 5]])
    ic_bounds = np.array([[0, 1], np.log10([conc_min, conc_max]), [0.1, 5]])
    prior = np.array([0.1, np.log10(np.median(df.concentration)), 2])

    gr_logistic_result = _fit(
        logistic, df.concentration, df[gr_value], prior, gr_bounds
    )
    gr_flat_result = _fit(
        _flat, df.concentration, df[gr_value], prior[[0]], gr_bounds[[0]]
    )
    gr_pval = _calculate_pval(gr_logistic_result, gr_flat_result, len(df.concentration))

    ic_logistic_result = _fit(
        logistic, df.concentration, df[ic_value], prior, ic_bounds
    )
    ic_flat_result = _fit(
        _flat, df.concentration, df[ic_value], prior[[0]], ic_bounds[[0]]
    )
    ic_pval = _calculate_pval(ic_logistic_result, ic_flat_result, len(df.concentration))

    gr_data = _summarize(
        df, gr_pval, gr_alpha, gr_logistic_result, gr_flat_result, gr_value
    )
    ic_data = _summarize(
        df, ic_pval, ic_alpha, ic_logistic_result, ic_flat_result, ic_value
    )

    return gr_data + ic_data


def gr_metrics(data, alpha=0.05, gr_value="GRvalue", ic_value="t1", keys=["Drug"]):
    gb = data.groupby(keys)
    metric_columns = [
        "GR50",
        "GRmax",
        "GR_AOC",
        "GEC50",
        "GRinf",
        "h_GR",
        "r2_GR",
        "pval_GR",
        "IC50",
        "Emax",
        "AUC",
        "EC50",
        "Einf",
        "h",
        "r2_rel_cell",
        "pval_rel_cell",
    ]
    data = [[k] + _metrics(v, alpha, gr_value, ic_value=ic_value) for k, v in gb]
    df = pd.DataFrame(data, columns=keys + metric_columns)
    return df
