"""
Detection of t0, ti, tc and computation of SusI, ti/tc, SupI and VI.

This is a faithful re-implementation of the algorithm described in the
Supplementary Methods of the SusI manuscript, written so that every threshold
is an explicit, overridable argument. That is what makes the sensitivity
analysis requested by Reviewer 2 possible: the same code runs the published
defaults and any perturbation of them.

Published defaults (Supplementary Methods / Table S1)
-----------------------------------------------------
smoothing              3-point moving average
t0 search window       20 h, starting 0.25 h after the start of measurement
t0 acceptance          a local maximum at which >= 60% of the following 10
                       points are lower
t0 monotonic shortcut  if >= 80% of the first 10 points after the search start
                       are decreasing, t0 = first time point
ti rule                first point after the post-t0 minimum at which >= 60% of
                       the following 7 points are higher AND the mean increase
                       over those 7 points is >= 0.008 OD
ti cap                 24 h if no ti is detected
tc rule                first point at which the control growth rate falls below
                       0.03 /h and stays below for >= 60% of the following 5
                       points; search starts 1.0 h after the start
tc/VI minimum time     1.0 h by default, 2.0 h for the PhiWec430 experiments

Open point to confirm against the author's R code
-------------------------------------------------
The manuscript states that the tc growth rate is obtained "using the central
difference method" and compared against a threshold of 0.03 h^-1. Two readings
are possible: dOD/dt (units OD/h, often written h^-1 informally) or the
specific growth rate d(ln OD)/dt (units h^-1 strictly). Both are implemented;
`rate_mode` selects between them. The default here is 'od' because 0.03 OD/h is
the value that behaves sensibly on curves plateauing near OD 1. This must be
checked against determine_tc() in susi_analysis.Rmd before any simulated result
is quoted alongside the experimental values.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np


# --------------------------------------------------------------------------
# Algorithm settings
# --------------------------------------------------------------------------

@dataclass
class DetectCfg:
    smooth_window: int = 3

    # t0
    t0_search_start: float = 0.25     # h
    t0_search_span: float = 20.0      # h
    t0_lookahead: int = 10            # points
    t0_frac_lower: float = 0.60       # fraction of look-ahead points that must be lower
    t0_mono_frac: float = 0.80        # monotonic-decrease shortcut threshold
    t0_mono_points: int = 10

    # ti
    ti_lookahead: int = 7             # points
    ti_frac_higher: float = 0.60
    ti_min_rise: float = 0.008        # OD, mean rise over the look-ahead points
    ti_cap: float = 24.0              # h

    # tc
    tc_min_time: float = 1.0          # h  (2.0 for the PhiWec430 experiments)
    tc_rate_thresh: float = 0.03      # per h
    tc_lookahead: int = 5             # points
    tc_frac_below: float = 0.60
    rate_mode: str = "od"             # 'od' -> dOD/dt, 'ln' -> d(ln OD)/dt

    # VI / SupI
    vi_min_time: float = 1.0          # h
    supi_window: float = 30.0         # h


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

def moving_average(y: np.ndarray, window: int = 3) -> np.ndarray:
    """Centred moving average with edge padding, preserving length."""
    if window <= 1:
        return np.asarray(y, dtype=float)
    pad = window // 2
    ypad = np.pad(np.asarray(y, dtype=float), pad, mode="edge")
    kernel = np.ones(window) / window
    return np.convolve(ypad, kernel, mode="valid")


def _central_rate(t: np.ndarray, y: np.ndarray, mode: str = "od") -> np.ndarray:
    v = np.log(np.clip(y, 1e-9, None)) if mode == "ln" else y
    r = np.gradient(v, t)          # central differences in the interior
    return r


# --------------------------------------------------------------------------
# Time-point detection
# --------------------------------------------------------------------------

def detect_t0(t, od, cfg: DetectCfg = DetectCfg()):
    """Lysis initiation time. Returns (t0, method) with method in
    {'monotonic', 'peak', 'failed'}."""
    y = moving_average(od, cfg.smooth_window)
    i_start = int(np.searchsorted(t, cfg.t0_search_start))
    i_stop = int(np.searchsorted(t, cfg.t0_search_start + cfg.t0_search_span))
    i_stop = min(i_stop, len(t) - 1)
    if i_start >= i_stop:
        return np.nan, "failed"

    # monotonic-decrease shortcut
    seg = y[i_start:i_start + cfg.t0_mono_points + 1]
    if len(seg) >= 2:
        dec = np.mean(np.diff(seg) < 0)
        if dec >= cfg.t0_mono_frac:
            return float(t[0]), "monotonic"

    # sequential peak detection
    for i in range(max(i_start, 1), i_stop):
        if not (y[i] > y[i - 1] and y[i] > y[i + 1]):
            continue
        ahead = y[i + 1:i + 1 + cfg.t0_lookahead]
        if len(ahead) == 0:
            continue
        if np.mean(ahead < y[i]) >= cfg.t0_frac_lower:
            return float(t[i]), "peak"
    return np.nan, "failed"


def detect_ti(t, od, t0, cfg: DetectCfg = DetectCfg()):
    """Turbidity-regrowth onset time (referred to as resistance emergence time
    in the original manuscript). Returns (ti, method) with method in
    {'detected', 'capped', 'failed'}."""
    if not np.isfinite(t0):
        return np.nan, "failed"
    y = moving_average(od, cfg.smooth_window)
    i0 = int(np.searchsorted(t, t0))
    i_cap = int(np.searchsorted(t, cfg.ti_cap))
    i_cap = min(i_cap, len(t) - 1)
    if i0 >= i_cap:
        return float(cfg.ti_cap), "capped"

    i_min = i0 + int(np.argmin(y[i0:i_cap + 1]))
    for i in range(i_min, i_cap + 1):
        ahead = y[i + 1:i + 1 + cfg.ti_lookahead]
        if len(ahead) < cfg.ti_lookahead:
            break
        frac_higher = np.mean(ahead > y[i])
        mean_rise = float(np.mean(ahead - y[i]))
        if frac_higher >= cfg.ti_frac_higher and mean_rise >= cfg.ti_min_rise:
            return float(t[i]), "detected"
    return float(cfg.ti_cap), "capped"


def detect_tc(t, od_control, cfg: DetectCfg = DetectCfg()):
    """Time at which the control culture reaches stationary phase."""
    y = moving_average(od_control, cfg.smooth_window)
    r = _central_rate(t, y, cfg.rate_mode)
    i_start = int(np.searchsorted(t, cfg.tc_min_time))
    for i in range(i_start, len(t)):
        if r[i] >= cfg.tc_rate_thresh:
            continue
        ahead = r[i + 1:i + 1 + cfg.tc_lookahead]
        if len(ahead) == 0:
            return float(t[i])
        if np.mean(ahead < cfg.tc_rate_thresh) >= cfg.tc_frac_below:
            return float(t[i])
    return float(t[-1])


# --------------------------------------------------------------------------
# Indices
# --------------------------------------------------------------------------

def _auc(t, y, a, b):
    """Trapezoidal integral of y over [a, b] with linear interpolation at the
    end points, matching the trapezoidal rule stated in the manuscript."""
    if not (np.isfinite(a) and np.isfinite(b)) or b <= a:
        return 0.0
    tt = np.concatenate(([a], t[(t > a) & (t < b)], [b]))
    yy = np.interp(tt, t, y)
    return float(np.trapezoid(yy, tt))


def susi(t, od_phage, od_control, t0, ti, tc):
    """SusI as defined in the manuscript, plus its exact decomposition.

    SusI = (integral of C from t0 to ti - integral of P from t0 to ti)
           / (integral of C from 0 to tc)

    The decomposition below is an algebraic identity, obtained by writing each
    integral as (interval length) x (time-average over that interval):

        SusI = (ti - t0)/tc  x  <C - P>_[t0,ti] / <C>_[0,tc]

    The first factor is the suppression duration relative to one control growth
    cycle; the second is the mean suppression depth relative to the mean control
    density. Note that the first factor is (ti - t0)/tc and NOT ti/tc: the two
    differ by t0/tc.
    """
    if not all(np.isfinite(v) for v in (t0, ti, tc)) or tc <= 0:
        return dict(SusI=np.nan, duration_term=np.nan, depth_term=np.nan,
                    numerator=np.nan, denom=np.nan)
    num_c = _auc(t, od_control, t0, ti)
    num_p = _auc(t, od_phage, t0, ti)
    den = _auc(t, od_control, 0.0, tc)
    if den <= 0:
        return dict(SusI=np.nan, duration_term=np.nan, depth_term=np.nan,
                    numerator=np.nan, denom=np.nan)
    value = (num_c - num_p) / den
    dur = (ti - t0) / tc
    mean_supp = (num_c - num_p) / (ti - t0) if ti > t0 else np.nan
    mean_ctrl = den / tc
    depth = mean_supp / mean_ctrl if mean_ctrl > 0 else np.nan
    return dict(SusI=value, duration_term=dur, depth_term=depth,
                numerator=num_c - num_p, denom=den)


def supi(t, od_phage, od_control, cfg: DetectCfg = DetectCfg()):
    """Suppression Index: percentage growth suppression over a fixed window."""
    a_c = _auc(t, od_control, 0.0, cfg.supi_window)
    a_p = _auc(t, od_phage, 0.0, cfg.supi_window)
    if a_c <= 0:
        return np.nan
    return 100.0 * (a_c - a_p) / a_c


def local_virulence(t, od_phage, od_control, tc, cfg: DetectCfg = DetectCfg()):
    """Local virulence vi = 1 - A_i/A_0 integrated from vi_min_time to tc."""
    a0 = _auc(t, od_control, cfg.vi_min_time, tc)
    ai = _auc(t, od_phage, cfg.vi_min_time, tc)
    if a0 <= 0:
        return np.nan
    return 1.0 - ai / a0


def virulence_index(log_mois, vis):
    """VI = area under the local-virulence/log(MOI) curve, normalised by the
    theoretical maximum, following Storms et al. (2020)."""
    lm = np.asarray(log_mois, dtype=float)
    v = np.asarray(vis, dtype=float)
    order = np.argsort(lm)
    lm, v = lm[order], v[order]
    ap = float(np.trapezoid(v, lm))
    amax = lm[-1] - lm[0]
    return ap / amax if amax > 0 else np.nan


# --------------------------------------------------------------------------
# One-stop analysis of a phage/control curve pair
# --------------------------------------------------------------------------

def analyse_pair(t, od_phage, od_control, cfg: DetectCfg = DetectCfg()):
    t0, t0_method = detect_t0(t, od_phage, cfg)
    ti, ti_method = detect_ti(t, od_phage, t0, cfg)
    tc = detect_tc(t, od_control, cfg)
    s = susi(t, od_phage, od_control, t0, ti, tc)
    out = dict(t0=t0, ti=ti, tc=tc,
               t0_method=t0_method, ti_method=ti_method,
               ti_over_tc=(ti / tc if np.isfinite(ti) and tc > 0 else np.nan),
               SupI=supi(t, od_phage, od_control, cfg),
               local_vi=local_virulence(t, od_phage, od_control, tc, cfg))
    out.update(s)
    return out
