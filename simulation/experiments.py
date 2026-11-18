"""
Simulation experiments answering specific reviewer comments on the SusI
manuscript.

Organising result
-----------------
Two dimensionless quantities control the two factors of SusI:

    SusI = (ti - t0)/tc  x  <C - P>/<C>
            duration            depth

* the duration factor is set by the frequency of phage-insensitive cells, f;
* the depth factor is set by the phage's effective reproductive number in the
  culture,

      R_phage = beta * (delta * N0) / (delta * N0 + omega)

  i.e. burst size times the probability that a free phage adsorbs to a host
  before it is lost. Depth saturates for R_phage greater than roughly 5.

This gives a mechanistic account of when SusI and ti/tc coincide and when they
separate: in a panel of strongly amplifying phages (R_phage >> 5) the depth
factor is pinned near its maximum, SusI reduces to the duration factor, and the
two metrics are nearly collinear. SusI carries extra information only when the
panel includes marginally amplifying phages, R_phage of order 1 to 5, which are
exactly the phages that produce incomplete ("blunt") lysis curves.

Experiments
-----------
D1  Matched suppression duration, different suppression depth.
    Phage archetypes spanning R_phage, each with f tuned so that the
    suppression duration matches a common target. ti/tc is then effectively
    identical across archetypes while SusI differs severalfold.
    Addresses: Reviewer 2, point 1 (what does the integrated area add).

D2  Map of the (R_phage, f) plane, and the induced (ti/tc, SusI) plane.
    Addresses: Reviewer 2, point 1 (robustness, added value).

D3  Resistance frequency against SusI, including the two-receptor cocktail
    case where the double-mutant frequency is the product of the singles.
    Addresses: Reviewer 1, point 1 (impact of resistance frequency on SusI).

D4  Sensitivity of phage ranking to the algorithm's detection thresholds.
    Addresses: Reviewer 2, point 1 (empirically chosen thresholds) and
    minor point 2 (24 h cap against 48 h measurement).

Usage:  python3 experiments.py [d1|d2|d3|d4|all]
"""

from __future__ import annotations

import sys
from dataclasses import asdict, replace
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from model import Params, simulate_pair
from detect import DetectCfg, analyse_pair

CFG = DetectCfg()


# --------------------------------------------------------------------------
# Phage archetypes
# --------------------------------------------------------------------------

def r_phage(p: Params) -> float:  # kept for backward compatibility
    """Effective reproductive number of the phage in the culture at inoculation."""
    aN = p.delta * p.N0
    return p.beta * aN / (aN + p.omega)


#: archetypes chosen to span the informative range of R_phage while keeping
#: every parameter inside a biologically defensible interval. omega is a lumped
#: phage-loss term (decay, aggregation, adsorption to debris and to insensitive
#: cells), not pure thermal inactivation, which is why values above a few per
#: hour are admissible.
ARCHETYPES = {
    "sharp":        dict(delta=1.0e-7, beta=150.0, tau=0.45, omega=0.05),
    "strong":       dict(delta=3.0e-8, beta=60.0,  tau=0.55, omega=0.50),
    "marginal":     dict(delta=3.0e-9, beta=50.0,  tau=0.70, omega=2.00),
    "blunt":        dict(delta=1.0e-9, beta=50.0,  tau=0.80, omega=4.00),
    "very blunt":   dict(delta=1.0e-9, beta=40.0,  tau=0.90, omega=6.00),
}


def make(archetype: str, **over) -> Params:
    kw = dict(ARCHETYPES[archetype])
    kw.update(over)
    return Params(**kw)


def measure(p: Params, cfg: DetectCfg = CFG) -> dict:
    t, odp, odc = simulate_pair(p)
    r = analyse_pair(t, odp, odc, cfg)
    m = (t >= r["t0"]) & (t <= r["ti"]) if np.isfinite(r["t0"]) else np.zeros_like(t, bool)
    r["mean_OD_supp"] = float(np.mean(odp[m])) if m.any() else np.nan
    r["min_OD"] = float(np.min(odp))
    r["R_phage"] = r_phage(p)
    r["f"] = p.f
    return r


# --------------------------------------------------------------------------
# D1
# --------------------------------------------------------------------------

def tune_f_for_duration(archetype: str, target: float,
                        lo: float = -9.0, hi: float = -2.5,
                        tol: float = 0.02, iters: int = 40,
                        cfg: DetectCfg = CFG):
    """Bisect log10(f) so that the suppression duration hits `target`.

    Duration falls monotonically with f, so the bisection is well posed. Returns
    (log10_f, measurement) or (None, bracket) if the target is unreachable.
    """
    d_lo = measure(make(archetype, f=10.0 ** lo), cfg)["duration_term"]
    d_hi = measure(make(archetype, f=10.0 ** hi), cfg)["duration_term"]
    if not np.isfinite(d_lo) or not np.isfinite(d_hi):
        return None, (d_lo, d_hi)
    if not (min(d_lo, d_hi) <= target <= max(d_lo, d_hi)):
        return None, (d_lo, d_hi)
    a, b, best = lo, hi, None
    for _ in range(iters):
        m = 0.5 * (a + b)
        r = measure(make(archetype, f=10.0 ** m), cfg)
        d = r["duration_term"]
        best = (m, r)
        if np.isfinite(d) and abs(d - target) < tol:
            break
        if (not np.isfinite(d)) or d > target:
            a = m
        else:
            b = m
    return best[0], best[1]


def d1_matched_duration(target: float = 1.10, archetypes=None):
    archetypes = archetypes or ["sharp", "strong", "marginal", "blunt"]
    recs = []
    for name in archetypes:
        log_f, res = tune_f_for_duration(name, target)
        if log_f is None:
            recs.append(dict(archetype=name, status=f"unreachable {res}"))
            continue
        recs.append(dict(archetype=name, status="ok", f=10.0 ** log_f,
                         R_phage=res["R_phage"], t0=res["t0"], ti=res["ti"],
                         tc=res["tc"], ti_over_tc=res["ti_over_tc"],
                         duration=res["duration_term"], depth=res["depth_term"],
                         SusI=res["SusI"], SupI=res["SupI"],
                         mean_OD_supp=res["mean_OD_supp"], min_OD=res["min_OD"]))
    return pd.DataFrame(recs)


# --------------------------------------------------------------------------
# D2
# --------------------------------------------------------------------------

def d2_map(n_f: int = 9):
    fs = np.logspace(-8.5, -3.0, n_f)
    recs = []
    for name in ARCHETYPES:
        for f in fs:
            p = make(name, f=f)
            try:
                r = measure(p)
            except RuntimeError:
                continue
            recs.append(dict(archetype=name, R_phage=r["R_phage"], f=f,
                             log10_f=np.log10(f), log10_R=np.log10(r["R_phage"]),
                             t0=r["t0"], ti=r["ti"], tc=r["tc"],
                             ti_method=r["ti_method"], ti_over_tc=r["ti_over_tc"],
                             SusI=r["SusI"], duration=r["duration_term"],
                             depth=r["depth_term"], SupI=r["SupI"],
                             mean_OD_supp=r["mean_OD_supp"]))
    return pd.DataFrame(recs)


# --------------------------------------------------------------------------
# D3
# --------------------------------------------------------------------------

def d3_frequency(archetypes=("sharp", "blunt"), n_f: int = 13):
    fs = np.logspace(-9, -3, n_f)
    recs = []
    for name in archetypes:
        for f in fs:
            r = measure(make(name, f=f))
            p = make(name, f=f)
            recs.append(dict(archetype=name, f=f, log10_f=np.log10(f),
                             founders_per_well=p.founders_per_well(),
                             t0=r["t0"], ti=r["ti"], ti_method=r["ti_method"],
                             ti_over_tc=r["ti_over_tc"], SusI=r["SusI"],
                             duration=r["duration_term"], depth=r["depth_term"],
                             SupI=r["SupI"]))
    return pd.DataFrame(recs)


def d3_cocktail(f_singles=(1e-4, 1e-5, 1e-6, 1e-7), archetype="strong"):
    """Single phage at f, same-receptor pair at f, different-receptor pair at f^2.

    A cocktail of phages binding independent receptors requires a double mutant,
    whose frequency is the product of the single-mutant frequencies. Representing
    the cocktail as one lytic agent carrying the product frequency is the minimal
    assumption that isolates the resistance-frequency effect from any change in
    lytic kinetics, so the resulting SusI increase is attributable to receptor
    independence alone.
    """
    recs = []
    for f1 in f_singles:
        for case, f_eff in (("single phage", f1),
                            ("same-receptor cocktail", f1),
                            ("different-receptor cocktail", f1 * f1)):
            p = make(archetype, f=f_eff)
            r = measure(p)
            recs.append(dict(case=case, f_single=f1, f_effective=f_eff,
                             founders_per_well=p.founders_per_well(),
                             ti=r["ti"], ti_method=r["ti_method"],
                             ti_over_tc=r["ti_over_tc"], SusI=r["SusI"],
                             duration=r["duration_term"], depth=r["depth_term"]))
    return pd.DataFrame(recs)


# --------------------------------------------------------------------------
# D4  threshold sensitivity
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# D5  partial suppression: lysis that does not go to completion
# --------------------------------------------------------------------------

#: Base condition for D5. A sub-saturating MOI with a slowly adsorbing phage
#: lets the culture grow to OD ~0.8 before the phage epidemic sweeps through,
#: so a substantial amount of biomass is lysed at once. The fraction of that
#: biomass that keeps scattering light (unlysed survivors frozen by resource
#: exhaustion, or debris from lysed cells; OD cannot distinguish the two) then
#: sets the level at which the decline bottoms out. Cryptic growth on the
#: contents of lysed cells (recycle = 0.35) is enabled so that regrowth after
#: ti is clearly visible; without it, lysis at high density leaves the
#: insensitive subpopulation almost no resource and the rise after ti is too
#: small to read off the plot.
D5_BASE = dict(moi=1e-5, delta=8.0e-9, beta=100.0, f=1.0e-7, omega=0.05,
               recycle=0.35)


def d5_partial_suppression(residuals=(0.0, 0.15, 0.30, 0.45, 0.60, 0.75),
                           base=None, cfg: DetectCfg = CFG):
    """Curves in which lysis clearly occurs but does not go to completion.

    This is the case the metric was designed for: between lysis onset (t0) and
    regrowth onset (ti), turbidity normally falls close to the blank, but for
    some phages it stops part-way down and stays there.

    Varying only the residual scattering fraction leaves t0, ti and tc
    untouched, so ti/tc is identical across the whole family by construction
    and every difference in SusI is carried by the depth term. That is the
    cleanest available demonstration that the integrated area is not redundant
    with ti/tc.

    Note for honesty in the write-up: SupI does track this difference (it falls
    from about 92% to about 42% across the family), so this scenario is an
    argument against ti/tc specifically, not against SupI. The manuscript's
    case for SusI over SupI rests on the two scenarios it already reports
    (delayed lysis onset, and secondary lysis after regrowth).
    """
    base = dict(D5_BASE if base is None else base)
    recs = []
    for kd in residuals:
        p = Params(kappa_debris=kd, **base)
        t, odp, odc = simulate_pair(p)
        r = analyse_pair(t, odp, odc, cfg)
        i0 = int(np.searchsorted(t, r["t0"]))
        od_at_t0 = float(np.median(odp[max(0, i0 - 1):i0 + 2]))
        m = (t >= r["t0"]) & (t <= r["ti"])
        plateau = float(np.median(odp[m])) if m.any() else np.nan
        recs.append(dict(residual_fraction=kd, t0=r["t0"], ti=r["ti"],
                         tc=r["tc"], ti_over_tc=r["ti_over_tc"],
                         od_at_t0=od_at_t0, plateau_OD=plateau,
                         plateau_pct_of_onset=100 * plateau / od_at_t0,
                         SusI=r["SusI"], duration=r["duration_term"],
                         depth=r["depth_term"], SupI=r["SupI"]))
    return pd.DataFrame(recs)


def d5_curves(residuals=(0.0, 0.30, 0.60), base=None):
    """Return (t, od_control, {label: od_phage}) for plotting."""
    base = dict(D5_BASE if base is None else base)
    curves, odc, tt = {}, None, None
    for kd in residuals:
        p = Params(kappa_debris=kd, **base)
        t, odp, oc = simulate_pair(p)
        tt, odc = t, oc
        curves[kd] = odp
    return tt, odc, curves


def _panel():
    """A reference panel spanning both controlling axes, used for ranking tests."""
    panel = []
    for name in ARCHETYPES:
        for f in (1e-7, 1e-6, 1e-5, 1e-4):
            panel.append((f"{name} f={f:.0e}", make(name, f=f)))
    return panel


def d4_threshold_sensitivity():
    panel = _panel()
    variants = {
        "published defaults": DetectCfg(),
        "ti_frac_higher 0.50": DetectCfg(ti_frac_higher=0.50),
        "ti_frac_higher 0.70": DetectCfg(ti_frac_higher=0.70),
        "ti_min_rise 0.004": DetectCfg(ti_min_rise=0.004),
        "ti_min_rise 0.016": DetectCfg(ti_min_rise=0.016),
        "ti_lookahead 5": DetectCfg(ti_lookahead=5),
        "ti_lookahead 10": DetectCfg(ti_lookahead=10),
        "t0_frac_lower 0.50": DetectCfg(t0_frac_lower=0.50),
        "t0_frac_lower 0.70": DetectCfg(t0_frac_lower=0.70),
        "ti cap 30 h": DetectCfg(ti_cap=30.0),
        "ti cap 36 h": DetectCfg(ti_cap=36.0),
        "ti cap 48 h": DetectCfg(ti_cap=48.0),
        "tc rate 0.02": DetectCfg(tc_rate_thresh=0.02),
        "tc rate 0.05": DetectCfg(tc_rate_thresh=0.05),
        "tc rate mode ln": DetectCfg(rate_mode="ln"),
        "smoothing 5-point": DetectCfg(smooth_window=5),
    }
    base_name = "published defaults"
    results, wide = [], {}
    for vname, cfg in variants.items():
        vals = {}
        for label, p in panel:
            try:
                vals[label] = measure(p, cfg)["SusI"]
            except RuntimeError:
                vals[label] = np.nan
        wide[vname] = vals

    base = pd.Series(wide[base_name])
    for vname, vals in wide.items():
        s = pd.Series(vals)
        ok = base.notna() & s.notna()
        rho = spearmanr(base[ok], s[ok]).statistic if ok.sum() > 2 else np.nan
        results.append(dict(variant=vname, n=int(ok.sum()),
                            spearman_vs_default=rho,
                            median_SusI=float(np.nanmedian(s)),
                            max_abs_change=float(np.nanmax(np.abs(s[ok] - base[ok]))),
                            n_missing=int(s.isna().sum())))
    return pd.DataFrame(results), pd.DataFrame(wide)


# --------------------------------------------------------------------------

def main(which="all"):
    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 60)
    ff = lambda v: f"{v:.4g}"

    if which in ("d1", "all"):
        print("\n=== D1  matched suppression duration, different depth ===")
        for target in (0.90, 1.10, 1.30):
            df = d1_matched_duration(target=target)
            print(f"\n-- target (ti-t0)/tc = {target}")
            print(df.to_string(index=False, float_format=ff))
            df.to_csv(f"out/d1_target_{target}.csv", index=False)

    if which in ("d2", "all"):
        print("\n=== D2  (R_phage, f) map ===")
        df = d2_map()
        df.to_csv("out/d2_map.csv", index=False)
        ok = df.dropna(subset=["SusI", "ti_over_tc", "depth", "duration"])
        print(f"{len(ok)} conditions")
        print("\nSpearman correlations:")
        for a, b in [("log10_R", "depth"), ("log10_R", "duration"),
                     ("log10_f", "depth"), ("log10_f", "duration"),
                     ("ti_over_tc", "SusI")]:
            print(f"  {a:12s} vs {b:10s} rho = {spearmanr(ok[a], ok[b]).statistic:+.3f}")
        print("\nSusI vs ti/tc correlation, by archetype:")
        for name, g in ok.groupby("archetype"):
            if len(g) > 2:
                print(f"  {name:12s} r = {np.corrcoef(g['ti_over_tc'], g['SusI'])[0,1]:.3f}"
                      f"   mean depth = {g['depth'].mean():.3f}  (n={len(g)})")
        strong = ok[ok.R_phage > 5]
        if len(strong) > 2:
            print(f"\n  restricted to R_phage > 5 (strongly amplifying panel): "
                  f"r = {np.corrcoef(strong['ti_over_tc'], strong['SusI'])[0,1]:.3f}  n={len(strong)}")
        allp = ok
        print(f"  whole panel including marginal phages: "
              f"r = {np.corrcoef(allp['ti_over_tc'], allp['SusI'])[0,1]:.3f}  n={len(allp)}")

    if which in ("d3", "all"):
        print("\n=== D3  resistance frequency against SusI ===")
        df = d3_frequency()
        print(df.to_string(index=False, float_format=ff))
        df.to_csv("out/d3_frequency.csv", index=False)
        print("\n=== D3b  cocktail product rule ===")
        df2 = d3_cocktail()
        print(df2.to_string(index=False, float_format=ff))
        df2.to_csv("out/d3_cocktail.csv", index=False)

    if which in ("d5", "all"):
        print("\n=== D5  partial suppression at matched t0, ti and tc ===")
        df = d5_partial_suppression()
        print(df.to_string(index=False, float_format=ff))
        df.to_csv("out/d5_partial_suppression.csv", index=False)

    if which in ("d4", "all"):
        print("\n=== D4  threshold sensitivity (ranking stability) ===")
        summ, wide = d4_threshold_sensitivity()
        print(summ.to_string(index=False, float_format=ff))
        summ.to_csv("out/d4_sensitivity_summary.csv", index=False)
        wide.to_csv("out/d4_sensitivity_values.csv")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "all")
