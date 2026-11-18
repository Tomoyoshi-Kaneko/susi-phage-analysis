"""
Well-level stochastic simulation of SusI.

Why this exists
---------------
The deterministic model seeds the phage-insensitive subpopulation as a
continuous concentration f * N0. In a 200 uL well holding 1.6e7 cells that is
not a safe approximation at low frequency: f = 1e-7 corresponds to an expected
1.6 insensitive founder cells per well, and f = 1e-8 to 0.16, so a substantial
fraction of wells contain none at all and cannot regrow by that route. In that
regime whether and when turbidity recovers is a chance event, and replicate
wells of the same condition give genuinely different answers.

This matters for the manuscript in three ways:

1. Reviewer 1 asked what effect phage resistance frequency has on SusI. The
   honest answer at low frequency is a distribution, not a number.
2. Figure S2 already colour-codes wells by detection outcome and reports the
   number of successful t0, ti and SusI determinations per panel, so
   well-to-well variability is visible in the published data. A deterministic
   model leaves that as unexplained noise; a stochastic one predicts it.
3. It removes a limitation of the deterministic model that would otherwise have
   to be conceded without an answer.

Method
------
Each well draws its number of insensitive founder cells from a Poisson
distribution with mean f * N0 * well_volume, and that integer is converted back
to a concentration for the ODE. Everything else is as in `model.py`.
Measurement noise is drawn independently per well.

The structure follows the batch-culture design used for the individual-based
stochastic model of population establishment in ammonia-oxidizing bacteria by
Ikeda, Fujitani and Tsuneda (bioRxiv 2026, doi 10.64898/2026.05.07.723170):
a fixed number of wells per simulated experiment, many independent experiments,
and results reported as a median with a prediction interval rather than a
single trajectory. Their wells-per-experiment count (12) happens to match the
technical replicate count used here.

Efficiency note
---------------
Given the parameters, the noise-free trajectory depends only on the integer
founder count, so trajectories are computed once per distinct count and cached;
independent measurement noise and the detection algorithm are then applied per
well. This turns ~1200 ODE solves per condition into a few dozen.
"""

from __future__ import annotations

from dataclasses import replace
import numpy as np
import pandas as pd

from model import Params, simulate
from detect import DetectCfg, analyse_pair

CFG = DetectCfg()


# --------------------------------------------------------------------------

def _true_curves(p: Params, founders: int):
    """Noise-free phage and control curves for a given integer founder count."""
    R0_conc = founders / p.well_volume          # cells per mL
    f_equiv = R0_conc / p.N0 if p.N0 > 0 else 0.0
    pp = replace(p, f=f_equiv, noise_sd=0.0, seed=None)
    t, _, od_phage, _ = simulate(pp, with_phage=True)
    _, _, od_ctrl, _ = simulate(pp, with_phage=False)
    return t, od_phage, od_ctrl


def simulate_wells(p: Params, n_wells: int = 12, n_experiments: int = 100,
                   cfg: DetectCfg = CFG, seed: int = 0,
                   max_founders_cached: int = 400):
    """Run n_experiments independent experiments of n_wells wells each.

    Returns a long DataFrame, one row per well.
    """
    rng = np.random.default_rng(seed)
    lam = p.f * p.N0 * p.well_volume
    counts = rng.poisson(lam, size=(n_experiments, n_wells))

    cache = {}
    rows = []
    for e in range(n_experiments):
        for w in range(n_wells):
            k = int(counts[e, w])
            key = min(k, max_founders_cached)
            if key not in cache:
                cache[key] = _true_curves(p, key)
            t, od_p_true, od_c_true = cache[key]
            od_p = od_p_true + rng.normal(0.0, p.noise_sd, size=t.shape)
            od_c = od_c_true + rng.normal(0.0, p.noise_sd, size=t.shape)
            r = analyse_pair(t, od_p, od_c, cfg)
            rows.append(dict(experiment=e, well=w, founders=k,
                             t0=r["t0"], ti=r["ti"], tc=r["tc"],
                             t0_method=r["t0_method"], ti_method=r["ti_method"],
                             ti_over_tc=r["ti_over_tc"], SusI=r["SusI"],
                             duration=r["duration_term"], depth=r["depth_term"],
                             SupI=r["SupI"]))
    df = pd.DataFrame(rows)
    df.attrs["lambda"] = lam
    df.attrs["n_cached_curves"] = len(cache)
    return df


# --------------------------------------------------------------------------

def summarise(df: pd.DataFrame) -> dict:
    """Well-level and experiment-level summary of one condition."""
    s = df.SusI.dropna()
    per_exp = df.groupby("experiment").SusI.mean()
    return dict(
        n_wells=len(df),
        lam=df.attrs.get("lambda", np.nan),
        p_zero_founders=float((df.founders == 0).mean()),
        p_ti_capped=float((df.ti_method == "capped").mean()),
        p_t0_failed=float((df.t0_method == "failed").mean()),
        SusI_median=float(s.median()),
        SusI_lo95=float(s.quantile(0.025)),
        SusI_hi95=float(s.quantile(0.975)),
        SusI_iqr=float(s.quantile(0.75) - s.quantile(0.25)),
        SusI_cv=float(s.std() / s.mean()) if s.mean() else np.nan,
        ti_median=float(df.ti.median()),
        ti_lo95=float(df.ti.quantile(0.025)),
        ti_hi95=float(df.ti.quantile(0.975)),
        exp_mean_spread=float(per_exp.max() - per_exp.min()),
    )


def frequency_scan(archetype_kw: dict, fs=None, n_wells: int = 12,
                   n_experiments: int = 100, cfg: DetectCfg = CFG,
                   seed: int = 0) -> pd.DataFrame:
    """Stochastic counterpart of the deterministic resistance-frequency sweep."""
    fs = fs if fs is not None else np.logspace(-8.5, -4.0, 10)
    recs = []
    for i, f in enumerate(fs):
        p = Params(f=f, **archetype_kw)
        df = simulate_wells(p, n_wells=n_wells, n_experiments=n_experiments,
                            cfg=cfg, seed=seed + i)
        rec = dict(f=f, log10_f=np.log10(f),
                   founders_expected=p.founders_per_well())
        rec.update(summarise(df))
        recs.append(rec)
    return pd.DataFrame(recs)


# --------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    from experiments import ARCHETYPES

    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 40)

    kw = dict(ARCHETYPES["strong"])
    n_exp = int(sys.argv[1]) if len(sys.argv) > 1 else 60

    print("=== stochastic resistance-frequency scan (archetype: strong) ===")
    print(f"12 wells per experiment, {n_exp} experiments per frequency\n")
    df = frequency_scan(kw, n_experiments=n_exp)
    cols = ["f", "founders_expected", "p_zero_founders", "p_ti_capped",
            "SusI_median", "SusI_lo95", "SusI_hi95", "SusI_cv",
            "ti_median", "ti_lo95", "ti_hi95"]
    print(df[cols].to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    df.to_csv("out/stoch_frequency_scan.csv", index=False)
