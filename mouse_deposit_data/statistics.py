"""
Statistical analysis for the SusI manuscript.

Produces everything the revised Statistical analysis section needs, and
everything Reviewer 2 asked for in minor points 8 and 9:

  * Kaplan-Meier survival estimates
  * log-rank tests between survival curves, with Holm correction for multiple
    comparisons
  * Fisher exact tests on day-7 survival, as a check that does not depend on
    the proportional-hazards assumption
  * Spearman correlation between SusI and day-7 survival, with a bootstrap
    confidence interval
  * Pearson and Spearman correlations between SusI and VI, SupI and ti/tc,
    each with sample size, P value and a bootstrap confidence interval

Inputs
------
survival_per_animal.csv   one row per animal: experiment, group, animal_id,
                          time_days, event (1 = died, 0 = censored at day 7)
susi_vs_survival.csv      SusI and day-7 survival per condition
final_indices.csv         optional; SusI, VI, SupI and ti/tc per condition,
                          produced by the re-analysis scripts

Usage
-----
    python3 statistics.py
    python3 statistics.py --indices ../verification/final_indices.csv

No third-party survival package is required; the log-rank test is implemented
directly so that the computation is visible and auditable.
"""

from __future__ import annotations

import argparse
import itertools
import numpy as np
import pandas as pd
from scipy.stats import chi2, fisher_exact, spearmanr, pearsonr

RNG = np.random.default_rng(20261116)
N_BOOT = 10000


# --------------------------------------------------------------------------
# survival
# --------------------------------------------------------------------------

def kaplan_meier(time, event):
    """Return (times, survival probability) as a step function."""
    time = np.asarray(time, float)
    event = np.asarray(event, int)
    ts = np.sort(np.unique(time[event == 1]))
    s, out = 1.0, [(0.0, 1.0)]
    for t in ts:
        at_risk = (time >= t).sum()
        died = ((time == t) & (event == 1)).sum()
        if at_risk > 0:
            s *= 1 - died / at_risk
        out.append((float(t), float(s)))
    return out


def logrank(t1, e1, t2, e2):
    """Two-sample log-rank test. Returns (chi2, p, observed1, expected1)."""
    t1, e1 = np.asarray(t1, float), np.asarray(e1, int)
    t2, e2 = np.asarray(t2, float), np.asarray(e2, int)
    times = np.sort(np.unique(np.concatenate([t1[e1 == 1], t2[e2 == 1]])))
    O1 = E1 = V = 0.0
    for t in times:
        n1, n2 = (t1 >= t).sum(), (t2 >= t).sum()
        n = n1 + n2
        d1 = ((t1 == t) & (e1 == 1)).sum()
        d2 = ((t2 == t) & (e2 == 1)).sum()
        d = d1 + d2
        if n < 2 or d == 0:
            continue
        O1 += d1
        E1 += d * n1 / n
        V += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
    if V <= 0:
        return np.nan, np.nan, O1, E1
    x2 = (O1 - E1) ** 2 / V
    return float(x2), float(1 - chi2.cdf(x2, 1)), float(O1), float(E1)


def holm(pvals):
    """Holm-Bonferroni adjusted P values, order preserved."""
    p = np.asarray(pvals, float)
    ok = ~np.isnan(p)
    idx = np.flatnonzero(ok)
    order = idx[np.argsort(p[idx])]
    m = len(order)
    adj = np.full_like(p, np.nan)
    running = 0.0
    for rank, i in enumerate(order):
        val = (m - rank) * p[i]
        running = max(running, val)
        adj[i] = min(1.0, running)
    return adj


def survival_tests(df, comparisons):
    rows = []
    for exp, g1, g2 in comparisons:
        a = df[(df.experiment == exp) & (df.group == g1)]
        b = df[(df.experiment == exp) & (df.group == g2)]
        if a.empty or b.empty:
            continue
        x2, p, o1, e1 = logrank(a.time_days, a.event, b.time_days, b.event)
        a_alive, a_dead = int((a.event == 0).sum()), int((a.event == 1).sum())
        b_alive, b_dead = int((b.event == 0).sum()), int((b.event == 1).sum())
        fe = fisher_exact([[a_alive, a_dead], [b_alive, b_dead]])
        rows.append(dict(experiment=exp, group_1=g1, group_2=g2,
                         survival_1=f"{a_alive}/{a_alive + a_dead}",
                         survival_2=f"{b_alive}/{b_alive + b_dead}",
                         chi2=x2, logrank_p=p, fisher_p=float(fe.pvalue)))
    out = pd.DataFrame(rows)
    if not out.empty:
        out["logrank_p_holm"] = holm(out.logrank_p.values)
        out["fisher_p_holm"] = holm(out.fisher_p.values)
    return out


# --------------------------------------------------------------------------
# correlations
# --------------------------------------------------------------------------

def boot_ci(x, y, method="pearson", n_boot=N_BOOT, alpha=0.05):
    """Percentile bootstrap confidence interval for a correlation coefficient."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    stats = []
    for _ in range(n_boot):
        i = RNG.integers(0, n, n)
        if len(np.unique(x[i])) < 2 or len(np.unique(y[i])) < 2:
            continue
        r = (pearsonr(x[i], y[i])[0] if method == "pearson"
             else spearmanr(x[i], y[i]).statistic)
        if np.isfinite(r):
            stats.append(r)
    if len(stats) < 100:
        return np.nan, np.nan
    lo, hi = np.percentile(stats, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(lo), float(hi)


def correlation_row(label, x, y, note=""):
    x, y = np.asarray(x, float), np.asarray(y, float)
    pr, pp = pearsonr(x, y)
    sr = spearmanr(x, y)
    plo, phi = boot_ci(x, y, "pearson")
    slo, shi = boot_ci(x, y, "spearman")
    return dict(comparison=label, n=len(x),
                pearson_r=pr, pearson_p=pp, pearson_CI_low=plo, pearson_CI_high=phi,
                spearman_rho=sr.statistic, spearman_p=sr.pvalue,
                spearman_CI_low=slo, spearman_CI_high=shi, note=note)


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-animal", default="survival_per_animal.csv")
    ap.add_argument("--susi-survival", default="susi_vs_survival.csv")
    ap.add_argument("--indices", default=None,
                    help="final_indices.csv from the re-analysis, optional")
    ap.add_argument("--out", default="statistics_results.xlsx")
    args = ap.parse_args()

    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 40)
    fmt = lambda v: f"{v:.4g}"

    animals = pd.read_csv(args.per_animal)

    # --- Kaplan-Meier ---
    km_rows = []
    for (exp, grp), g in animals.groupby(["experiment", "group"], sort=False):
        for t, s in kaplan_meier(g.time_days, g.event):
            km_rows.append(dict(experiment=exp, group=grp, time_days=t,
                                survival_probability=s))
    km = pd.DataFrame(km_rows)

    # --- log-rank comparisons ---
    # Every phage-treated group against the untreated control in its own
    # experiment, plus the cocktail against each of its constituent phages.
    comparisons = []
    for exp, g in animals.groupby("experiment", sort=False):
        groups = [x for x in g.group.unique()
                  if not x.startswith("Untreated") and not x.startswith("Meropenem")]
        for grp in groups:
            comparisons.append((exp, grp, "Untreated (1x PBS)"))
    comparisons += [("experiment_1", "T6 + PhiWec461", "T6"),
                    ("experiment_1", "T6 + PhiWec461", "PhiWec461")]

    tests = survival_tests(animals, comparisons)
    print("=== log-rank and Fisher exact tests, Holm-corrected ===")
    print(tests.to_string(index=False, float_format=fmt))

    # --- pooled cocktail versus untreated across both experiments ---
    a = animals[animals.group == "T6 + PhiWec461"]
    b = animals[animals.group == "Untreated (1x PBS)"]
    x2, p, _, _ = logrank(a.time_days, a.event, b.time_days, b.event)
    fe = fisher_exact([[int((a.event == 0).sum()), int((a.event == 1).sum())],
                       [int((b.event == 0).sum()), int((b.event == 1).sum())]])
    pooled = pd.DataFrame([dict(
        comparison="T6 + PhiWec461 vs untreated, both experiments pooled",
        survival_1=f"{int((a.event == 0).sum())}/{len(a)}",
        survival_2=f"{int((b.event == 0).sum())}/{len(b)}",
        chi2=x2, logrank_p=p, fisher_p=float(fe.pvalue))])
    print("\n=== pooled comparison ===")
    print(pooled.to_string(index=False, float_format=fmt))

    # --- correlations ---
    cor_rows = []
    ss = pd.read_csv(args.susi_survival)
    cor_rows.append(correlation_row(
        "SusI vs day-7 survival", ss.SusI, ss.day7_survival_percent,
        note="Conditions are not independent: several cocktails share constituent "
             "phages. Reported as a descriptive association, not as a model."))

    if args.indices:
        idx = pd.read_csv(args.indices)
        for strain in idx.strain.unique():
            g = idx[idx.strain == strain]
            for col, lab in (("VI", "VI"), ("SupI", "SupI"), ("titc", "ti/tc")):
                if col in g:
                    cor_rows.append(correlation_row(
                        f"{strain}: SusI vs {lab}", g.SusI, g[col],
                        note="Cocktails share constituent phages; points are not "
                             "fully independent."))
    cors = pd.DataFrame(cor_rows)
    print("\n=== correlations with bootstrap confidence intervals "
          f"({N_BOOT} resamples) ===")
    print(cors.to_string(index=False, float_format=fmt))

    with pd.ExcelWriter(args.out, engine="openpyxl") as xl:
        tests.to_excel(xl, sheet_name="survival_tests", index=False)
        pooled.to_excel(xl, sheet_name="pooled_test", index=False)
        cors.to_excel(xl, sheet_name="correlations", index=False)
        km.to_excel(xl, sheet_name="kaplan_meier", index=False)
    print(f"\nwritten to {args.out}")


if __name__ == "__main__":
    main()
