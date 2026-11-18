"""
Working figures for the SusI simulation study.

These are diagnostic figures for reading the results, not manuscript artwork.
For the manuscript they should be rebuilt in R with the project's existing
conventions (unified cross-experiment palette, graduated line widths with the
thicker curve drawn first).

Colour choices follow the data's job rather than taste:
* the five phage archetypes are ordered by R_phage, which is a magnitude, so
  they take steps from a single-hue sequential blue ramp (light = low R_phage);
* the two-archetype comparisons are an identity contrast, so they take
  categorical slots 1 and 2 (blue, orange), a pair validated for
  colour-vision-deficiency separation (worst adjacent protan dE 24.7).

Usage:  python3 figures.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from model import Params, simulate_pair
from detect import DetectCfg, analyse_pair
from experiments import ARCHETYPES, make, measure, d2_map, d3_frequency, \
    d3_cocktail, d4_threshold_sensitivity

# sequential blue ramp, steps 150/250/400/550/700 (light = low magnitude)
SEQ = ["#b7d3f6", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
CAT = ["#2a78d6", "#eb6834"]
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8c8b85"
SURFACE = "#fcfcfb"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.8,
    "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2,
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "axes.labelsize": 10, "axes.titlesize": 11,
    "legend.frameon": False, "legend.fontsize": 9,
    "grid.color": "#e4e3de", "grid.linewidth": 0.6,
    "font.size": 10, "lines.linewidth": 2.0,
    "svg.fonttype": "none",
})


#: output formats. SVG is included so that figures can be edited in Illustrator,
#: Inkscape or Affinity: text stays as text, so labels can be moved, re-wrapped
#: or resized, and panels can be split apart without re-rendering.
FORMATS = ("png", "svg")


def _save(fig, path):
    """Save one figure in every format in FORMATS. `path` may carry any suffix."""
    import os
    base = os.path.splitext(path)[0]
    written = []
    for ext in FORMATS:
        out = f"{base}.{ext}"
        fig.savefig(out, dpi=170, bbox_inches="tight",
                    format=ext, transparent=False)
        written.append(out)
    plt.close(fig)
    return written[0]


ORDER = ["very blunt", "blunt", "marginal", "strong", "sharp"]


def _tidy(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, axis="y", alpha=0.7)
    ax.set_axisbelow(True)


# --------------------------------------------------------------------------

#: the three panels shown in the manuscript supplement. The full five-way
#: panel is still available by passing names=ORDER. Panels are identified by
#: R_phage, a defined quantity, rather than by an invented category name, so
#: that no figure depends on vocabulary introduced by another figure.
FIG1_PANELS = ["very blunt", "marginal", "sharp"]
FIG1_COLS = [SEQ[0], SEQ[2], SEQ[4]]


def fig1_curves(path="figures/fig1_archetype_curves.png", names=None):
    """Turbidity curves across the range of phage reproductive number, at a
    common resistance frequency, with the detected time points marked."""
    names = names or FIG1_PANELS
    cols = FIG1_COLS if names is FIG1_PANELS else SEQ
    fig, axes = plt.subplots(1, len(names), figsize=(3.0 * len(names), 3.4),
                             sharey=True)
    axes = np.atleast_1d(axes)
    for ax, name, col in zip(axes, names, cols):
        p = make(name, f=1e-6)
        t, odp, odc = simulate_pair(p)
        r = analyse_pair(t, odp, odc, DetectCfg())
        ax.plot(t, odc, color=MUTED, lw=2.0, label="control")
        ax.plot(t, odp, color=col, lw=2.0, label="phage")
        for v, ls, lab in ((r["t0"], ":", "t0"), (r["ti"], "--", "ti"),
                           (r["tc"], "-.", "tc")):
            if np.isfinite(v):
                ax.axvline(v, color=INK2, ls=ls, lw=0.9, alpha=0.65)
                ax.text(v, 1.30, lab, ha="center", va="bottom",
                        fontsize=8, color=INK2)
        ax.set_title(f"$R_{{phage}}$ = {p.r_phage():.2g}", fontsize=10)
        ax.set_xlabel("time (h)")
        ax.set_xlim(0, 36)
        ax.set_ylim(0, 1.45)
        ax.annotate(f"SusI = {r['SusI']:.2f}\n$t_i/t_c$ = {r['ti_over_tc']:.2f}",
                    xy=(0.97, 0.05), xycoords="axes fraction",
                    ha="right", va="bottom", fontsize=9, color=INK)
        _tidy(ax)
    axes[0].set_ylabel("OD")
    axes[0].legend(loc="center right", fontsize=8)
    fig.suptitle("Simulated turbidity across the range of phage reproductive "
                 "number, at a common resistance frequency (f = 1e-6)",
                 fontsize=11, y=1.04)
    fig.tight_layout()
    return _save(fig, path)


def fig2_plane(df, path="figures/fig2_susi_vs_titc.png"):
    """The (ti/tc, SusI) plane. Strongly amplifying phages fall on a tight
    line; marginally amplifying phages fall far off it."""
    from matplotlib.colors import LinearSegmentedColormap, LogNorm
    cmap = LinearSegmentedColormap.from_list("seq_blue", SEQ)
    fig, ax = plt.subplots(figsize=(7.8, 5.0))
    ok = df.dropna(subset=["SusI", "ti_over_tc"])
    sc = ax.scatter(ok.ti_over_tc, ok.SusI, s=58, c=ok.R_phage, cmap=cmap,
                    norm=LogNorm(vmin=ok.R_phage.min(), vmax=ok.R_phage.max()),
                    edgecolor=SURFACE, linewidth=1.2, zorder=3)
    cb = fig.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label("$R_{phage}$", fontsize=10)
    cb.outline.set_visible(False)
    strong = ok[ok.R_phage > 5]
    if len(strong) > 2:
        b, a = np.polyfit(strong.ti_over_tc, strong.SusI, 1)
        xs = np.linspace(strong.ti_over_tc.min(), strong.ti_over_tc.max(), 20)
        ax.plot(xs, a + b * xs, color=INK2, ls="--", lw=1.2, zorder=2)
        r = np.corrcoef(strong.ti_over_tc, strong.SusI)[0, 1]
        ax.annotate(f"strongly amplifying phages\n($R_{{phage}}$ > 5):  r = {r:.3f}",
                    xy=(0.04, 0.95), xycoords="axes fraction", va="top",
                    fontsize=9, color=INK)
    # No whole-panel correlation coefficient is reported. Its value depends on
    # how many marginally amplifying phages the panel happens to contain, so it
    # is a property of the panel and not a result. The claim is qualitative:
    # outside the amplifying regime the relationship does not hold.
    blunt = ok[ok.R_phage < 1.5]
    if not blunt.empty:
        bx, by = blunt.ti_over_tc.mean(), blunt.SusI.mean()
        ax.annotate("$R_{phage}$ < 1: no decline is detected,\n"
                    "so $t_i$ takes the 24 h cap and $t_i/t_c$\n"
                    "is the highest in the panel while\nSusI is near zero",
                    xy=(bx, by), xytext=(bx - 0.78, by + 1.15),
                    fontsize=9, color=INK, ha="center",
                    arrowprops=dict(arrowstyle="->", color=INK2, lw=1.0,
                                    shrinkA=0, shrinkB=6))
    ax.set_xlabel("$t_i/t_c$")
    ax.set_ylabel("SusI")
    ax.set_title("The agreement between SusI and $t_i/t_c$ holds only\n"
                 "where the phage amplifies strongly", loc="left")
    _tidy(ax)
    ax.grid(True, axis="x", alpha=0.7)
    fig.tight_layout()
    return _save(fig, path)


def fig3_frequency(df, path="figures/fig3_susi_vs_frequency.png"):
    """SusI against resistance frequency, for a strongly and a marginally
    amplifying phage: two distinct regimes."""
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.1))
    for ax, ycol, ylab in ((axes[0], "SusI", "SusI"),
                           (axes[1], "duration", "factor value")):
        for (name, col) in zip(("sharp", "blunt"), CAT):
            g = df[df.archetype == name].sort_values("log10_f")
            # label by R_phage, not by an invented category name
            lab = f"$R_{{phage}}$ = {make(name).r_phage():.2g}"
            if ycol == "SusI":
                ax.plot(g.log10_f, g.SusI, "o-", color=col, ms=6,
                        markeredgecolor=SURFACE, markeredgewidth=1.0, label=lab)
            else:
                ax.plot(g.log10_f, g.duration, "o-", color=col, ms=6,
                        markeredgecolor=SURFACE, markeredgewidth=1.0,
                        label=f"{lab}: duration")
                ax.plot(g.log10_f, g.depth, "s--", color=col, ms=5, alpha=0.75,
                        markeredgecolor=SURFACE, markeredgewidth=1.0,
                        label=f"{lab}: depth")
        ax.set_xlabel("$\\log_{10}$ frequency of insensitive cells")
        ax.set_ylabel(ylab)
        ax.legend(loc="best")
        _tidy(ax)
    axes[0].set_title("SusI against resistance frequency", loc="left")
    axes[1].set_title("decomposition into duration and depth", loc="left")
    fig.tight_layout()
    return _save(fig, path)


def fig4_cocktail(df, path="figures/fig4_cocktail_product_rule.png"):
    """Different-receptor cocktails raise SusI because the double-mutant
    frequency is the product of the single-mutant frequencies."""
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    cases = ["single phage", "same-receptor cocktail",
             "different-receptor cocktail"]
    colours = [SEQ[1], SEQ[2], SEQ[4]]
    fs = sorted(df.f_single.unique())
    x = np.arange(len(fs))
    w = 0.26
    for i, (case, col) in enumerate(zip(cases, colours)):
        g = df[df.case == case].set_index("f_single").loc[fs]
        ax.bar(x + (i - 1) * w, g.SusI.values, width=w - 0.02, color=col,
               edgecolor=SURFACE, linewidth=1.5, label=case, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels([f"$10^{{{int(np.log10(v))}}}$" for v in fs])
    ax.set_xlabel("single-mutant frequency")
    ax.set_ylabel("SusI")
    ax.set_title("Receptor independence raises SusI", loc="left")
    ax.legend(loc="upper left")
    _tidy(ax)
    fig.tight_layout()
    return _save(fig, path)


def fig5_sensitivity(summ, path="figures/fig5_threshold_sensitivity.png"):
    """Rank stability and absolute-scale stability under perturbation of every
    detection threshold."""
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0), sharey=True)
    s = summ[summ.variant != "published defaults"].copy()
    s = s.sort_values("spearman_vs_default")
    y = np.arange(len(s))

    axes[0].barh(y, s.spearman_vs_default, color=CAT[0], height=0.62,
                 edgecolor=SURFACE, linewidth=1.2, zorder=3)
    axes[0].axvline(1.0, color=INK2, lw=0.9, ls="--")
    axes[0].set_xlim(0.0, 1.05)
    axes[0].set_xlabel("Spearman rank correlation with published defaults")
    axes[0].set_title("ranking is preserved", loc="left")

    axes[1].barh(y, s.max_abs_change, color=CAT[1], height=0.62,
                 edgecolor=SURFACE, linewidth=1.2, zorder=3)
    axes[1].set_xlabel("largest absolute change in SusI")
    axes[1].set_title("absolute values are not always preserved", loc="left")

    axes[0].set_yticks(y)
    axes[0].set_yticklabels(s.variant, fontsize=9)
    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, axis="x", alpha=0.7)
        ax.set_axisbelow(True)
    fig.tight_layout()
    return _save(fig, path)


def main():
    made = [fig1_curves()]
    d2 = d2_map()
    d2.to_csv("out/d2_map.csv", index=False)
    made.append(fig2_plane(d2))
    d3 = d3_frequency()
    d3.to_csv("out/d3_frequency.csv", index=False)
    made.append(fig3_frequency(d3))
    d3b = d3_cocktail()
    d3b.to_csv("out/d3_cocktail.csv", index=False)
    made.append(fig4_cocktail(d3b))
    summ, wide = d4_threshold_sensitivity()
    summ.to_csv("out/d4_sensitivity_summary.csv", index=False)
    wide.to_csv("out/d4_sensitivity_values.csv")
    made.append(fig5_sensitivity(summ))
    made.append(fig6_partial())
    # fig7 is not part of the manuscript figure set: it answers no reviewer
    # point. It is kept in the module and can be regenerated by calling
    # fig7_stochastic() directly.
    for m in made:
        print("wrote", m)




def fig6_partial(path="figures/fig6_partial_suppression.png"):
    """The case the metric is for: lysis occurs but does not go to completion.

    t0, ti and tc are identical across the family, so ti/tc cannot tell these
    curves apart at all; SusI separates them 2.3-fold.
    """
    from experiments import d5_partial_suppression, d5_curves
    df = d5_partial_suppression()
    residuals = (0.0, 0.30, 0.60)
    t, odc, curves = d5_curves(residuals)
    cols = [SEQ[4], SEQ[2], SEQ[1]]

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6),
                             gridspec_kw={"width_ratios": [1.25, 1]})

    ax = axes[0]
    ax.plot(t, odc, color=MUTED, lw=2.0, label="control")
    for kd, col in zip(residuals, cols):
        row = df[np.isclose(df.residual_fraction, kd)].iloc[0]
        ax.plot(t, curves[kd], color=col, lw=2.0,
                label=f"plateau {row.plateau_pct_of_onset:.0f}% of onset  "
                      f"(SusI {row.SusI:.2f})")
    row = df.iloc[0]
    for v, ls, lab in ((row.t0, ":", "$t_0$"), (row.ti, "--", "$t_i$"),
                       (row.tc, "-.", "$t_c$")):
        ax.axvline(v, color=INK2, ls=ls, lw=0.9, alpha=0.65)
        ax.text(v, 1.30, lab, ha="center", va="bottom", fontsize=9, color=INK2)
    ax.set_xlim(0, 30)
    ax.set_ylim(0, 1.42)
    ax.set_xlabel("time (h)")
    ax.set_ylabel("OD")
    ax.set_title(f"Lysis that does not go to completion\n"
                 f"$t_0$, $t_i$, $t_c$ identical, so $t_i/t_c$ = "
                 f"{row.ti_over_tc:.2f} for every curve", loc="left")
    ax.legend(loc="upper right", bbox_to_anchor=(0.995, 0.78), fontsize=8.5)
    _tidy(ax)

    ax = axes[1]
    x = df.plateau_pct_of_onset
    ax.plot(x, df.SusI, "o-", color=CAT[0], ms=7,
            markeredgecolor=SURFACE, markeredgewidth=1.2, label="SusI")
    ax.plot(x, df.ti_over_tc, "s--", color=MUTED, ms=6,
            markeredgecolor=SURFACE, markeredgewidth=1.2, label="$t_i/t_c$")
    ax.plot(x, df.depth, "^-", color=CAT[1], ms=6, alpha=0.9,
            markeredgecolor=SURFACE, markeredgewidth=1.2, label="depth term")
    ax.set_xlabel("plateau turbidity (% of OD at lysis onset)")
    ax.set_ylabel("index value")
    ax.set_title("$t_i/t_c$ is flat; SusI is not", loc="left")
    ax.legend(loc="center left")
    _tidy(ax)

    fig.tight_layout()
    return _save(fig, path)


def fig7_stochastic(path="figures/fig7_stochastic_wells.png", n_exp=60):
    """Well-level stochastic behaviour: variability, endpoint-capped wells, and
    recovery of the insensitive-cell frequency from replicate structure."""
    import pandas as pd
    from model import Params
    from stochastic import simulate_wells, summarise
    from experiments import ARCHETYPES

    kw = dict(ARCHETYPES["strong"])
    fs = np.logspace(-8.5, -4.5, 9)
    recs, wells = [], {}
    for i, f in enumerate(fs):
        p = Params(f=f, **kw)
        d = simulate_wells(p, n_wells=12, n_experiments=n_exp, seed=100 + i)
        wells[f] = d
        r = dict(f=f, lam=p.founders_per_well())
        r.update(summarise(d))
        recs.append(r)
    df = pd.DataFrame(recs)
    df.to_csv("out/stoch_figure_scan.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.3))

    # A: SusI median with 95% prediction interval across replicate wells
    ax = axes[0]
    x = np.log10(df.f)
    ax.fill_between(x, df.SusI_lo95, df.SusI_hi95, color=CAT[0], alpha=0.18,
                    linewidth=0, label="95% of replicate wells")
    ax.plot(x, df.SusI_median, "o-", color=CAT[0], ms=6,
            markeredgecolor=SURFACE, markeredgewidth=1.2, label="median")
    ax.set_xlabel("$\\log_{10}$ frequency of insensitive cells")
    ax.set_ylabel("SusI")
    ax.set_title("Replicate wells disagree most\nwhere founders are few", loc="left")
    ax.legend(loc="upper right")
    _tidy(ax)

    # B: fraction of wells reaching the endpoint, against the Poisson prediction
    ax = axes[1]
    lam = np.logspace(-1.6, 1.2, 100)
    ax.plot(lam, np.exp(-lam), "-", color=INK2, lw=1.6,
            label="$e^{-\\lambda}$ (no founder cell)")
    ax.plot(df.lam, df.p_ti_capped, "o", color=CAT[1], ms=8,
            markeredgecolor=SURFACE, markeredgewidth=1.2,
            label="simulated wells", zorder=3)
    ax.set_xscale("log")
    ax.set_xlabel("$\\lambda$ = expected founder cells per well")
    ax.set_ylabel("fraction of wells with no $t_i$ by 24 h")
    ax.set_title("Endpoint-capped wells follow\nthe Poisson zero probability", loc="left")
    ax.legend(loc="lower left")
    _tidy(ax)
    ax.grid(True, axis="x", alpha=0.7)

    # C: recovery of f from the capped fraction, on a grid dense enough to
    #    cover the window where the estimate is informative (0.1 < lambda < 5)
    ax = axes[2]
    fs_c = np.logspace(-8.2, -6.4, 10)
    rec_c = []
    for i, f in enumerate(fs_c):
        p = Params(f=f, **kw)
        d = simulate_wells(p, n_wells=12, n_experiments=n_exp, seed=300 + i)
        rec_c.append(dict(f=f, lam=p.founders_per_well(),
                          p_ti_capped=float((d.ti_method == "capped").mean())))
    ok = pd.DataFrame(rec_c)
    ok = ok[(ok.p_ti_capped > 0.02) & (ok.p_ti_capped < 0.95)].copy()
    ok["f_hat"] = -np.log(ok.p_ti_capped) / (Params().N0 * Params().well_volume)
    ok.to_csv("out/stoch_recovery.csv", index=False)
    lims = [min(ok.f.min(), ok.f_hat.min()) * 0.6, max(ok.f.max(), ok.f_hat.max()) * 1.6]
    ax.plot(lims, lims, "--", color=INK2, lw=1.2, label="exact recovery")
    ax.plot(ok.f, ok.f_hat, "o", color=CAT[0], ms=8,
            markeredgecolor=SURFACE, markeredgewidth=1.2,
            label="estimate from well counts", zorder=3)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*lims)
    ax.set_ylim(*lims)
    ax.set_xlabel("true frequency used in the simulation")
    ax.set_ylabel("frequency estimated from wells")
    ax.set_title("Resistance frequency is recoverable\nfrom replicate wells alone", loc="left")
    ax.legend(loc="upper left")
    _tidy(ax)
    ax.grid(True, axis="x", alpha=0.7)

    fig.tight_layout()
    return _save(fig, path), df


if __name__ == "__main__":
    main()
