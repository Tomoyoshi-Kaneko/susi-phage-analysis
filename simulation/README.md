# SusI simulation

Mechanistic simulation of phage–host population dynamics, used to characterise
how the Sustainability Index (SusI) behaves under lysis kinetics that are not
represented in the experimental phage panel of the accompanying study.

The simulation generates synthetic turbidity curves whose underlying biology is
known exactly, then analyses them with the **same detection algorithm and the
same thresholds** that were applied to the measured growth curves. This makes
it possible to ask what SusI measures, and when it differs from the simpler
time ratio tᵢ/t꜀, without confounding the answer with experimental noise.

---

## Install

Python 3.10 or newer. Work in a virtual environment so that nothing is
installed system-wide:

```bash
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
```

To remove everything afterwards, run `deactivate` and delete the `.venv`
folder.

## Run

```bash
python3 experiments.py all     # numerical results -> out/*.csv
python3 stochastic.py 60       # well-level stochastic model
python3 figures.py             # figures -> figures/*.png and figures/*.svg
```

Individual experiments can be run as `python3 experiments.py d1` and so on.

## Files

| file | contents |
|---|---|
| `model.py` | the ODE model: state variables, parameters, integration |
| `detect.py` | t₀, tᵢ and t꜀ detection, and the SusI, SupI and VI calculations, with every threshold exposed as an argument |
| `experiments.py` | the numbered experiments, the phage archetypes, and the reproductive number |
| `stochastic.py` | well-level model in which the number of insensitive founder cells per well is drawn from a Poisson distribution |
| `figures.py` | all figures, written as both PNG and SVG |
| `DESIGN.md` | **design record**: equations and the reasoning behind each modelling choice, the full parameter table with sources, the calibration procedure, the results, and an explicit list of what the simulation does **not** show |

Figures are written in PNG and SVG. The SVG keeps text as text rather than
converting it to outlines, so labels can be moved, re-wrapped or resized, and
multi-panel figures can be split apart, in Illustrator, Inkscape or Affinity
without re-rendering. Set `FORMATS` at the top of `figures.py` to change this.

## Which figures appear in the paper

`figures.py` produces both the figures used in the paper and several
diagnostics kept for transparency. The mapping is:

| generated file | used as |
|---|---|
| `fig6_partial_suppression.png` | **Fig. 4A, 4B** |
| `fig3_susi_vs_frequency.png` | **Fig. 4C** |
| `fig2_susi_vs_titc.png` | **Fig. 4D** |
| `fig4_cocktail_product_rule.png` | **Supplementary** |
| `fig5_threshold_sensitivity.png` | **Supplementary** |
| `fig7_stochastic_wells.png` | **Supplementary** |
| `fig1_archetype_curves.png` | diagnostic only, not in the paper |

The published figures were redrawn in R to match the conventions of the other
figures in the paper; the versions here are the ones produced directly by the
analysis.

## Model in brief

Susceptible cells, reversibly tolerant cells, infected cells in ten sequential
stages giving a realistic latent period, heritably insensitive cells, free
phage, and light-scattering debris, growing on a fast and a slow resource.
Burst size and lysis rate are tied to the host's physiological state, because
phage replication needs an actively growing host.

Two quantities control the two factors of SusI:

* the **frequency of insensitive cells** sets how long suppression lasts;
* the phage's **effective reproductive number** in the culture, burst size
  times the probability that a free phage adsorbs before it is lost, sets how
  deep the suppression is. Depth saturates above a reproductive number of about
  five.

`DESIGN.md` records two earlier versions of the model that failed, and why.
Both failures would have produced a wrong answer if they had gone unnoticed,
and they are kept in the record for that reason.

## Changing parameters

Everything is exposed through two dataclasses; nothing is hard-coded.

```python
from model import Params, simulate_pair
from detect import DetectCfg, analyse_pair

p = Params(delta=5e-9, beta=40, tau=0.8, omega=2.0, f=1e-5, moi=0.01)
t, od_phage, od_control = simulate_pair(p)

cfg = DetectCfg(ti_min_rise=0.004, ti_cap=36.0)
r = analyse_pair(t, od_phage, od_control, cfg)
print(r["SusI"], r["duration_term"], r["depth_term"], r["ti_over_tc"])
```

`ARCHETYPES` and `D5_BASE` in `experiments.py` hold the phage parameter sets
used for the figures.

## Note

`detect.py` is a reimplementation of the detection algorithm with every
threshold as an argument, which is what makes the sensitivity analysis
possible. It is not byte-identical to the R implementation used for the
experimental analysis: it uses a central difference for the control growth rate
and excludes the candidate point from the look-ahead window, where the R code
uses a forward difference and includes it. These differences do not affect any
conclusion reported here, but simulated and measured index values should not be
quoted side by side without first aligning the two.

The experimental analysis itself is in `../susi_analysis.Rmd`.
