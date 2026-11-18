# SusI: Sustainability Index for Phage Lysis Kinetics

Code and data accompanying our manuscript introducing the **Sustainability Index (SusI)**, a metric for the duration *and* depth of phage-mediated bacterial growth suppression before turbidity regrowth.

This repository contains three things:

| Folder | What it holds |
|---|---|
| `susi_analysis.Rmd`, `data/` | The R code that calculates SusI, VI, SupI and ti/tc from plate-reader turbidity data, and every raw dataset reported in the manuscript |
| `simulation/` | A Python model of phage–host population dynamics, used to test how SusI behaves under conditions that were not measured experimentally |
| `mouse_deposit_data/` | Survival data from the mouse systemic infection experiments, and the survival analysis |

---

## 1. Indices

| Index | What it captures | Time window | Reference |
|---|---|---|---|
| **SusI** | Duration *and* magnitude of suppression, normalized to the control's growth cycle | t0 → ti, normalized by 0 → tc | This work |
| VI | Initial infection and killing dynamics | 0 → tc | Storms et al., 2020 |
| SupI | Overall growth suppression | **24 h in this study** (Kim et al. use 30 h) | Kim et al., 2024 |
| ti/tc | Simple ratio of suppression duration to control growth-cycle length | ti relative to tc | This work (simplified) |

Three time points are detected automatically from each turbidity curve:

- **t0** — lysis onset
- **ti** — the start of turbidity regrowth
- **tc** — the time at which the phage-free control reaches stationary phase

SusI is the exact product of two terms:

```
SusI = (ti - t0)/tc  x  <C - P>[t0,ti] / <C>[0,tc]
       \__________/     \____________________/
        duration              depth
```

where `C` is the control curve and `P` the phage-treated curve, and `<.>` denotes a time average over the interval shown. Note that the duration term is `(ti - t0)/tc`, **not** `ti/tc`; the two differ by `t0/tc`.

---

## 2. R analysis

### Requirements

R 4.0 or later, plus:

| Package | Used for |
|---|---|
| openxlsx | Reading/writing .xlsx files |
| dplyr | Data wrangling |
| tidyr | Data reshaping |
| ggplot2 | Plotting |
| gridExtra | Combining plots |

```r
install.packages(c("openxlsx", "dplyr", "tidyr", "ggplot2", "gridExtra"))
```

### Data layout

```
data/
├── MG1655/        # 15 datasets: 5 single phages + 10 two-phage cocktails
└── ESBL1064/      # 10 datasets: 4 single phages + 6 two-phage cocktails
```

Each workbook has two sheets: a wide-format OD time course, and a sheet mapping
wells to sample labels. Both strains contain a dataset named `T6`, so the strain
folder is part of the file's identity.

### Quick start

Open `susi_analysis.Rmd` in RStudio, run the function-definition chunks, then:

```r
phage_name <- "T1"
strain     <- "MG1655"
file_path  <- build_file_path(phage_name, strain)

for (method in c("biological_replicates", "individual_wells", "overall_mean")) {
  results <- process_phage_data_v20(
    file_path               = file_path,
    phage_name              = phage_name,
    calculation_method      = method,
    tc_min_time             = 1.0,
    n_forward_points_t0     = 10,
    n_forward_points_ti     = 7,
    threshold_percentage_t0 = 0.6,
    threshold_percentage_ti = 0.6,
    min_increase_threshold  = 0.008,
    t0_sustainability       = 20,
    time_limit_hours        = 24,
    output_file             = paste0(phage_name, "_results.txt"),
    create_plots            = (method == "individual_wells"),
    create_od_curves        = (method == "biological_replicates")
  )
}
```

**The values reported in the manuscript come from `calculation_method = "biological_replicates"`.** The other two methods are provided for comparison; they aggregate the wells differently and give slightly different numbers.

### Parameters

The table below separates the **function default** from the **value used in this study**. They are not the same, and the values used are what reproduce the published numbers.

| Parameter | Function default | **Used in this study** | Meaning |
|---|---|---|---|
| `t0_detection` | `"auto"` | `"auto"` | Detection mode for lysis onset |
| `t0_min_time` | `0.25` | `0.25` | Earliest time (h) t0 may be detected |
| `t0_sustainability` | `5.0` | **`20`** | Hours of sustained decline required to confirm t0 |
| `n_forward_points_t0` | `5` | **`10`** | Forward time points checked when confirming t0 |
| `threshold_percentage_t0` | `0.6` | `0.6` | Fraction of forward points that must be lower to confirm t0 |
| `n_forward_points_ti` | `5` | **`7`** | Forward time points checked when confirming ti |
| `threshold_percentage_ti` | `0.6` | `0.6` | Fraction of forward points that must be higher to confirm ti |
| `min_increase_threshold` | `0.05` | **`0.008`** | Minimum mean OD rise used to detect regrowth (ti) |
| `min_lysis_time` | `2.0` | `2.0` | Minimum interval (h) required between t0 and ti |
| `tc_slope_threshold` | `0.03` | `0.03` | Slope threshold for detecting stationary phase (tc) |
| `tc_min_time` | `1.0` | `1.0`, **`2.0` for ΦWec430** | Earliest time (h) tc may be detected |
| `vi_min_time` | `1.0` | `1.0`, **`2.0` for ΦWec430** | Same, for the VI integration window |
| `time_limit_hours` | `NULL` | **`24`** | Upper time limit (h) applied to **all** calculations, not only ti |

Two details that are easy to miss when re-running the code:

- `time_limit_hours = 24` truncates every calculation at 24 h, not only the detection of ti.
- **VI is integrated from 0 to tc, not to `time_limit_hours`.** Integrating VI to 24 h gives values three to four times too low.

#### Why ΦWec430 needs `tc_min_time = 2.0`

ΦWec430's control culture has a long, nearly flat lag phase. With the default
`tc_min_time = 1.0` the automatic detection locks onto that flat region and
reports tc = 1.0 h, which is far too early and drives VI towards zero. Raising
the search start to 2.0 h gives tc = 7.50 h, in line with every other dataset.
If you see a suspiciously early tc together with a near-zero VI, this is the
cause.

### Output

- a text summary per dataset
- turbidity curves and index plots
- a diagnostic Excel workbook giving per-well t0/ti detection outcomes

---

## 3. Simulation

`simulation/` contains a deterministic model of phage and bacteria growing
together in batch culture, written in Python. It was built to answer questions
the experimental panel cannot: how SusI responds to the frequency of
phage-insensitive cells, how far the index depends on the detection thresholds,
and what happens when lysis stops part-way rather than going to completion.

```
simulation/
├── model.py          # ODE model: two resources, Erlang latent period, debris
├── detect.py         # t0 / ti / tc detection and the index calculations
├── experiments.py    # the parameter sets behind each figure
├── stochastic.py     # per-well Poisson founder-cell model
├── figures.py        # regenerates every figure, as PNG and SVG
├── figures/
└── requirements.txt
```

```bash
cd simulation
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python3 figures.py          # writes figures/ as .png and .svg
```

SVG output keeps text as text, so labels can be edited in a vector editor.

**What the model is for, and what it is not.** It is a tool for asking how the
index behaves, not a fitted description of any particular phage. Burst size and
adsorption rate were not measured for the phages used in this study, so the
model is calibrated to reproduce the shape and timing of the measured control
curves rather than the properties of individual phages. No simulated number
should be read as a prediction about a named phage.

---

## 4. Mouse survival data

```
mouse_deposit_data/
├── mouse_survival_deposit.xlsx   # README sheet + the tables below
├── survival_daily.csv            # alive counts, day 0-7, by experiment and group
├── survival_per_animal.csv       # per-animal time-to-event, for Kaplan-Meier
├── susi_vs_survival.csv          # SusI against 7-day survival
├── build_deposit.py              # how the tables were derived from the source file
└── statistics.py                 # Kaplan-Meier, log-rank, Holm, bootstrap CIs
```

```bash
cd mouse_deposit_data
python3 statistics.py --per-animal --susi-survival --out statistics_results.xlsx
```

Only the two phage treatment experiments reported in the manuscript are
deposited. Individual animals are not identifiable in the source records;
per-animal rows are reconstructed from the daily counts, which is sufficient for
survival analysis.

`statistics.py` implements the log-rank test directly rather than calling an
external survival package, so the calculation can be read and checked.

---

## 5. Citation

Manuscript citation pending. Please cite this repository via its Zenodo DOI:
<https://doi.org/10.5281/zenodo.21890229>

See `CITATION.cff`.

---

## 6. License

Code: MIT License (see `LICENSE`).
Data in `data/` and `mouse_deposit_data/`: CC BY 4.0.

---

## 7. Contact

Please open a GitHub issue, or contact the corresponding author.
