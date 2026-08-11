# SusI: Sustainability Index for Phage Lysis Kinetics

R code to calculate and visualize four metrics of phage-mediated bacterial growth suppression from 96-well plate turbidity (OD) time-course data:

- **SusI** (Sustainability Index) — duration *and* depth of lysis suppression, from lysis onset (t0) to resistance emergence (ti)
- **VI** (Virulence Index) — early infection dynamics, following Storms et al. (2020)
- **SupI** (Suppression Index) — overall growth suppression over a fixed observation window, following Kim et al. (2024)
- **ti/tc** — a simple time-ratio metric, provided for comparison

This repository accompanies our manuscript introducing SusI. A sample dataset (`data/T1_rawdata.xlsx`) is included so the code can be run immediately.

## What each index measures

| Index | What it captures | Time window used | Reference |
|---|---|---|---|
| SusI | Duration *and* magnitude of suppression, normalized to the control's growth cycle | t0 → ti | This work |
| VI | Initial infection/killing dynamics | Early fixed window | Storms et al., 2020 |
| SupI | Overall growth suppression | Fixed 30 h window | Kim et al., 2024 |
| ti/tc | Simple ratio of suppression duration to control growth-cycle length | t0 → ti (relative to tc) | This work (simplified) |

## Repository contents

```
.
├── susi_analysis.Rmd    # All calculation + plotting functions, plus quick-start examples
├── data/
│   └── T1_rawdata.xlsx  # Example dataset: phage T1 vs. E. coli MG1655
├── README.md
└── LICENSE
```

## Requirements

- R ≥ 4.0
- R packages:

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

## Quick start

1. Clone or download this repository.
2. Open `susi_analysis.Rmd` in RStudio.
3. Run every chunk from top to bottom **once** — this loads all calculation and plotting functions into your R session (it does not analyze anything by itself).
4. Scroll to the **"Quick start: process_phage_data_v20()"** chunk near the end. It is already pointed at the bundled `data/T1_rawdata.xlsx` sample and will run out of the box.
5. To analyze your own data, either rename your file to match the pattern `build_file_path()` expects (`<phage_name>_rawdata.xlsx` inside `data/`), or just set `file_path` directly to your file.

## Input data format

Each phage/dataset is one `.xlsx` file with two sheets:

**Sheet `R`** — wide-format OD time course. First column is `Time` (seconds); every other column is one well.

| Time | A1 | A2 | A3 | ... |
|---|---|---|---|---|
| 0 | 0.236 | 0.237 | 0.236 | ... |
| 900 | 0.240 | 0.237 | 0.236 | ... |
| 1800 | 0.282 | 0.276 | 0.272 | ... |

**Sheet `name`** — maps each well to a sample/condition label. In the bundled example, wells are labeled `Ct` (phage-free control), `10^0`...`10^-5` (MOI conditions), and `free` (phage-only, no bacteria); 12 wells per condition (3 biological replicates × 4 technical replicates each).

| num | sample |
|---|---|
| A1 | Ct |
| A2 | Ct |
| ... | ... |
| B1 | 10^0 |

## Main functions

| Function | Purpose |
|---|---|
| `process_phage_data_v20()` | Full pipeline for one phage/file: calculates SusI, SupI, VI, and ti/tc, and (optionally) generates all five figures |
| `diagnose_all_mois_enhanced()` | Runs t0/ti detection across every MOI condition in a file, classifies each well as success/failure, and writes a diagnostic Excel report + plots — useful for tuning detection thresholds before a full run |
| `build_file_path()` | Convenience helper: builds a file path from a phage name (edit `base_dir`/`suffix` to match your own file layout) |

Both `process_phage_data_v20()` and `diagnose_all_mois_enhanced()` support three calculation methods, set via `calculation_method`:

- `"biological_replicates"` — mean curve per biological replicate, then averaged (default; matches the paper)
- `"individual_wells"` — every well analyzed independently, no averaging
- `"overall_mean"` — single mean curve across all wells of a condition

## Key parameters

`process_phage_data_v20()` and `diagnose_all_mois_enhanced()` share the parameters below, which control t0/ti/tc detection and therefore the calculated index values (these are the values reported for the paper; see the "Quick start" chunks in `susi_analysis.Rmd` for the exact settings used for each dataset).

| Parameter | Default | Meaning |
|---|---|---|
| `t0_detection` | `"auto"` | Detection mode for lysis onset (t0) |
| `t0_min_time` | `0.25` | Earliest time (h) t0 is allowed to be detected |
| `t0_sustainability` | `5.0` | How many hours of sustained decline are required to confirm t0 |
| `n_forward_points_t0` | `5` | Number of forward time points checked when confirming t0 |
| `threshold_percentage_t0` | `0.6` | Fraction of forward points that must be lower to confirm t0 |
| `n_forward_points_ti` | `5` | Number of forward time points checked when confirming ti (resistance emergence) |
| `threshold_percentage_ti` | `0.6` | Fraction of forward points that must be higher to confirm ti |
| `min_increase_threshold` | `0.05` | Minimum OD increase used to detect the start of regrowth (ti) |
| `min_lysis_time` | `2.0` | Minimum duration (h) required between t0 and ti |
| `tc_slope_threshold` | `0.03` | Slope threshold used to detect when the control reaches stationary phase (tc) |
| `tc_min_time` | `1.0` | Earliest time (h) tc is allowed to be detected — raise this for a long flat lag phase (see note below) |
| `vi_min_time` | `1.0` | Same idea as `tc_min_time`, but for the VI integration window |
| `time_limit_hours` | `NULL` | Upper time limit (h) for all calculations; `NULL` = use the full time course |

Plot-appearance-only arguments (`title_size`, `axis_text_size`, `legend_text_size`, `show_legend`, `show_title`, the various `*_ylim` ranges, etc.) change how figures look but do not affect any calculated value. For the complete, always-up-to-date argument list, run this in R after loading the functions:

```r
args(process_phage_data_v20)
formals(process_phage_data_v20)   # same, but shows default values too
```

## A note on parameter tuning (tc_min_time / vi_min_time)

If a control (phage-free) culture has a long, nearly flat lag phase, the automatic tc/VI-window detection can lock onto that flat region and report a falsely early tc (and a VI close to 0). If you see this, increase `tc_min_time` (and `vi_min_time` for VI) so the search window starts later — e.g. `tc_min_time = 2.0`.

## Output

- `process_phage_data_v20()` appends a text summary to `output_file`, and (if `create_plots`/`create_od_curves` are `TRUE`) prints the growth-curve and three-level SusI/SupI/VI/ti-tc plots.
- `diagnose_all_mois_enhanced()` writes an Excel workbook (`output_excel`) with:

| Sheet | Contents |
|---|---|
| Summary | One row per MOI condition, with detection/success counts |
| Parameters | The detection parameters used for that run |
| `<MOI>` (one per condition) | Per-well diagnostic detail |
| Non_Success | Only the wells where t0 detection did not succeed |

Diagnostic plot colors: **green** = SusI calculated successfully, **red** = t0 detection failed, **orange** = ti detection failed, **purple** = other calculation failure, **dark purple** = no results generated.

## Citation

If you use this code, please cite the associated manuscript (citation to be added on publication — see `CITATION.cff`, and check back here for a preprint/DOI link). A citable DOI for this repository itself is also planned via Zenodo.

## License

Released under the MIT License — see [LICENSE](LICENSE).

## Contact

Questions or issues: please open a GitHub issue, or contact the corresponding author of the associated manuscript (Tsuneda Laboratory, Graduate School of Advanced Science and Engineering, Waseda University).
