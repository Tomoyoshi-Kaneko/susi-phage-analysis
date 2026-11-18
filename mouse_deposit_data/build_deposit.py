"""
Build a deposit-grade workbook of the mouse survival data.

Source: the working file supplied by D. Tanaka. Only the two phage treatment
experiments reported in the manuscript are included here; dose-finding and
technique-establishment runs from the same file are pilot work and are not
deposited. This script rewrites
them as one tidy table plus a per-animal table suitable for Kaplan-Meier and
log-rank analysis, and records what each experiment was for.

Counts were transcribed from the "alive" block of each sheet. Per-animal records
are reconstructed from the daily counts: if a group falls from a to b alive
between day d-1 and day d, then a-b animals have time = d and event = 1, and
animals still alive on day 7 are censored at 7. Individual animals are not
identifiable in the source data, which does not matter for survival analysis.
"""

from __future__ import annotations
import numpy as np
import pandas as pd

# experiment id -> (title, purpose, note, {group: [alive counts, day 0..7]})
EXPERIMENTS = {
    "experiment_1": dict(
        title="Phage treatment experiment 1",
        purpose="Single phages versus the T6 + PhiWec461 cocktail",
        note="Reported in the manuscript as the first of the two phage administration experiments.",
        groups={
            "Untreated (1x PBS)": [5, 0, 0, 0, 0, 0, 0, 0],
            "Meropenem 30 mg/kg": [5, 5, 5, 5, 5, 5, 5, 5],
            "T6": [5, 1, 0, 0, 0, 0, 0, 0],
            "PhiWec461": [5, 0, 0, 0, 0, 0, 0, 0],
            "T6 + PhiWec461": [5, 5, 5, 5, 5, 5, 5, 5],
        }),
    "experiment_2": dict(
        title="Phage treatment experiment 2",
        purpose="Five phage conditions spanning a range of SusI",
        note="Reported in the manuscript as the second phage administration experiment.",
        groups={
            "Untreated (1x PBS)": [5, 0, 0, 0, 0, 0, 0, 0],
            "Meropenem 30 mg/kg": [5, 5, 5, 5, 5, 5, 5, 5],
            "T6 + PhiWec461": [5, 5, 5, 5, 5, 5, 5, 5],
            "T6 + PhiWec464": [5, 4, 4, 3, 3, 3, 3, 3],
            "PhiWec430 + PhiWec464": [5, 3, 3, 3, 3, 3, 3, 3],
            "PhiWec430 + PhiWec461": [5, 4, 4, 4, 4, 4, 4, 4],
            "PhiWec430": [5, 2, 2, 2, 2, 2, 2, 2],
        }),
}

# SusI (MOI = 1, biological_replicates method) against day-7 survival
SUSI_SURVIVAL = [
    ("T6 + PhiWec461", 2.776, 0.790, 100, "experiment_2"),
    ("T6 + PhiWec464", 2.187, 0.066, 60, "experiment_2"),
    ("PhiWec430 + PhiWec464", 1.504, 0.306, 60, "experiment_2"),
    ("PhiWec430 + PhiWec461", 1.272, 0.237, 80, "experiment_2"),
    ("PhiWec430", 0.101, 0.028, 40, "experiment_2"),
    ("PhiWec461", 0.291, 0.032, 0, "experiment_1"),
    ("T6", 0.351, 0.008, 0, "experiment_1"),
]


def long_table():
    rows = []
    for exp, meta in EXPERIMENTS.items():
        for group, counts in meta["groups"].items():
            n0 = counts[0]
            for day, alive in enumerate(counts):
                rows.append(dict(experiment=exp, group=group, day=day,
                                 n_alive=alive, n_total=n0,
                                 survival_percent=round(100 * alive / n0, 1)))
    return pd.DataFrame(rows)


def per_animal_table():
    rows = []
    for exp, meta in EXPERIMENTS.items():
        for group, counts in meta["groups"].items():
            n0 = counts[0]
            aid = 0
            for day in range(1, len(counts)):
                deaths = counts[day - 1] - counts[day]
                for _ in range(deaths):
                    aid += 1
                    rows.append(dict(experiment=exp, group=group,
                                     animal_id=f"{aid:02d}", time_days=day, event=1))
            for _ in range(counts[-1]):
                aid += 1
                rows.append(dict(experiment=exp, group=group,
                                 animal_id=f"{aid:02d}", time_days=7, event=0))
    return pd.DataFrame(rows)


def readme_table():
    rows = []
    for exp, meta in EXPERIMENTS.items():
        n_per = sorted({c[0] for c in meta["groups"].values()})
        rows.append(dict(experiment=exp, title=meta["title"], purpose=meta["purpose"],
                         groups="; ".join(meta["groups"]),
                         n_per_group=", ".join(str(n) for n in n_per),
                         follow_up_days=7, note=meta["note"]))
    return pd.DataFrame(rows)


def susi_table():
    return pd.DataFrame(SUSI_SURVIVAL,
                        columns=["condition", "SusI", "SusI_SE",
                                 "day7_survival_percent", "source_experiment"])


def main():
    long_df = long_table()
    animal_df = per_animal_table()
    with pd.ExcelWriter("mouse_survival_deposit.xlsx", engine="openpyxl") as xl:
        readme_table().to_excel(xl, sheet_name="README", index=False)
        long_df.to_excel(xl, sheet_name="survival_daily", index=False)
        animal_df.to_excel(xl, sheet_name="survival_per_animal", index=False)
        susi_table().to_excel(xl, sheet_name="susi_vs_survival", index=False)
    long_df.to_csv("survival_daily.csv", index=False)
    animal_df.to_csv("survival_per_animal.csv", index=False)
    susi_table().to_csv("susi_vs_survival.csv", index=False)
    print(f"survival_daily: {len(long_df)} rows")
    print(f"survival_per_animal: {len(animal_df)} rows "
          f"({animal_df.event.sum()} deaths, {(animal_df.event == 0).sum()} censored)")
    print("\nday-7 survival by experiment and group:")
    d7 = long_df[long_df.day == 7][["experiment", "group", "n_alive", "n_total",
                                     "survival_percent"]]
    print(d7.to_string(index=False))


if __name__ == "__main__":
    main()
