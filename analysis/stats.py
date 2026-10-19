"""Part 3: compare responders with non-responders.

For each population: a two-sided Mann-Whitney U test on the relative
frequencies, Benjamini-Hochberg correction across the populations, Cliff's
delta as the effect size, and Welch's t-test alongside for reference. The
comparison runs on every sample (the primary analysis, as the brief asks), on
per-subject means, per timepoint, and on each subject's change from baseline.
The README explains the reasoning behind these choices.
"""
import sqlite3
import warnings

import numpy as np
import pandas as pd
from scipy import stats

from .db import query

TIME = "time_from_treatment_start"

COHORT_SQL = """
SELECT
    f.sample, f.population, f.percentage, f.count, f.total_count,
    sub.subject_key, sub.subject_id, s.time_from_treatment_start,
    sub.response, sub.sex, sub.project_id
FROM sample_frequency f
JOIN sample  s   ON s.sample_id = f.sample
JOIN subject sub ON sub.subject_key = s.subject_key
WHERE LOWER(sub.condition) = LOWER(?)
  AND LOWER(sub.treatment) = LOWER(?)
  AND UPPER(s.sample_type) = UPPER(?)
  AND sub.response IN ('yes', 'no')
"""

RESULT_COLUMNS = [
    "population", "unit", "n_responders", "n_non_responders",
    "median_responders", "median_non_responders", "median_difference",
    "mean_responders", "mean_non_responders", "cliffs_delta",
    "mannwhitney_u", "p_value", "welch_t", "welch_p_value",
    "p_adjusted_bh", "significant", "direction",
]


def response_cohort(condition: str = "melanoma", treatment: str = "miraclib", sample_type: str = "PBMC",
                    conn: sqlite3.Connection | None = None) -> pd.DataFrame:
    """Per-sample relative frequencies for the cohort, with a "group" column (responder or non-responder)."""
    df = query(COHORT_SQL, (condition, treatment, sample_type), conn=conn)
    df["group"] = np.where(df["response"] == "yes", "responder", "non-responder")
    return df


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """P(a > b) - P(a < b), taken straight from the Mann-Whitney U statistic."""
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    u, _ = stats.mannwhitneyu(a, b, alternative="two-sided")
    return float(2 * u / (len(a) * len(b)) - 1)


def bh_adjust(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values.

    A NaN p-value (a population with no variation, say) stays NaN and does not
    count toward the number of tests.
    """
    p = np.asarray(p, dtype=float)
    adjusted = np.full(p.shape, np.nan)
    tested = np.flatnonzero(np.isfinite(p))
    n = len(tested)
    if n == 0:
        return adjusted
    order = np.argsort(p[tested])
    ranked = p[tested][order] * n / (np.arange(n) + 1)
    # Walk from the largest p-value down, taking a running minimum, so no
    # adjusted p-value is larger than the one ranked above it.
    running_min = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted[tested[order]] = np.clip(running_min, 0, 1)
    return adjusted


def direction_label(delta: float, positive: str, negative: str) -> str:
    """Put the sign of Cliff's delta into words."""
    if pd.isna(delta) or delta == 0:
        return "no difference"
    return positive if delta > 0 else negative


def _empty() -> pd.DataFrame:
    # Same columns, no rows, so callers can concat, filter and write without special cases.
    return pd.DataFrame(columns=RESULT_COLUMNS)


def compare_groups(cohort: pd.DataFrame, unit: str = "sample", alpha: float = 0.05) -> pd.DataFrame:
    """Test every population, responders against non-responders.

    unit="sample" treats every sample as one observation; unit="subject"
    averages each subject's samples first. A population is skipped when either
    group has fewer than two observations.
    """
    cohort = cohort.dropna(subset=["percentage"])  # A sample with zero cells in total has NULL percentages.
    if cohort.empty:
        return _empty()
    df = cohort
    if unit == "subject":
        df = cohort.groupby(["subject_key", "population", "group"], as_index=False)["percentage"].mean()

    rows = []
    for population, sub in df.groupby("population", sort=False):
        responders = sub.loc[sub["group"] == "responder", "percentage"].to_numpy()
        non_responders = sub.loc[sub["group"] == "non-responder", "percentage"].to_numpy()
        if len(responders) < 2 or len(non_responders) < 2:
            continue
        with warnings.catch_warnings():
            # SciPy warns when a group is constant; the p-value is NaN in that case anyway.
            warnings.simplefilter("ignore", RuntimeWarning)
            mw = stats.mannwhitneyu(responders, non_responders, alternative="two-sided")
            welch = stats.ttest_ind(responders, non_responders, equal_var=False)
        rows.append({
            "population": population,
            "unit": unit,
            "n_responders": int(len(responders)),
            "n_non_responders": int(len(non_responders)),
            "median_responders": float(np.median(responders)),
            "median_non_responders": float(np.median(non_responders)),
            "median_difference": float(np.median(responders) - np.median(non_responders)),
            "mean_responders": float(responders.mean()),
            "mean_non_responders": float(non_responders.mean()),
            "cliffs_delta": cliffs_delta(responders, non_responders),
            "mannwhitney_u": float(mw.statistic),
            "p_value": float(mw.pvalue),
            "welch_t": float(welch.statistic),
            "welch_p_value": float(welch.pvalue),
        })

    res = pd.DataFrame(rows)
    if res.empty:
        return _empty()
    res["p_adjusted_bh"] = bh_adjust(res["p_value"].to_numpy())
    res["significant"] = res["p_adjusted_bh"] < alpha
    # The direction follows Cliff's delta, which is what the rank test orders on,
    # so it always agrees with the p-value next to it.
    res["direction"] = [direction_label(delta, "higher in responders", "lower in responders")
                        for delta in res["cliffs_delta"]]
    return res.sort_values("p_adjusted_bh").reset_index(drop=True)


def compare_by_timepoint(cohort: pd.DataFrame, alpha: float = 0.05) -> pd.DataFrame:
    """Run compare_groups at each timepoint separately, so each subject contributes one sample to each test."""
    frames = []
    for timepoint, sub in cohort.groupby(TIME, sort=True):
        res = compare_groups(sub, unit="sample", alpha=alpha)
        if not res.empty:
            res.insert(1, TIME, int(timepoint))
            frames.append(res)
    return pd.concat(frames, ignore_index=True) if frames else _empty()


# Used when the cohort looks fine but, after pairing each subject's samples, a
# group is left with fewer than two subjects.
TOO_FEW_SUBJECTS = "too few subjects have both a time 0 sample and a later one"


def baseline_skip_reason(cohort: pd.DataFrame) -> str | None:
    """Why the change-from-baseline test cannot run on this cohort, or None when it can."""
    if cohort.empty:
        return "no sample matches this cohort"
    if cohort[TIME].notna().sum() == 0:
        return "no sample in this cohort has a timepoint"
    times = cohort[TIME].dropna()
    if not (times == 0).any():
        return "no sample in this cohort is at time 0, so there is no baseline"
    if not (times > 0).any():
        return "no sample in this cohort is later than time 0, so there is nothing to compare baseline with"
    return None


def change_from_baseline(cohort: pd.DataFrame, alpha: float = 0.05) -> pd.DataFrame:
    """Per subject, the change from the time 0 sample to that subject's last sample, compared between groups.

    Baseline is always time 0, as the brief defines it, never just the earliest
    timepoint in the cohort. Subjects without a time 0 sample or without a
    later sample are left out.
    """
    if baseline_skip_reason(cohort) is not None:
        return _empty()
    keys = ["subject_key", "population", "group"]
    baseline = cohort[cohort[TIME] == 0].groupby(keys, as_index=False)["percentage"].mean()
    later = cohort[cohort[TIME] > 0]
    last_time = later.groupby("subject_key")[TIME].transform("max")
    final = (later[later[TIME] == last_time]
             .groupby(keys, as_index=False)
             .agg(percentage=("percentage", "mean"), final_time=(TIME, "first")))
    change = baseline.merge(final, on=keys, suffixes=("_baseline", "_final"))
    change["percentage"] = change["percentage_final"] - change["percentage_baseline"]

    res = compare_groups(change[keys + ["percentage"]], unit="subject", alpha=alpha)
    if res.empty:
        return res
    res["unit"] = "change_from_baseline"
    res.insert(1, "baseline_time", 0)
    res.insert(2, "final_time", int(change["final_time"].max()))
    res["direction"] = [direction_label(delta, "increases more in responders", "decreases more in responders")
                        for delta in res["cliffs_delta"]]
    return res


def run(condition: str = "melanoma", treatment: str = "miraclib", sample_type: str = "PBMC",
        conn: sqlite3.Connection | None = None) -> dict[str, pd.DataFrame]:
    """Every Part 3 table, keyed by name. Used by both the pipeline and the dashboard."""
    cohort = response_cohort(condition, treatment, sample_type, conn=conn)
    pooled = pd.concat(
        [compare_groups(cohort, unit="sample"), compare_groups(cohort, unit="subject")],
        ignore_index=True,
    )
    return {
        "response_cohort": cohort,
        "response_stats": pooled,
        "response_by_timepoint": compare_by_timepoint(cohort),
        "response_trajectory": change_from_baseline(cohort),
    }


def summarize(results: pd.DataFrame, unit: str = "sample") -> str:
    """One sentence about the significant populations, for the terminal and the dashboard."""
    if results.empty:
        return "Not enough data to compare responders and non-responders."
    r = results[results["unit"] == unit]
    if r.empty:
        return "Not enough data to compare responders and non-responders."
    sig = r[r["significant"].astype(bool)]
    if sig.empty:
        tested = r.dropna(subset=["p_adjusted_bh"])
        if tested.empty:
            return "No population could be tested (the frequencies do not vary between samples)."
        closest = tested.sort_values("p_adjusted_bh").iloc[0]
        return ("No population shows a significant difference after FDR correction "
                f"(closest: {closest.population}, FDR p = {closest.p_adjusted_bh:.3f}).")
    is_change = unit == "change_from_baseline"
    label = "median change" if is_change else "median"
    fmt = "{:+.2f}" if is_change else "{:.2f}"
    parts = [
        f"{row.population} {row.direction} "
        f"({label} {fmt.format(row.median_responders)}% vs {fmt.format(row.median_non_responders)}%, "
        f"Cliff's delta {row.cliffs_delta:+.2f}, FDR p = {row.p_adjusted_bh:.3g})"
        for row in sig.itertuples()
    ]
    return "Significant: " + "; ".join(parts) + "."
