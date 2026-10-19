"""Parts 2 to 4 against cell_counts.db.

Prints each result and also writes it to outputs/*.csv and to a result_* table
in the database so the dashboard can read it. `make pipeline` runs this after
load_data.py.
"""
import sys

import pandas as pd

from analysis import overview, stats, subsets
from analysis.db import OUTPUT_DIR, connect, write_table

pd.set_option("display.width", 140)
pd.set_option("display.max_columns", 20)

SHOW = ["population", "unit", "n_responders", "n_non_responders", "median_responders",
        "median_non_responders", "cliffs_delta", "p_value", "p_adjusted_bh", "significant"]
SHOW_TRAJ = ["population", "n_responders", "n_non_responders", "median_responders",
             "median_non_responders", "cliffs_delta", "p_value", "p_adjusted_bh", "significant"]


def section(title: str) -> None:
    """Print a banner between the parts."""
    print(f"\n{'=' * 80}\n{title}\n{'=' * 80}")


def fmt(value: float) -> str:
    """Four significant figures, for the printed tables."""
    return f"{value:.4g}"


def main() -> None:
    """Run Parts 2 to 4, print them, and save them to outputs/ and the database."""
    OUTPUT_DIR.mkdir(exist_ok=True)
    try:
        conn = connect()
    except FileNotFoundError as e:
        sys.exit(str(e))
    try:
        section("Part 2: relative frequency of each population per sample")
        freq = overview.frequency_table(conn)
        write_table(freq, "result_frequency", conn)
        freq.to_csv(OUTPUT_DIR / "part2_frequency_summary.csv", index=False)
        print(freq.head(10).to_string(index=False))
        print(f"... {len(freq)} rows ({freq['sample'].nunique()} samples x "
              f"{freq['population'].nunique()} populations)")

        section("Part 3: responders vs non-responders (melanoma, miraclib, PBMC)")
        part3 = stats.run(conn=conn)
        for name, df in part3.items():
            write_table(df, f"result_{name}", conn)
            df.to_csv(OUTPUT_DIR / f"part3_{name}.csv", index=False)
        cohort = part3["response_cohort"]
        results = part3["response_stats"]
        if results.empty:
            print("Not enough responders and non-responders in this cohort to run the comparison.")
        else:
            print(results[SHOW].to_string(index=False, float_format=fmt))
        print("\nAll timepoints pooled, sample level:  " + stats.summarize(results, "sample"))
        print("All timepoints pooled, subject level: " + stats.summarize(results, "subject"))
        traj = part3["response_trajectory"]
        if traj.empty:
            reason = stats.baseline_skip_reason(cohort) or stats.TOO_FEW_SUBJECTS
            print(f"\nChange from baseline not tested: {reason}.")
        else:
            print(f"\nChange from time 0 to each subject's last sample (up to day {traj['final_time'].iloc[0]}):")
            print(traj[SHOW_TRAJ].to_string(index=False, float_format=fmt))
            print("\n" + stats.summarize(traj, "change_from_baseline"))

        section("Part 4: melanoma PBMC baseline samples, miraclib-treated")
        tables = subsets.run(conn=conn)
        for name, df in tables.items():
            write_table(df, f"result_{name}", conn)
            df.to_csv(OUTPUT_DIR / f"part4_{name}.csv", index=False)
        base = tables["baseline_samples"]
        print(f"{len(base)} samples from {base['subject_key'].nunique()} subjects\n")
        for name in ("baseline_by_project", "baseline_by_response", "baseline_by_sex"):
            print(tables[name].to_string(index=False), "\n")

        conn.commit()
    finally:
        conn.close()
    print("Results written to outputs/ and to result_* tables in the database.")


if __name__ == "__main__":
    main()
