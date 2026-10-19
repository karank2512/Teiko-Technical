"""Part 2 against cell_counts.db.

Prints the frequency table and also writes it to outputs/ and to a result_*
table in the database. `make pipeline` runs this after load_data.py.
"""
import sys

import pandas as pd

from analysis import overview
from analysis.db import OUTPUT_DIR, connect, write_table

pd.set_option("display.width", 140)
pd.set_option("display.max_columns", 20)


def section(title: str) -> None:
    """Print a banner between the parts."""
    print(f"\n{'=' * 80}\n{title}\n{'=' * 80}")


def main() -> None:
    """Run Part 2, print it, and save it to outputs/ and the database."""
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

        conn.commit()
    finally:
        conn.close()
    print("Results written to outputs/ and to result_* tables in the database.")


if __name__ == "__main__":
    main()
