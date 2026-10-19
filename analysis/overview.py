"""Part 2: relative frequency of each population in each sample."""
import sqlite3

import pandas as pd

from .db import query


def frequency_table(conn: sqlite3.Connection | None = None) -> pd.DataFrame:
    """The Part 2 summary table: sample, total_count, population, count, percentage.

    The calculation itself is the sample_frequency view in schema.sql.
    """
    return query(
        "SELECT sample, total_count, population, count, percentage "
        "FROM sample_frequency ORDER BY sample, population",
        conn=conn,
    )
