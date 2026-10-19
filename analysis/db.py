"""Shared database helpers: paths, connections and small query wrappers."""
import sqlite3
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DB_PATH = ROOT / "cell_counts.db"
OUTPUT_DIR = ROOT / "outputs"


def connect(db_path: Path = DB_PATH) -> sqlite3.Connection:
    """Open the database with foreign keys switched on."""
    if not db_path.exists():
        raise FileNotFoundError(f"{db_path.name} not found, run `python load_data.py` (or `make pipeline`) first")
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def query(sql: str, params: tuple = (), conn: sqlite3.Connection | None = None) -> pd.DataFrame:
    """Run a query and return the rows as a DataFrame. Opens and closes a connection if none is given."""
    if conn is not None:
        return pd.read_sql_query(sql, conn, params=params)
    conn = connect()
    try:
        return pd.read_sql_query(sql, conn, params=params)
    finally:
        conn.close()


def write_table(df: pd.DataFrame, name: str, conn: sqlite3.Connection) -> None:
    """Replace a table in the database with the contents of a DataFrame."""
    df.to_sql(name, conn, if_exists="replace", index=False)


def populations(conn: sqlite3.Connection | None = None) -> list[str]:
    """Population names in the order they appear in the CSV."""
    return query("SELECT population FROM population ORDER BY rowid", conn=conn)["population"].tolist()
