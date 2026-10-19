"""Build cell_counts.db from cell-count.csv.

Usage: python load_data.py

The CSV has one wide row per sample. This script splits it into the project,
subject, sample and cell_count tables described in schema.sql, then prints a
summary line and a note for anything it skipped or could not interpret.
"""
import csv
import os
import sqlite3
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "cell-count.csv"
DB_PATH = ROOT / "cell_counts.db"
SCHEMA_PATH = ROOT / "schema.sql"

# Metadata columns, with the alternative names the brief uses for some of them
# (sample_id, indication, gender). Every other column that holds whole-number
# counts is treated as a cell population.
METADATA_COLUMNS = {
    "project": ("project",),
    "subject": ("subject", "subject_id"),
    "condition": ("condition", "indication"),
    "age": ("age",),
    "sex": ("sex", "gender"),
    "treatment": ("treatment",),
    "response": ("response",),
    "sample": ("sample", "sample_id"),
    "sample_type": ("sample_type",),
    "time_from_treatment_start": ("time_from_treatment_start",),
}
ID_COLUMNS = ("project", "subject", "sample")
REQUIRED_COLUMNS = ("subject", "sample")

SEX_VALUES = {"m": "M", "male": "M", "f": "F", "female": "F"}
RESPONSE_VALUES = {"yes": "yes", "y": "yes", "true": "yes", "1": "yes", "responder": "yes",
                   "no": "no", "n": "no", "false": "no", "0": "no", "non-responder": "no"}
MISSING_VALUES = {"", "na", "n/a", "nan", "null"}

MAX_IDS_IN_NOTE = 5  # How many ids to list in a note before saying "and N more".


def clean(cell):
    """Strip a cell and turn the usual spellings of "no value" into None."""
    cell = (cell or "").strip()
    if cell.lower() in MISSING_VALUES:
        return None
    return cell


def clean_id(cell):
    """Strip an identifier. Only a blank cell counts as missing, so a subject called "NA" is kept."""
    cell = (cell or "").strip()
    return cell or None


def to_int(text):
    """Parse a whole number written like 123 or 123.0. Anything else, including 1.5, gives None."""
    try:
        number = float(text)
    except (TypeError, ValueError):
        return None
    if not number.is_integer():  # Also False for inf and nan.
        return None
    return int(number)


def display_name(population):
    """b_cell -> B Cell, cd8_t_cell -> CD8 T Cell, nk_cell -> NK Cell."""
    return population.replace("_", " ").title().replace("Cd", "CD").replace("Nk", "NK")


def plural(n, word):
    """"1 sample", "2 samples"."""
    return f"{n} {word}" if n == 1 else f"{n} {word}s"


def shorten(ids):
    """Format a set of ids for a note, e.g. "s1, s2, s3, s4, s5 and 12 more"."""
    ids = sorted(ids)
    if len(ids) <= MAX_IDS_IN_NOTE:
        return ", ".join(ids)
    return ", ".join(ids[:MAX_IDS_IN_NOTE]) + f" and {len(ids) - MAX_IDS_IN_NOTE} more"


def read_csv(csv_path):
    """Read the file into a list of dicts. Returns (header, rows)."""
    with open(csv_path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        header = reader.fieldnames or []
    return header, rows


def looks_like_counts(values):
    """True when most of the cells, blanks included, are whole numbers of zero or more.

    "Most" rather than "all", so that a few bad or blank cells do not throw away
    a whole population column; those cells are reported per sample when the
    counts are read. Blank cells count against the column, so a mostly empty
    extra column is not mistaken for a population.
    """
    values = list(values)
    if not values:
        return False
    good = 0
    for value in values:
        number = to_int(value)
        if number is not None and number >= 0:
            good += 1
    return good > len(values) / 2


def resolve_columns(header, rows):
    """Match the CSV columns to metadata fields and find the population columns.

    Returns (mapping, populations, ignored): mapping is {field: column name},
    populations is the list of count columns, ignored is every other column.
    """
    names = []
    for name in header:
        if name and name.strip() and name not in names:
            names.append(name)
    by_lower = {name.strip().lower(): name for name in names}

    mapping = {}
    for field, accepted in METADATA_COLUMNS.items():
        for name in accepted:
            if name in by_lower:
                mapping[field] = by_lower[name]
                break
    missing = [field for field in REQUIRED_COLUMNS if field not in mapping]
    if missing:
        sys.exit(f"cell-count.csv is missing required columns: {missing}")

    populations, ignored = [], []
    for name in names:
        if name in mapping.values():
            continue
        if looks_like_counts(clean(row.get(name)) for row in rows):
            populations.append(name)
        else:
            ignored.append(name)
    if not populations:
        sys.exit("No cell population columns (whole-number counts) found in cell-count.csv.")
    return mapping, populations, ignored


def normalize_rows(rows, mapping, populations):
    """Rebuild each CSV row as a dict keyed by field name, with cleaned values.

    A metadata field is None when its column is absent or the cell is blank.
    "counts" holds the raw cell for each population.
    """
    out = []
    for raw in rows:
        row = {}
        for field in METADATA_COLUMNS:
            column = mapping.get(field)
            cell = raw.get(column) if column else None
            row[field] = clean_id(cell) if field in ID_COLUMNS else clean(cell)
        row["counts"] = {population: clean(raw.get(population)) for population in populations}
        out.append(row)
    return out


def assign_subject_keys(rows, notes):
    """Add project_id and subject_key to every row.

    Subject ids are normally unique across the whole file. If the same id turns
    up in two projects it has to be two different people, so in that case
    subjects are keyed by project:subject instead.
    """
    projects_of = {}
    for row in rows:
        if row["project"] is not None:
            projects_of.setdefault(row["subject"], set()).add(row["project"])
    scoped = any(len(projects) > 1 for projects in projects_of.values())
    if scoped:
        notes.append("Subject ids repeat across projects, so subjects are keyed by project:subject.")

    for row in rows:
        project = row["project"]
        if project is None:
            # A row with no project takes the one its subject is otherwise filed under.
            known = projects_of.get(row["subject"], set())
            project = list(known)[0] if len(known) == 1 else "unknown"
        row["project_id"] = project
        row["subject_key"] = f"{project}:{row['subject']}" if scoped else row["subject"]


def lookup(field, raw, table, unrecognized):
    """Map a sex or response cell through its lookup table, counting values the table does not know."""
    if raw is None:
        return None
    value = table.get(raw.lower())
    if value is None:
        unrecognized[(field, raw)] += 1
    return value


def parse_int(field, raw, unrecognized):
    """Parse a whole-number cell such as age or timepoint, counting values that are not whole numbers."""
    if raw is None:
        return None
    value = to_int(raw)
    if value is None:
        unrecognized[(field, raw)] += 1
    return value


def add_unrecognized_notes(unrecognized, notes):
    """One note per distinct value that a lookup or parse did not understand."""
    for (field, raw), n in sorted(unrecognized.items()):
        notes.append(f"Unrecognized {field} value {raw!r} stored as NULL ({plural(n, 'row')}).")


def build_subjects(rows, notes):
    """One record per subject.

    A subject's facts come from the first row that has them, so a later row can
    fill in a blank. If two rows give different values the first one is kept
    and the disagreement is reported.
    """
    unrecognized = Counter()
    conflicts = Counter()
    subjects = {}
    for row in rows:
        record = {
            "subject_key": row["subject_key"],
            "subject_id": row["subject"],
            "project_id": row["project_id"],
            "condition": row["condition"].lower() if row["condition"] else None,
            "age": parse_int("age", row["age"], unrecognized),
            "sex": lookup("sex", row["sex"], SEX_VALUES, unrecognized),
            "treatment": row["treatment"].lower() if row["treatment"] else None,
            "response": lookup("response", row["response"], RESPONSE_VALUES, unrecognized),
        }
        key = record["subject_key"]
        if key not in subjects:
            subjects[key] = record
            continue
        kept = subjects[key]
        for field, value in record.items():
            if value is None:
                continue
            if kept[field] is None:
                kept[field] = value
            elif kept[field] != value:
                conflicts[(key, field, kept[field], value)] += 1

    for record in subjects.values():
        if record["condition"] is None:
            record["condition"] = "unknown"
        if record["treatment"] is None:
            record["treatment"] = "none"

    for (key, field, kept, other), n in sorted(conflicts.items()):
        notes.append(f"Subject {key} has more than one {field}: kept {kept!r}, ignored {other!r} "
                     f"({plural(n, 'row')}).")
    add_unrecognized_notes(unrecognized, notes)
    return subjects


def build_samples(rows, notes):
    """One record per sample id, keeping the first row for a repeated id. Returns (samples, rows kept)."""
    unrecognized = Counter()
    samples, kept = [], []
    seen = set()
    duplicates = 0
    no_timepoint = 0
    for row in rows:
        if row["sample"] in seen:
            duplicates += 1
            continue
        seen.add(row["sample"])
        kept.append(row)
        if row["time_from_treatment_start"] is None:
            no_timepoint += 1
        samples.append({
            "sample_id": row["sample"],
            "subject_key": row["subject_key"],
            "sample_type": (row["sample_type"] or "unknown").upper(),
            "time_from_treatment_start": parse_int("timepoint", row["time_from_treatment_start"], unrecognized),
        })
    if duplicates:
        notes.append(f"Skipped {plural(duplicates, 'row')} repeating an earlier sample id.")
    if no_timepoint:
        notes.append(f"No timepoint for {plural(no_timepoint, 'sample')}, stored as NULL.")
    add_unrecognized_notes(unrecognized, notes)
    return samples, kept


def build_counts(rows, populations, notes):
    """Map (sample_id, population) to its count.

    A sample with a blank count, or a count that is not a non-negative whole
    number, is reported and gets no percentages (see the sample_frequency view).
    """
    counts = {}
    blank = {population: set() for population in populations}
    invalid = {population: set() for population in populations}
    example = {}
    for row in rows:
        for population, raw in row["counts"].items():
            if raw is None:
                blank[population].add(row["sample"])
                continue
            value = to_int(raw)
            if value is None or value < 0:
                invalid[population].add(row["sample"])
                example.setdefault(population, raw)
                continue
            counts[(row["sample"], population)] = value
    for population in populations:
        if blank[population]:
            notes.append(f"Left out of the frequency table because the {population} count is blank: "
                         f"{plural(len(blank[population]), 'sample')} ({shorten(blank[population])}).")
        if invalid[population]:
            notes.append(f"Left out of the frequency table because the {population} count is not a "
                         f"non-negative whole number, for example {example[population]!r}: "
                         f"{plural(len(invalid[population]), 'sample')} ({shorten(invalid[population])}).")
    return counts


def write_database(db_path, populations, subjects, samples, counts):
    """Build the database in a temporary file and swap it in once everything is written.

    A failed or interrupted load then cannot leave a half-written cell_counts.db behind.
    """
    temp_path = db_path.with_name(f"{db_path.stem}.{os.getpid()}.tmp")
    conn = sqlite3.connect(temp_path)
    written = False
    try:
        conn.executescript(SCHEMA_PATH.read_text())
        conn.executemany(
            "INSERT INTO population (population, display_name) VALUES (?, ?)",
            [(population, display_name(population)) for population in populations])
        conn.executemany(
            "INSERT INTO project (project_id) VALUES (?)",
            [(project,) for project in sorted({s["project_id"] for s in subjects.values()})])
        conn.executemany(
            "INSERT INTO subject (subject_key, subject_id, project_id, condition, age, sex, treatment, response) "
            "VALUES (:subject_key, :subject_id, :project_id, :condition, :age, :sex, :treatment, :response)",
            list(subjects.values()))
        conn.executemany(
            "INSERT INTO sample (sample_id, subject_key, sample_type, time_from_treatment_start) "
            "VALUES (:sample_id, :subject_key, :sample_type, :time_from_treatment_start)",
            samples)
        conn.executemany(
            "INSERT INTO cell_count (sample_id, population, count) VALUES (?, ?, ?)",
            [(sample, population, n) for (sample, population), n in counts.items()])
        conn.commit()
        written = True
    finally:
        conn.close()
        if written:
            os.replace(temp_path, db_path)
        else:
            temp_path.unlink(missing_ok=True)


def load(csv_path: Path = CSV_PATH, db_path: Path = DB_PATH) -> None:
    """Read the CSV, build every table and write the database."""
    if not csv_path.exists():
        sys.exit(f"Input file not found: {csv_path}")

    header, raw_rows = read_csv(csv_path)
    mapping, populations, ignored = resolve_columns(header, raw_rows)
    rows = normalize_rows(raw_rows, mapping, populations)

    notes = []
    if ignored:
        notes.append(f"Ignored columns that are not whole-number counts: {ignored}.")
    usable = [row for row in rows if row["sample"] and row["subject"]]
    if len(usable) < len(rows):
        notes.append(f"Skipped {plural(len(rows) - len(usable), 'row')} with no sample or subject id.")

    assign_subject_keys(usable, notes)
    subjects = build_subjects(usable, notes)
    samples, kept = build_samples(usable, notes)
    counts = build_counts(kept, populations, notes)
    write_database(db_path, populations, subjects, samples, counts)

    print(f"Loaded {plural(len(samples), 'sample')}, {plural(len(subjects), 'subject')} and "
          f"{plural(len(populations), 'population')} ({', '.join(populations)}) into {db_path.name}.")
    for note in notes:
        print("  Note:", note)


if __name__ == "__main__":
    load()
