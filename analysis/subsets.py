"""Part 4: baseline melanoma PBMC samples from miraclib patients, and a few counts over them."""
import sqlite3

import pandas as pd

from .db import query

BASELINE_SQL = """
SELECT
    s.sample_id AS sample, sub.subject_key, sub.subject_id, sub.project_id AS project,
    sub.condition, sub.age, sub.sex, sub.treatment, sub.response, s.sample_type,
    s.time_from_treatment_start
FROM sample s
JOIN subject sub ON sub.subject_key = s.subject_key
WHERE LOWER(sub.condition) = LOWER(?)
  AND LOWER(sub.treatment) = LOWER(?)
  AND UPPER(s.sample_type) = UPPER(?)
  AND s.time_from_treatment_start = ?
ORDER BY s.sample_id
"""


def baseline_samples(condition: str = "melanoma", treatment: str = "miraclib", sample_type: str = "PBMC",
                     timepoint: int = 0, conn: sqlite3.Connection | None = None) -> pd.DataFrame:
    """Every sample matching the condition, treatment, sample type and timepoint, with its subject's details."""
    return query(BASELINE_SQL, (condition, treatment, sample_type, timepoint), conn=conn)


def samples_per_project(baseline: pd.DataFrame) -> pd.DataFrame:
    """Number of samples from each project."""
    return (baseline.groupby("project").size().rename("samples").reset_index()
            .sort_values("project").reset_index(drop=True))


def subjects_by_response(baseline: pd.DataFrame) -> pd.DataFrame:
    """Distinct subjects per response. Subjects with no recorded response get their own row."""
    out = (baseline.groupby("response", dropna=False)["subject_key"].nunique()
           .rename("subjects").reset_index())
    out["response"] = out["response"].map({"yes": "responder", "no": "non-responder"}).fillna("not recorded")
    order = {"responder": 0, "non-responder": 1, "not recorded": 2}
    return out.sort_values("response", key=lambda s: s.map(order)).reset_index(drop=True)


def subjects_by_sex(baseline: pd.DataFrame) -> pd.DataFrame:
    """Distinct subjects per sex. Subjects with no recorded sex count as "unknown"."""
    out = baseline.groupby("sex", dropna=False)["subject_key"].nunique().rename("subjects").reset_index()
    out["sex"] = out["sex"].map({"M": "male", "F": "female"}).fillna("unknown")
    order = {"male": 0, "female": 1, "unknown": 2}
    return out.sort_values("sex", key=lambda s: s.map(order)).reset_index(drop=True)


def run(condition: str = "melanoma", treatment: str = "miraclib", sample_type: str = "PBMC",
        timepoint: int = 0, conn: sqlite3.Connection | None = None) -> dict[str, pd.DataFrame]:
    """Every Part 4 table, keyed by name."""
    baseline = baseline_samples(condition, treatment, sample_type, timepoint, conn=conn)
    return {
        "baseline_samples": baseline,
        "baseline_by_project": samples_per_project(baseline),
        "baseline_by_response": subjects_by_response(baseline),
        "baseline_by_sex": subjects_by_sex(baseline),
    }
