-- Schema for cell_counts.db.
--
-- The CSV is one wide row per sample. Here it is split into project -> subject
-- -> sample, with the counts stored long (one row per sample and population),
-- so nothing downstream has to know how many populations there are.

PRAGMA foreign_keys = ON;

DROP VIEW  IF EXISTS sample_frequency;
DROP TABLE IF EXISTS cell_count;
DROP TABLE IF EXISTS population;
DROP TABLE IF EXISTS sample;
DROP TABLE IF EXISTS subject;
DROP TABLE IF EXISTS project;

CREATE TABLE project (
    project_id   TEXT PRIMARY KEY
);

-- One row per patient. Condition, treatment and response are patient-level facts.
CREATE TABLE subject (
    subject_key  TEXT PRIMARY KEY,   -- subject_id, or project:subject_id if ids repeat across projects
    subject_id   TEXT NOT NULL,      -- the id as written in the source file
    project_id   TEXT NOT NULL REFERENCES project(project_id),
    condition    TEXT NOT NULL,      -- melanoma, carcinoma, healthy
    age          INTEGER,
    sex          TEXT CHECK (sex IN ('M', 'F') OR sex IS NULL),
    treatment    TEXT,               -- miraclib, phauximab, none
    response     TEXT CHECK (response IN ('yes', 'no') OR response IS NULL)
);

CREATE TABLE sample (
    sample_id                 TEXT PRIMARY KEY,
    subject_key               TEXT NOT NULL REFERENCES subject(subject_key),
    sample_type               TEXT NOT NULL,   -- PBMC, WB
    time_from_treatment_start INTEGER          -- NULL if the row had no timepoint
);

CREATE TABLE population (
    population   TEXT PRIMARY KEY,
    display_name TEXT NOT NULL
);

CREATE TABLE cell_count (
    sample_id    TEXT NOT NULL REFERENCES sample(sample_id),
    population   TEXT NOT NULL REFERENCES population(population),
    count        INTEGER NOT NULL CHECK (count >= 0),
    PRIMARY KEY (sample_id, population)
);

CREATE INDEX idx_subject_lookup ON subject(condition, treatment, response, sex);
CREATE INDEX idx_sample_lookup  ON sample(subject_key, sample_type, time_from_treatment_start);

-- Part 2 lives here as a view so it cannot drift from the raw counts. A sample
-- is only included once it has a count for every population; otherwise its
-- total would be summed over fewer populations and its percentages would come
-- out too high. The loader prints a note naming any sample that is left out.
CREATE VIEW sample_frequency AS
WITH totals AS (
    SELECT sample_id, SUM(count) AS total_count, COUNT(*) AS populations_counted
    FROM cell_count
    GROUP BY sample_id
)
SELECT
    c.sample_id                                    AS sample,
    t.total_count,
    c.population,
    c.count,
    ROUND(100.0 * c.count / t.total_count, 4)      AS percentage
FROM cell_count c
JOIN totals t USING (sample_id)
WHERE t.populations_counted = (SELECT COUNT(*) FROM population);
