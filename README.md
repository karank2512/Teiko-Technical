# Loblaw Bio cell count analysis

Take-home for the Loblaw Bio trial data: load `cell-count.csv` into SQLite, answer Parts 2 to 4 of the brief, and put the results on a small Streamlit dashboard.

**Dashboard:** http://localhost:8501 after `make pipeline` and then `make dashboard` (in a Codespace, open the forwarded port 8501).

## Running it

Python 3.10 or newer. Everything runs from the repo root.

```bash
make setup       # python3 -m venv .venv, then pip install -r requirements.txt into it
make pipeline    # python load_data.py && python run_analysis.py
make dashboard   # streamlit run dashboard.py
```

The Makefile runs everything with `.venv/bin/python`, which `make setup` creates (and which `make pipeline` or `make dashboard` will create for you if it is missing), because the Pythons that Homebrew and apt install refuse to `pip install` outside a virtualenv. Pass `PYTHON=...` to run the scripts with a different interpreter. Tested with Python 3.13, pandas 3.0, numpy 2.5, scipy 1.18, streamlit 1.66 and plotly 7.1; `requirements.txt` accepts anything from pandas 2.2 and streamlit 1.51 up to, but not including, the next major versions.

`make pipeline` prints the Part 2 table, the Part 3 statistics and the Part 4 counts, and writes the same things to `outputs/*.csv` and to `result_*` tables in `cell_counts.db`. The database and `outputs/` are build products and are not committed. If the dashboard is started before the pipeline has run (which is what happens on a hosted deploy) it runs the pipeline itself first.

## Part 1: the database

The CSV is one wide row per sample, which is fine for a spreadsheet but awkward once the questions are about subjects. I split it into `project`, `subject` and `sample` tables and store the counts long in `cell_count` (one row per sample and population) with a small `population` lookup. Subject-level facts (condition, age, sex, treatment, response) live once on the subject row. The Part 2 summary is a view, `sample_frequency`, computed from the counts on read so it cannot drift from the raw data. Foreign keys are on and the columns the queries filter by are indexed.

Keying subjects by id alone would merge two patients if two projects ever reused ids, so the loader checks for that and switches to `project:subject` keys if it happens.

The loader needs a `sample` column and a `subject` column; the other metadata columns are optional, and the population columns are whatever other columns hold whole-number counts, so it is not tied to this particular file. It accepts the column names the brief uses (`sample_id`, `indication`, `gender`) as well as the ones in the file, normalizes the usual spellings of sex and response, and reads `NA`-style cells as blank (in the id columns only an empty cell counts as missing, so a subject called `NA` is kept). Anything it skips or cannot interpret is printed as a note: unrecognized values (including timepoints that are not whole numbers) are stored as NULL, a sample with a blank or invalid count is kept but left out of the frequency table rather than given percentages over a partial total, a repeated sample id keeps its first row, and a subject whose rows disagree keeps the first value, with the note naming the subject and the field. If a cohort is too small to test, the analysis returns empty tables and says so instead of falling over.

## Part 2: frequencies

Per sample, sum the populations and express each as a percentage of the total. One row per sample and population with `sample, total_count, population, count, percentage`. It is the `sample_frequency` view, `outputs/part2_frequency_summary.csv`, and the first dashboard tab (filterable, downloadable).

## Part 3: responders vs non-responders

Cohort: melanoma, miraclib, PBMC only, response recorded.

For each population I compare relative frequency between responders and non-responders with a two-sided Mann-Whitney U test (frequencies are bounded and a bit skewed, so I would rather not lean on normality), with Welch's t-test reported alongside, Benjamini-Hochberg correction across the populations, and Cliff's delta plus the median difference as effect sizes. The direction reported for each population follows the sign of Cliff's delta, so it always agrees with the test. The brief asks for the comparison on the per-sample table so that is the primary analysis, but every subject contributes three timepoints, so I repeat it on per-subject means as a check, and I also test each subject's change from its time 0 sample to its last sample (day 14 for every subject in this file), because a drug effect should appear over time rather than at time 0. A result counts as significant at FDR-adjusted p < 0.05. Two details worth knowing: the tests use the view's percentages, which are rounded to four decimals, and a subject without both a time 0 sample and a later one is left out of the change-from-baseline test only.

What the data says, for the file provided:

- Pooled across timepoints, nothing survives FDR correction. CD4 T cells come closest (higher in responders, FDR p = 0.067).
- At baseline the two groups are indistinguishable on every population, which is what you would hope for.
- Over treatment they separate. The B cell fraction drops in responders and stays flat in non-responders; the change from baseline is significant (FDR p = 0.031, Cliff's delta -0.12). CD4 T cells trend the other way but do not reach significance.
- Effects are small. A falling B cell fraction looks like a candidate early marker of response, worth checking prospectively, not something to act on yet.

The second dashboard tab has the boxplots per population, both statistics tables, the median trajectories by timepoint, and lets you change the cohort filters and re-run.

## Part 4: baseline subset

Melanoma PBMC samples at `time_from_treatment_start = 0` from miraclib patients, then samples per project and distinct subjects by response and by sex. `outputs/part4_*.csv` and the third tab, where the filters can be changed.

## Hosting

The dashboard is a local server on purpose, so it always reflects the database the pipeline just built. It also runs as-is on Streamlit Community Cloud (`dashboard.py`, `requirements.txt`, no other config).
