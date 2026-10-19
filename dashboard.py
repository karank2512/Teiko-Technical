"""Streamlit dashboard. Everything on the page comes from cell_counts.db.

The default views read the result_* tables the pipeline wrote; changing a
filter re-queries the database with the same analysis functions.
"""
import subprocess
import sys

import pandas as pd
import plotly.express as px
import streamlit as st

from analysis import db, overview, stats, subsets
from analysis.db import DB_PATH, query

st.set_page_config(page_title="Loblaw Bio cell counts", layout="wide")

TIME = "time_from_treatment_start"
COLORS = {"responder": "#2a78d6", "non-responder": "#eb6834"}
BLUE = "#2a78d6"
HIGHLIGHT = "background-color: rgba(42,120,214,0.12)"
LAYOUT = dict(
    template="plotly_white",
    font=dict(family="Inter, system-ui, sans-serif", size=13),
    margin=dict(l=40, r=20, t=50, b=40),
    legend=dict(orientation="h", yanchor="top", y=-0.18, xanchor="left", x=0),
)

# Columns of the Part 3 result tables that hold numbers. SQLite stores NaN as
# NULL, which comes back as None, so these are converted before formatting.
NUMERIC_RESULT_COLUMNS = [
    "n_responders", "n_non_responders", "median_responders", "median_non_responders",
    "median_difference", "mean_responders", "mean_non_responders", "cliffs_delta",
    "mannwhitney_u", "p_value", "welch_t", "welch_p_value", "p_adjusted_bh",
]
STATS_FORMAT = {
    "median_responders": "{:.2f}", "median_non_responders": "{:.2f}", "median_difference": "{:+.2f}",
    "cliffs_delta": "{:+.3f}", "p_value": "{:.3g}", "p_adjusted_bh": "{:.3g}", "welch_p_value": "{:.3g}",
}
TRAJECTORY_FORMAT = {
    "median_change_responders": "{:+.2f}", "median_change_non_responders": "{:+.2f}",
    "cliffs_delta": "{:+.3f}", "p_value": "{:.3g}", "p_adjusted_bh": "{:.3g}",
}


@st.cache_data(show_spinner=False)
def load_table(sql, params=()):
    """Run a query against the database, cached per statement."""
    return query(sql, params)


@st.cache_data(show_spinner=False)
def frequency():
    """The Part 2 table."""
    return overview.frequency_table()


@st.cache_data(show_spinner="Running the comparison for this cohort...")
def live_part3(condition, treatment, sample_type):
    """Part 3 recomputed for a cohort other than the brief's."""
    return stats.run(condition, treatment, sample_type)


@st.cache_data(show_spinner=False)
def live_part4(condition, treatment, sample_type, timepoint):
    """Part 4 recomputed for filters other than the brief's."""
    return subsets.baseline_samples(condition, treatment, sample_type, timepoint)


@st.cache_resource
def last_seen():
    """One dict shared by every session, used to notice when the database file changes."""
    return {}


@st.cache_data(show_spinner=False)
def filter_options():
    """Every value the filter dropdowns can offer, read from the database."""
    return {
        "condition": load_table("SELECT DISTINCT condition FROM subject ORDER BY 1")["condition"].tolist(),
        "treatment": load_table("SELECT DISTINCT treatment FROM subject ORDER BY 1")["treatment"].tolist(),
        "sample_type": load_table("SELECT DISTINCT sample_type FROM sample ORDER BY 1")["sample_type"].tolist(),
        "timepoint": [int(t) for t in load_table(
            "SELECT DISTINCT time_from_treatment_start AS t FROM sample "
            "WHERE time_from_treatment_start IS NOT NULL ORDER BY 1")["t"]],
        "population": db.populations(),
    }


def table_exists(name):
    """True if the pipeline has written this result table."""
    return not load_table("SELECT name FROM sqlite_master WHERE type='table' AND name=?", (name,)).empty


def tidy_results(df):
    """Make a Part 3 result table look the same whether it was computed live or read back from SQLite."""
    for col in NUMERIC_RESULT_COLUMNS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col])
    if "significant" in df.columns:
        df["significant"] = df["significant"].astype(bool)
    return df


def pick(options, preferred):
    """Index of the preferred option in a dropdown, or 0 if it is not there."""
    for i, o in enumerate(options):
        if str(o).lower() == str(preferred).lower():
            return i
    return 0


def download(df, label, name):
    """A download button for a table as CSV."""
    st.download_button(label, df.to_csv(index=False).encode(), file_name=name, mime="text/csv")


def highlight_significant(row):
    """Row style for the statistics tables: shaded when the result is significant."""
    return [HIGHLIGHT if row["significant"] else "" for _ in row]


def facet_titles(fig, names):
    """Replace Plotly's "population=b_cell" facet titles with the display names."""
    fig.for_each_annotation(lambda a: a.update(text=names.get(a.text.split("=")[-1], a.text.split("=")[-1])))


if not DB_PATH.exists():
    # Hosted deploys only run `streamlit run`, so build the database on first start.
    with st.spinner("No database yet, running the pipeline..."):
        for script in ("load_data.py", "run_analysis.py"):
            proc = subprocess.run([sys.executable, str(DB_PATH.parent / script)],
                                  capture_output=True, text=True, cwd=DB_PATH.parent)
            if proc.returncode != 0:
                st.error(f"`{script}` failed:\n\n```\n{proc.stderr[-2000:]}\n```")
                st.stop()

# The caches above are keyed on arguments only, so drop them when the pipeline rebuilds the database.
stamp = DB_PATH.stat().st_mtime_ns
if last_seen().get("db") != stamp:
    st.cache_data.clear()
    last_seen()["db"] = stamp

opts = filter_options()
if not opts["population"] or not opts["condition"]:
    st.warning("The database has no populations or no subjects, so there is nothing to show. "
               "Check the notes printed by `make pipeline`.")
    st.stop()
pop_names = dict(load_table("SELECT population, display_name FROM population").itertuples(index=False, name=None))

st.title("Loblaw Bio: miraclib and immune cell populations")
st.caption("Immune cell counts per sample, the responder comparison, and the baseline subset from the brief.")

kpi = load_table(
    "SELECT (SELECT COUNT(*) FROM sample) AS samples, (SELECT COUNT(*) FROM subject) AS subjects, "
    "(SELECT COUNT(*) FROM project) AS projects, (SELECT COUNT(*) FROM population) AS populations"
).iloc[0]
for col, label, val in zip(st.columns(4), ["Samples", "Subjects", "Projects", "Cell populations"],
                           [kpi.samples, kpi.subjects, kpi.projects, kpi.populations]):
    col.metric(label, f"{val:,}")

tab2, tab3, tab4, tab_schema = st.tabs([
    "Part 2: Frequency overview", "Part 3: Responders vs non-responders",
    "Part 4: Baseline subset", "Database schema",
])


# ---- Part 2
with tab2:
    st.subheader("Relative frequency of each population in each sample")
    st.write("Counts are summed per sample and each population is shown as a percentage of that total.")
    freq = frequency()

    if freq.empty:
        st.warning("No sample has a count for every population, so the frequency table is empty. "
                   "The notes printed by `make pipeline` say which samples were left out.")
    else:
        c1, c2 = st.columns([2, 1])
        pops = c1.multiselect("Populations", opts["population"], default=opts["population"], key="p2pops")
        search = c2.text_input("Filter by sample id", placeholder="e.g. sample00042")
        shown = freq[freq["population"].isin(pops)]
        if search:
            shown = shown[shown["sample"].str.contains(search.strip(), case=False, regex=False, na=False)]
        st.dataframe(shown, width="stretch", hide_index=True, height=360)
        st.caption(f"{len(shown):,} rows shown of {len(freq):,}.")
        download(freq, "Download full summary table (CSV)", "part2_frequency_summary.csv")

        left, right = st.columns(2)
        with left:
            st.markdown("**Distribution across all samples**")
            fig = px.box(freq[freq["population"].isin(pops)], x="population", y="percentage", points=False,
                         labels={"percentage": "% of total cells", "population": ""},
                         category_orders={"population": opts["population"]})
            fig.update_traces(marker_color=BLUE, line_width=1.5, fillcolor="rgba(42,120,214,0.15)")
            fig.update_layout(**LAYOUT, showlegend=False, height=380)
            st.plotly_chart(fig, width="stretch")
        with right:
            st.markdown("**One sample**")
            chosen = st.selectbox("Sample", freq["sample"].unique().tolist(), index=0)
            one = freq[freq["sample"] == chosen].sort_values("percentage")
            fig = px.bar(one, x="percentage", y="population", orientation="h", text="percentage",
                         labels={"percentage": "% of total cells", "population": ""})
            fig.update_traces(marker_color=BLUE, texttemplate="%{text:.1f}%", textposition="outside", width=0.55)
            fig.update_layout(**LAYOUT, height=380, xaxis_range=[0, one["percentage"].max() * 1.35])
            st.plotly_chart(fig, width="stretch")
            st.caption(f"Total cells in {chosen}: {int(one['total_count'].iloc[0]):,}")


# ---- Part 3
with tab3:
    st.subheader("Do responders and non-responders differ in immune cell composition?")
    with st.expander("Cohort filters", expanded=False):
        a, b, c = st.columns(3)
        cond = a.selectbox("Condition", opts["condition"], index=pick(opts["condition"], "melanoma"), key="p3c")
        trt = b.selectbox("Treatment", opts["treatment"], index=pick(opts["treatment"], "miraclib"), key="p3t")
        styp = c.selectbox("Sample type", opts["sample_type"], index=pick(opts["sample_type"], "PBMC"), key="p3s")
        st.caption("Defaults are the brief's cohort: melanoma, miraclib, PBMC, response recorded.")

    is_default = (cond.lower(), trt.lower(), styp.upper()) == ("melanoma", "miraclib", "PBMC")
    if is_default and table_exists("result_response_stats"):
        cohort = load_table("SELECT * FROM result_response_cohort")
        results = load_table("SELECT * FROM result_response_stats")
        traj = load_table("SELECT * FROM result_response_trajectory") if table_exists("result_response_trajectory") else pd.DataFrame()
        source = "Read from the tables `make pipeline` wrote."
    else:
        part3 = live_part3(cond, trt, styp)
        cohort, results, traj = part3["response_cohort"], part3["response_stats"], part3["response_trajectory"]
        source = "Computed live from the database for these filters."
    results, traj = tidy_results(results), tidy_results(traj)

    n_r = cohort.loc[cohort["group"] == "responder", "subject_key"].nunique() if not cohort.empty else 0
    n_n = cohort.loc[cohort["group"] == "non-responder", "subject_key"].nunique() if not cohort.empty else 0

    if cohort.empty:
        st.warning("No samples with a recorded response match these filters.")
    elif results.empty:
        st.warning(f"This cohort has {n_r} responder(s) and {n_n} non-responder(s), which is too few "
                   "to run the comparison. At least two samples per group are required.")
    else:
        k1, k2, k3 = st.columns(3)
        k1.metric("Responders (subjects)", f"{n_r:,}")
        k2.metric("Non-responders (subjects)", f"{n_n:,}")
        k3.metric("Samples compared", f"{cohort['sample'].nunique():,}")
        st.caption(source)

        st.markdown("#### Conclusion")
        text = "**All timepoints pooled (per sample):** " + stats.summarize(results, "sample")
        if not traj.empty:
            text += "\n\n**Change from baseline (per subject):** " + stats.summarize(traj, "change_from_baseline")
        st.info(text)

        st.markdown("#### Relative frequency by population")
        unit_choice = st.radio("Unit", ["Every sample", "Per-subject mean"], horizontal=True)
        plot_df = cohort
        if unit_choice == "Per-subject mean":
            plot_df = cohort.groupby(["subject_key", "population", "group"], as_index=False)["percentage"].mean()
        fig = px.box(plot_df, x="group", y="percentage", color="group", facet_col="population",
                     facet_col_wrap=5, points="outliers", color_discrete_map=COLORS,
                     category_orders={"group": ["responder", "non-responder"], "population": opts["population"]},
                     labels={"percentage": "% of total cells", "group": ""})
        fig.update_traces(line_width=1.5, marker=dict(size=4, opacity=0.5), width=0.5)
        fig.update_yaxes(matches=None, showticklabels=True)
        fig.update_xaxes(showticklabels=False, title="")
        facet_titles(fig, pop_names)
        fig.update_layout(**LAYOUT, height=440, legend_title_text="")
        st.plotly_chart(fig, width="stretch")

        st.markdown("#### Statistical tests")
        st.write("Two-sided Mann-Whitney U per population, Benjamini-Hochberg adjusted. "
                 "Positive Cliff's delta means responders tend to be higher, and the direction column "
                 "follows its sign. Welch's t is shown for comparison.")
        unit_key = "sample" if unit_choice == "Every sample" else "subject"
        cols = ["population", "n_responders", "n_non_responders", "median_responders", "median_non_responders",
                "median_difference", "cliffs_delta", "p_value", "p_adjusted_bh", "welch_p_value",
                "significant", "direction"]
        tbl = results[results["unit"] == unit_key][cols].reset_index(drop=True)
        st.dataframe(
            tbl.style.format(STATS_FORMAT, na_rep="").apply(highlight_significant, axis=1),
            width="stretch", hide_index=True,
        )
        download(results, "Download all test results (CSV)", "part3_response_stats.csv")

        if cohort[TIME].nunique() > 1:
            st.markdown("#### Does the difference emerge over treatment?")
            st.write("Median relative frequency per group at each timepoint, and a test on each "
                     "subject's change from time 0 to that subject's last sample.")
            med = (cohort.dropna(subset=[TIME])
                   .groupby(["population", "group", TIME], as_index=False)["percentage"]
                   .median())
            fig = px.line(med, x=TIME, y="percentage", color="group", facet_col="population",
                          facet_col_wrap=5, markers=True, color_discrete_map=COLORS,
                          category_orders={"group": ["responder", "non-responder"], "population": opts["population"]},
                          labels={"percentage": "median % of total cells",
                                  TIME: "days from treatment start", "group": ""})
            fig.update_traces(line_width=2, marker_size=8)
            fig.update_yaxes(matches=None, showticklabels=True)
            fig.update_xaxes(tickvals=sorted(med[TIME].unique()))
            facet_titles(fig, pop_names)
            fig.update_layout(**LAYOUT, height=400, legend_title_text="")
            fig.update_layout(legend=dict(y=-0.3))
            st.plotly_chart(fig, width="stretch")

            if traj.empty:
                reason = stats.baseline_skip_reason(cohort) or stats.TOO_FEW_SUBJECTS
                st.caption(f"The change-from-baseline test was not run: {reason}.")
            else:
                tcols = ["population", "n_responders", "n_non_responders", "median_responders",
                         "median_non_responders", "cliffs_delta", "p_value", "p_adjusted_bh", "significant", "direction"]
                st.markdown(f"**Change from time 0 to each subject's last sample (up to day "
                            f"{int(traj['final_time'].iloc[0])})** (percentage-point changes)")
                st.dataframe(
                    traj[tcols].rename(columns={"median_responders": "median_change_responders",
                                                "median_non_responders": "median_change_non_responders"})
                    .style.format(TRAJECTORY_FORMAT, na_rep="").apply(highlight_significant, axis=1),
                    width="stretch", hide_index=True,
                )
                download(traj, "Download trajectory results (CSV)", "part3_response_trajectory.csv")


# ---- Part 4
with tab4:
    st.subheader("Baseline samples from treated patients")
    a, b, c, d = st.columns(4)
    cond4 = a.selectbox("Condition", opts["condition"], index=pick(opts["condition"], "melanoma"), key="p4c")
    trt4 = b.selectbox("Treatment", opts["treatment"], index=pick(opts["treatment"], "miraclib"), key="p4t")
    styp4 = c.selectbox("Sample type", opts["sample_type"], index=pick(opts["sample_type"], "PBMC"), key="p4s")
    t4 = d.selectbox("Time from treatment start", opts["timepoint"], index=pick(opts["timepoint"], 0), key="p4time")

    is_default4 = (cond4.lower(), trt4.lower(), styp4.upper(), t4) == ("melanoma", "miraclib", "PBMC", 0)
    if t4 is None:
        baseline = pd.DataFrame()  # The selectbox had nothing to offer: no sample has a timepoint.
    elif is_default4 and table_exists("result_baseline_samples"):
        baseline = load_table("SELECT * FROM result_baseline_samples")
        st.caption("Read from the tables `make pipeline` wrote.")
    else:
        baseline = live_part4(cond4, trt4, styp4, t4)
        st.caption("Computed live from the database for these filters.")

    if baseline.empty:
        st.warning("No samples match these filters." if t4 is not None else
                   "No sample in this database has a timepoint, so there is no baseline subset to show.")
    else:
        m1, m2 = st.columns(2)
        m1.metric("Matching samples", f"{len(baseline):,}")
        m2.metric("Distinct subjects", f"{baseline['subject_key'].nunique():,}")

        def small_bar(df, x, y, title):
            if df.empty:
                st.markdown(f"**{title}**")
                st.caption("No data for this selection.")
                return
            fig = px.bar(df, x=x, y=y, text=y, labels={x: "", y: y})
            fig.update_traces(marker_color=BLUE, textposition="outside", width=0.5)
            fig.update_layout(**LAYOUT, title=title, height=300, yaxis_range=[0, df[y].max() * 1.2])
            st.plotly_chart(fig, width="stretch")
            st.dataframe(df, hide_index=True, width="stretch")

        g1, g2, g3 = st.columns(3)
        with g1:
            small_bar(subsets.samples_per_project(baseline), "project", "samples", "Samples per project")
        with g2:
            small_bar(subsets.subjects_by_response(baseline), "response", "subjects", "Subjects by response")
        with g3:
            small_bar(subsets.subjects_by_sex(baseline), "sex", "subjects", "Subjects by sex")

        st.markdown("#### Sample list")
        st.dataframe(baseline, width="stretch", hide_index=True, height=320)
        download(baseline, "Download sample list (CSV)", "part4_baseline_samples.csv")


# ---- Schema
with tab_schema:
    st.subheader("How the data is modeled")
    st.write("One wide CSV row per sample becomes project, subject and sample tables plus a long cell_count "
             "table. The Part 2 summary is a view over the counts.")
    st.code((DB_PATH.parent / "schema.sql").read_text(), language="sql")
    tables = load_table("SELECT name, type FROM sqlite_master WHERE type IN ('table','view') "
                        "AND name NOT LIKE 'sqlite_%' ORDER BY type, name")
    tables["rows"] = [int(load_table(f'SELECT COUNT(*) AS n FROM "{name}"')["n"].iloc[0]) for name in tables["name"]]
    st.dataframe(tables, hide_index=True, width="stretch")
