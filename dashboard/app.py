import sys
from contextlib import closing
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import load_data
import plotly.express as px
import streamlit as st
from analysis import stats, subsets
from analysis.db import DB_PATH, connect, query

st.set_page_config(page_title="Miraclib immune profiling", layout="wide")

if not DB_PATH.exists():
    # build db from bundled csv
    load_data.main()

POPULATION_ORDER = list(load_data.POPULATIONS)
ANY = "Any"


def run(fn, *args, **kwargs):
    with closing(connect()) as conn:
        return fn(conn, *args, **kwargs)


@st.cache_data
def frequencies_with_metadata():
    return run(
        query,
        """
        SELECT sample_id AS sample, total_count, population, count, percentage,
               project_id AS project, subject_id AS subject, condition, sex, treatment,
               COALESCE(response, 'n/a') AS response, sample_type, time_from_treatment_start
        FROM sample_population_frequencies
        ORDER BY sample_id, population_id
        """,
    )


@st.cache_data
def responder_cohort():
    return run(stats.cohort_frequencies)


@st.cache_data
def responder_tests(time_point):
    cohort = responder_cohort()
    if time_point is not None:
        cohort = cohort[cohort["time_from_treatment_start"] == time_point]
    return cohort, stats.compare_responders(cohort)


@st.cache_data
def repeat_measure_tables():
    cohort = responder_cohort()
    return stats.mixed_model(cohort), stats.compare_by_timepoint(cohort)


@st.cache_data
def subset_tables(filters):
    filters = dict(filters)
    return (
        run(subsets.subset_samples, filters),
        run(subsets.samples_per_project, filters),
        run(subsets.subjects_by_response, filters),
        run(subsets.subjects_by_sex, filters),
    )


@st.cache_data
def mean_b_cells():
    return run(subsets.mean_b_cells_melanoma_male_responders_baseline)


def multiselect_filter(df, column, label, container):
    options = sorted(df[column].unique().tolist())
    chosen = container.multiselect(label, options, default=options)
    return df[df[column].isin(chosen)]


def stats_table(df, extra_columns=()):
    columns = [*extra_columns, "population", "n_responders", "n_non_responders", "median_responders",
               "median_non_responders", "median_difference", "rank_biserial", "p_value", "q_value",
               "significant"]
    st.dataframe(
        df[columns],
        hide_index=True,
        use_container_width=True,
        column_config={
            "median_responders": st.column_config.NumberColumn("median % (resp.)", format="%.2f"),
            "median_non_responders": st.column_config.NumberColumn("median % (non-resp.)", format="%.2f"),
            "median_difference": st.column_config.NumberColumn("median diff (pp)", format="%+.2f"),
            "rank_biserial": st.column_config.NumberColumn("rank-biserial r", format="%+.3f"),
            "p_value": st.column_config.NumberColumn("p", format="%.4f"),
            "q_value": st.column_config.NumberColumn("q (BH)", format="%.4f"),
        },
    )


st.title("Miraclib trial: immune cell populations")
st.caption(f"Database: `{DB_PATH.name}`")

tab_overview, tab_response, tab_subsets = st.tabs(
    ["Cell frequencies", "Responders vs non-responders", "Baseline subset"]
)

with tab_overview:
    st.subheader("Relative frequency by sample")
    data = frequencies_with_metadata()

    cols = st.columns(5)
    filtered = multiselect_filter(data, "project", "Project", cols[0])
    filtered = multiselect_filter(filtered, "condition", "Condition", cols[1])
    filtered = multiselect_filter(filtered, "treatment", "Treatment", cols[2])
    filtered = multiselect_filter(filtered, "sample_type", "Sample type", cols[3])
    filtered = multiselect_filter(filtered, "time_from_treatment_start", "Day", cols[4])
    sample_search = st.text_input("Sample id contains", "")
    if sample_search:
        filtered = filtered[filtered["sample"].str.contains(sample_search.strip(), case=False, regex=False)]

    m1, m2, m3 = st.columns(3)
    m1.metric("Samples", f"{filtered['sample'].nunique():,}")
    m2.metric("Subjects", f"{filtered['subject'].nunique():,}")
    m3.metric("Rows", f"{len(filtered):,}")

    show_metadata = st.toggle("Show metadata columns", value=False)
    table_columns = ["sample", "total_count", "population", "count", "percentage"]
    if show_metadata:
        table_columns += ["project", "subject", "condition", "sex", "treatment", "response",
                          "sample_type", "time_from_treatment_start"]
    st.dataframe(
        filtered[table_columns],
        hide_index=True,
        use_container_width=True,
        height=380,
        column_config={"percentage": st.column_config.NumberColumn("percentage", format="%.2f")},
    )
    st.download_button(
        "Download CSV",
        filtered[table_columns].to_csv(index=False),
        file_name="frequency_summary_filtered.csv",
        mime="text/csv",
    )

    st.markdown("#### Mean composition")
    group_by = st.selectbox(
        "Group by",
        ["condition", "treatment", "time_from_treatment_start", "sample_type", "response", "project", "sex"],
    )
    if filtered.empty:
        st.info("No rows match filters.")
    else:
        composition = (
            filtered.groupby([group_by, "population"], as_index=False)["percentage"].mean()
            .assign(**{group_by: lambda d: d[group_by].astype(str)})
        )
        st.plotly_chart(
            px.bar(
                composition, x=group_by, y="percentage", color="population",
                category_orders={"population": POPULATION_ORDER},
                labels={"percentage": "Mean relative frequency (%)"},
                template="plotly_white",
            ),
            use_container_width=True,
        )

with tab_response:
    st.subheader("Melanoma, miraclib, PBMC")
    time_label = st.radio(
        "Timepoints", ["All (days 0, 7, 14)", "Day 0", "Day 7", "Day 14"], horizontal=True
    )
    time_point = None if time_label.startswith("All") else int(time_label.split()[-1])
    cohort, comparison = responder_tests(time_point)

    st.caption(
        f"{cohort['sample_id'].nunique()} samples, {cohort['subject_id'].nunique()} subjects. "
        "Two-sided Mann-Whitney; q adjusts for five populations."
    )
    st.plotly_chart(
        stats.responder_boxplot(cohort, comparison, title=f"Responders vs non-responders ({time_label})"),
        use_container_width=True,
    )

    significant = comparison[comparison["significant"]]
    if not significant.empty:
        lines = [
            f"{r.population}: median diff {r.median_difference:+.2f} pp, q = {r.q_value:.3g}"
            for r in significant.itertuples()
        ]
        st.success("q < 0.05: " + "; ".join(lines))
    else:
        st.warning("Nothing significant at q < 0.05 after correction.")
    nominal = comparison[(comparison["p_value"] < stats.ALPHA) & ~comparison["significant"]]
    if not nominal.empty:
        lines = [
            f"{r.population}: p = {r.p_value:.3g}, q = {r.q_value:.3g}"
            for r in nominal.itertuples()
        ]
        st.info("p < 0.05 but q ≥ 0.05: " + "; ".join(lines))

    mixed, by_time = repeat_measure_tables()
    if time_point is None:
        mixed_sig = mixed[mixed["significant"]]
        if not mixed_sig.empty:
            lines = [
                f"{r.population}: {r.responder_effect:+.2f} pp (q = {r.q_value:.3g})"
                for r in mixed_sig.itertuples()
            ]
            st.info("Mixed model (all days): " + "; ".join(lines))

    stats_table(comparison)

    with st.expander("Mixed model and by-day tests"):
        st.markdown(
            "Each subject has samples at days 0, 7, and 14. "
            "Mixed model includes day; by-day table adjusts q over 15 tests."
        )
        st.dataframe(
            mixed, hide_index=True, use_container_width=True,
            column_config={
                "responder_effect": st.column_config.NumberColumn("effect (pp)", format="%+.3f"),
                "ci_low": st.column_config.NumberColumn("95% CI low", format="%+.3f"),
                "ci_high": st.column_config.NumberColumn("95% CI high", format="%+.3f"),
                "p_value": st.column_config.NumberColumn("p", format="%.4f"),
                "q_value": st.column_config.NumberColumn("q (BH)", format="%.4f"),
                "subject_variance": st.column_config.NumberColumn("subject variance", format="%.2e"),
            },
        )
        stats_table(by_time, extra_columns=("time_from_treatment_start",))

with tab_subsets:
    st.subheader("Melanoma, miraclib, PBMC at day 0")
    filters = tuple(subsets.BASELINE_FILTERS.items())
    baseline, per_project, by_response, by_sex = subset_tables(filters)

    st.metric("Samples in subset", len(baseline))
    c1, c2, c3 = st.columns(3)
    with c1:
        st.markdown("**By project**")
        st.dataframe(per_project[["project", "n_samples"]], hide_index=True, use_container_width=True)
    with c2:
        st.markdown("**By response (subjects)**")
        st.dataframe(by_response[["response", "n_subjects"]], hide_index=True, use_container_width=True)
    with c3:
        st.markdown("**By sex (subjects)**")
        st.dataframe(by_sex[["sex", "n_subjects"]], hide_index=True, use_container_width=True)

    mean_b, n_b = mean_b_cells()
    st.metric(
        "Mean B cells: melanoma, male, responder, day 0 (all treatments/types)",
        f"{mean_b:.2f}",
        help=f"n = {n_b} samples",
    )

    with st.expander("Sample list"):
        st.dataframe(baseline, hide_index=True, use_container_width=True)

    st.markdown("#### Other filters")
    meta = frequencies_with_metadata().drop_duplicates("sample")
    e = st.columns(4)
    choice = {
        "condition": e[0].selectbox("Condition", [ANY, *sorted(meta["condition"].unique())], index=0),
        "treatment": e[1].selectbox("Treatment", [ANY, *sorted(meta["treatment"].unique())], index=0),
        "sample_type": e[2].selectbox("Sample type", [ANY, *sorted(meta["sample_type"].unique())], index=0),
        "time_from_treatment_start": e[3].selectbox(
            "Day", [ANY, *sorted(meta["time_from_treatment_start"].unique().tolist())], index=0
        ),
    }
    custom = tuple((k, None if v == ANY else v) for k, v in choice.items())
    samples_c, project_c, response_c, sex_c = subset_tables(custom)
    st.metric("Matching samples", len(samples_c))
    x1, x2, x3 = st.columns(3)
    x1.dataframe(project_c, hide_index=True, use_container_width=True)
    x2.dataframe(response_c, hide_index=True, use_container_width=True)
    x3.dataframe(sex_c, hide_index=True, use_container_width=True)
