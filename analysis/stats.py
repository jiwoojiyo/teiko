import warnings

import pandas as pd
import plotly.graph_objects as go
import statsmodels.formula.api as smf
from plotly.subplots import make_subplots
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from statsmodels.tools.sm_exceptions import ConvergenceWarning

from analysis.db import query

ALPHA = 0.05
RESPONSE_LABELS = {"yes": "Responder", "no": "Non-responder"}
RESPONSE_COLORS = {"yes": "#1b9e77", "no": "#d95f02"}

COHORT_SQL = """
SELECT sample_id, subject_id, project_id, time_from_treatment_start, response,
       population_id, population, percentage
FROM sample_population_frequencies
WHERE condition = ? AND treatment = ? AND sample_type = ?
  AND response IN ('yes', 'no')
ORDER BY population_id, sample_id
"""


def cohort_frequencies(conn, condition="melanoma", treatment="miraclib", sample_type="PBMC"):
    return query(conn, COHORT_SQL, (condition, treatment, sample_type))


def _bh_adjust(results):
    results["q_value"] = multipletests(results["p_value"], method="fdr_bh")[1]
    results["significant"] = results["q_value"] < ALPHA
    return results


def _mann_whitney(population, group):
    responders = group.loc[group["response"] == "yes", "percentage"]
    non_responders = group.loc[group["response"] == "no", "percentage"]
    u, p = mannwhitneyu(responders, non_responders, alternative="two-sided")
    n1, n2 = len(responders), len(non_responders)
    return {
        "population": population,
        "n_responders": n1,
        "n_non_responders": n2,
        "median_responders": responders.median(),
        "median_non_responders": non_responders.median(),
        "median_difference": responders.median() - non_responders.median(),
        # positive when responders higher
        "rank_biserial": 2 * u / (n1 * n2) - 1,
        "u_statistic": u,
        "p_value": p,
    }


def compare_responders(df):
    rows = [_mann_whitney(pop, group) for pop, group in df.groupby("population", sort=False)]
    return _bh_adjust(pd.DataFrame(rows))


def compare_by_timepoint(df):
    rows = [
        {"time_from_treatment_start": time, **_mann_whitney(pop, group)}
        for (time, pop), group in df.groupby(["time_from_treatment_start", "population"], sort=False)
    ]
    return _bh_adjust(pd.DataFrame(rows).sort_values("time_from_treatment_start", kind="stable"))


def mixed_model(df):
    rows = []
    for population, group in df.groupby("population", sort=False):
        data = group.assign(responder=(group["response"] == "yes").astype(int))
        model = smf.mixedlm(
            "percentage ~ responder + C(time_from_treatment_start)", data, groups=data["subject_id"]
        )
        # singular covariance warnings are ok
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            warnings.filterwarnings("ignore", message="Random effects covariance is singular")
            fit = model.fit(reml=True)
        ci_low, ci_high = fit.conf_int().loc["responder"]
        rows.append({
            "population": population,
            "n_samples": len(data),
            "n_subjects": data["subject_id"].nunique(),
            "responder_effect": fit.params["responder"],
            "ci_low": ci_low,
            "ci_high": ci_high,
            "p_value": fit.pvalues["responder"],
            "subject_variance": float(fit.cov_re.iloc[0, 0]),
            "converged": fit.converged,
        })
    return _bh_adjust(pd.DataFrame(rows))


def responder_boxplot(df, stats=None, title="Melanoma, miraclib, PBMC: responders vs non-responders"):
    populations = list(df["population"].unique())
    titles = populations
    if stats is not None:
        by_pop = stats.set_index("population")
        titles = [
            f"{pop}<br><sup>p = {by_pop.at[pop, 'p_value']:.2g}, q = {by_pop.at[pop, 'q_value']:.2g}</sup>"
            for pop in populations
        ]

    fig = make_subplots(rows=1, cols=len(populations), subplot_titles=titles, horizontal_spacing=0.04)
    for col, pop in enumerate(populations, start=1):
        subset = df[df["population"] == pop]
        for code, label in RESPONSE_LABELS.items():
            fig.add_trace(
                go.Box(
                    y=subset.loc[subset["response"] == code, "percentage"],
                    name=label,
                    legendgroup=code,
                    showlegend=col == 1,
                    marker=dict(color=RESPONSE_COLORS[code], size=3, opacity=0.35),
                    boxpoints="all",
                    jitter=0.4,
                    pointpos=0,
                ),
                row=1,
                col=col,
            )
    fig.update_xaxes(showticklabels=False)
    fig.update_yaxes(title_text="Relative frequency (%)", row=1, col=1)
    fig.update_layout(
        title=title,
        height=520,
        legend=dict(orientation="h", y=-0.08),
        margin=dict(t=110),
        template="plotly_white",
    )
    return fig
