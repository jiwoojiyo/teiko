#!/usr/bin/env python3
"""Run analyses and write outputs/."""
from contextlib import closing

from analysis import stats, subsets
from analysis.db import ROOT, connect
from analysis.overview import frequency_summary

OUTPUT_DIR = ROOT / "outputs"


def write_csv(df, name):
    path = OUTPUT_DIR / name
    df.to_csv(path, index=False)
    return path.relative_to(ROOT)


def run_part2(conn):
    frequencies = frequency_summary(conn)
    path = write_csv(frequencies, "frequency_summary.csv")
    print(f"Part 2: {len(frequencies)} rows ({frequencies['sample'].nunique()} samples) -> {path}")


def run_part3(conn):
    cohort = stats.cohort_frequencies(conn)
    comparison = stats.compare_responders(cohort)
    write_csv(comparison, "responder_stats.csv")
    mixed = stats.mixed_model(cohort)
    write_csv(mixed, "responder_stats_mixed_model.csv")
    by_time = stats.compare_by_timepoint(cohort)
    write_csv(by_time, "responder_stats_by_timepoint.csv")

    fig = stats.responder_boxplot(cohort, comparison)
    fig.write_html(OUTPUT_DIR / "responder_boxplot.html", include_plotlyjs="cdn")
    try:
        fig.write_image(OUTPUT_DIR / "responder_boxplot.png", width=1400, height=560)
    except Exception as exc:  # png export is optional
        print(f"  warning: PNG export failed ({exc}); HTML boxplot still written")

    n_samples = cohort["sample_id"].nunique()
    n_subjects = cohort["subject_id"].nunique()
    print(f"Part 3: melanoma/miraclib/PBMC, {n_samples} samples from {n_subjects} subjects")
    print("  Mann-Whitney (q = BH, 5 populations):")
    for row in comparison.itertuples():
        print(
            f"    {row.population:<11} median diff {row.median_difference:+.2f} pp  "
            f"p={row.p_value:.4f}  q={row.q_value:.4f}{'  *' if row.significant else ''}"
        )
    print("  Mixed model (random subject, day fixed):")
    for row in mixed.itertuples():
        print(
            f"    {row.population:<11} effect {row.responder_effect:+.2f} pp "
            f"[{row.ci_low:+.2f}, {row.ci_high:+.2f}]  p={row.p_value:.4f}  q={row.q_value:.4f}"
            f"{'  *' if row.significant else ''}"
        )
    return comparison, mixed, by_time


def run_part4(conn):
    baseline = subsets.subset_samples(conn)
    write_csv(baseline, "baseline_samples.csv")
    per_project = subsets.samples_per_project(conn)
    write_csv(per_project, "baseline_samples_per_project.csv")
    by_response = subsets.subjects_by_response(conn)
    write_csv(by_response, "baseline_subjects_by_response.csv")
    by_sex = subsets.subjects_by_sex(conn)
    write_csv(by_sex, "baseline_subjects_by_sex.csv")
    mean_b, n_b = subsets.mean_b_cells_melanoma_male_responders_baseline(conn)

    print(f"Part 4: melanoma/miraclib/PBMC baseline samples: {len(baseline)}")
    print("  samples per project: " + ", ".join(f"{r.project}={r.n_samples}" for r in per_project.itertuples()))
    print("  subjects by response: " + ", ".join(f"{r.response}={r.n_subjects}" for r in by_response.itertuples()))
    print("  subjects by sex: " + ", ".join(f"{r.sex}={r.n_subjects}" for r in by_sex.itertuples()))
    print(f"  mean B cells, melanoma male responders at time 0 (n={n_b} samples): {mean_b:.2f}")
    return baseline, per_project, by_response, by_sex, mean_b, n_b


def write_summary(part3, part4):
    comparison, mixed, by_time = part3
    baseline, per_project, by_response, by_sex, mean_b, n_b = part4

    def table(df, columns):
        header = "| " + " | ".join(columns) + " |\n|" + "---|" * len(columns) + "\n"
        return header + "".join("| " + " | ".join(str(v) for v in row) + " |\n" for row in df[columns].values)

    fmt = comparison.assign(
        median_difference=comparison["median_difference"].round(3),
        rank_biserial=comparison["rank_biserial"].round(3),
        p_value=comparison["p_value"].map("{:.4f}".format),
        q_value=comparison["q_value"].map("{:.4f}".format),
    )
    fmt_mixed = mixed.assign(
        responder_effect=mixed["responder_effect"].round(3),
        ci=[f"[{lo:.3f}, {hi:.3f}]" for lo, hi in zip(mixed["ci_low"], mixed["ci_high"])],
        p_value=mixed["p_value"].map("{:.4f}".format),
        q_value=mixed["q_value"].map("{:.4f}".format),
    )
    sig_baseline = by_time[(by_time["time_from_treatment_start"] == 0) & by_time["significant"]]

    text = (
        "# Run summary\n\n"
        "## Responders vs non-responders (melanoma, miraclib, PBMC)\n\n"
        + table(fmt, ["population", "median_difference", "p_value", "q_value", "significant"])
        + "\nMixed model (responder + day, random subject):\n\n"
        + table(fmt_mixed, ["population", "responder_effect", "ci", "p_value", "q_value", "significant"])
        + f"\nSignificant at day 0 only: {', '.join(sig_baseline['population']) or 'none'}\n\n"
        "## Baseline subset (same filters, day 0)\n\n"
        f"Samples: {len(baseline)}\n\n"
        + table(per_project, ["project", "n_samples"])
        + "\n" + table(by_response, ["response", "n_subjects"])
        + "\n" + table(by_sex, ["sex", "n_subjects"])
        + f"\nMean B cells — melanoma, male, responder, day 0 (all types/treatments, n={n_b}): {mean_b:.2f}\n"
    )
    (OUTPUT_DIR / "summary.md").write_text(text)


def main():
    OUTPUT_DIR.mkdir(exist_ok=True)
    with closing(connect()) as conn:
        run_part2(conn)
        part3 = run_part3(conn)
        part4 = run_part4(conn)
    write_summary(part3, part4)
    print(f"All outputs written to {OUTPUT_DIR.relative_to(ROOT)}/")


if __name__ == "__main__":
    main()
