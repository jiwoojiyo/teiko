from analysis.db import query

BASELINE_FILTERS = {
    "condition": "melanoma",
    "treatment": "miraclib",
    "sample_type": "PBMC",
    "time_from_treatment_start": 0,
}

FILTERABLE_COLUMNS = {
    "project_id": "su.project_id",
    "condition": "su.condition",
    "sex": "su.sex",
    "treatment": "su.treatment",
    "response": "su.response",
    "sample_type": "sa.sample_type",
    "time_from_treatment_start": "sa.time_from_treatment_start",
}

SUBSET_CTE = """
WITH subset AS (
    SELECT sa.sample_id, sa.subject_id, su.project_id, su.condition, su.age, su.sex,
           su.treatment, su.response, sa.sample_type, sa.time_from_treatment_start
    FROM samples AS sa
    JOIN subjects AS su ON su.subject_id = sa.subject_id
    WHERE {where}
)
"""


def _where(filters):
    clauses, params = [], []
    for column, value in filters.items():
        if column not in FILTERABLE_COLUMNS:
            raise ValueError(f"cannot filter on {column!r}")
        if value is None:
            continue
        sql_column = FILTERABLE_COLUMNS[column]
        if isinstance(value, (list, tuple, set)):
            values = list(value)
            if not values:
                clauses.append("0")
                continue
            clauses.append(f"{sql_column} IN ({', '.join('?' * len(values))})")
            params.extend(values)
        else:
            clauses.append(f"{sql_column} = ?")
            params.append(value)
    return " AND ".join(clauses) or "1", params


def _subset_query(conn, select_sql, filters):
    where, params = _where(filters)
    return query(conn, SUBSET_CTE.format(where=where) + select_sql, params)


def subset_samples(conn, filters=BASELINE_FILTERS):
    return _subset_query(conn, "SELECT * FROM subset ORDER BY sample_id", filters)


def samples_per_project(conn, filters=BASELINE_FILTERS):
    return _subset_query(
        conn,
        """
        SELECT project_id AS project, COUNT(*) AS n_samples, COUNT(DISTINCT subject_id) AS n_subjects
        FROM subset GROUP BY project_id ORDER BY project_id
        """,
        filters,
    )


def subjects_by_response(conn, filters=BASELINE_FILTERS):
    return _subset_query(
        conn,
        """
        SELECT CASE response WHEN 'yes' THEN 'responder' WHEN 'no' THEN 'non-responder'
                             ELSE 'not applicable' END AS response,
               COUNT(DISTINCT subject_id) AS n_subjects, COUNT(*) AS n_samples
        FROM subset GROUP BY subset.response ORDER BY subset.response DESC
        """,
        filters,
    )


def subjects_by_sex(conn, filters=BASELINE_FILTERS):
    return _subset_query(
        conn,
        """
        SELECT CASE sex WHEN 'M' THEN 'male' WHEN 'F' THEN 'female' END AS sex,
               COUNT(DISTINCT subject_id) AS n_subjects, COUNT(*) AS n_samples
        FROM subset GROUP BY subset.sex ORDER BY subset.sex DESC
        """,
        filters,
    )


AVG_B_CELL_SQL = """
SELECT AVG(cc.count) AS mean_b_cell, COUNT(*) AS n_samples
FROM cell_counts AS cc
JOIN populations AS p ON p.population_id = cc.population_id
JOIN samples AS sa ON sa.sample_id = cc.sample_id
JOIN subjects AS su ON su.subject_id = sa.subject_id
WHERE p.name = 'b_cell'
  AND su.condition = 'melanoma'
  AND su.sex = 'M'
  AND su.response = 'yes'
  AND sa.time_from_treatment_start = 0
"""


def mean_b_cells_melanoma_male_responders_baseline(conn):
    row = query(conn, AVG_B_CELL_SQL).iloc[0]
    return float(row["mean_b_cell"]), int(row["n_samples"])
