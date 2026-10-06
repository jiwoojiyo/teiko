from analysis.db import query

FREQUENCY_SQL = """
SELECT sample_id AS sample, total_count, population, count, percentage
FROM sample_population_frequencies
ORDER BY sample_id, population_id
"""


def frequency_summary(conn):
    return query(conn, FREQUENCY_SQL)
