#!/usr/bin/env python3
"""Load cell-count.csv into cell_counts.db."""
import csv
import os
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "cell-count.csv"
SCHEMA_PATH = ROOT / "schema.sql"
DB_PATH = ROOT / "cell_counts.db"

POPULATIONS = ("b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte")
SUBJECT_COLUMNS = ("project", "condition", "age", "sex", "treatment", "response")
SAMPLE_COLUMNS = ("subject", "sample", "sample_type", "time_from_treatment_start")


def read_csv(path):
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        missing = set(SUBJECT_COLUMNS + SAMPLE_COLUMNS + POPULATIONS) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path.name} is missing columns: {sorted(missing)}")
        return list(reader)


def split_rows(rows):
    subjects = {}
    samples = []
    counts = []
    for line_no, row in enumerate(rows, start=2):
        subject = (
            row["subject"],
            row["project"],
            row["condition"],
            int(row["age"]),
            row["sex"],
            row["treatment"],
            row["response"] or None,
        )
        existing = subjects.setdefault(row["subject"], subject)
        if existing != subject:
            raise ValueError(
                f"line {line_no}: subject {row['subject']} has attributes {subject[1:]} "
                f"that conflict with an earlier row {existing[1:]}"
            )

        samples.append(
            (row["sample"], row["subject"], row["sample_type"], int(row["time_from_treatment_start"]))
        )
        for population_id, name in enumerate(POPULATIONS, start=1):
            counts.append((row["sample"], population_id, int(row[name])))

    projects = sorted({(s[1],) for s in subjects.values()})
    return projects, list(subjects.values()), samples, counts


def build_database(db_path, projects, subjects, samples, counts):
    conn = sqlite3.connect(db_path)
    try:
        conn.execute("PRAGMA foreign_keys = ON")
        conn.executescript(SCHEMA_PATH.read_text())
        with conn:
            conn.executemany("INSERT INTO projects VALUES (?)", projects)
            conn.executemany(
                "INSERT INTO populations (population_id, name) VALUES (?, ?)",
                list(enumerate(POPULATIONS, start=1)),
            )
            conn.executemany("INSERT INTO subjects VALUES (?, ?, ?, ?, ?, ?, ?)", subjects)
            conn.executemany("INSERT INTO samples VALUES (?, ?, ?, ?)", samples)
            conn.executemany("INSERT INTO cell_counts VALUES (?, ?, ?)", counts)
        return {
            table: conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
            for table in ("projects", "subjects", "samples", "populations", "cell_counts")
        }
    finally:
        conn.close()


def main():
    rows = read_csv(CSV_PATH)
    projects, subjects, samples, counts = split_rows(rows)

    # write db via temp file
    tmp_path = DB_PATH.with_suffix(".db.tmp")
    tmp_path.unlink(missing_ok=True)
    try:
        table_sizes = build_database(tmp_path, projects, subjects, samples, counts)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise
    os.replace(tmp_path, DB_PATH)

    print(f"Loaded {len(rows)} rows from {CSV_PATH.name} into {DB_PATH.name}")
    for table, n in table_sizes.items():
        print(f"  {table:<12} {n:>7}")


if __name__ == "__main__":
    main()
