-- tables for project through counts

CREATE TABLE projects (
    project_id TEXT PRIMARY KEY
);

-- one row per patient
CREATE TABLE subjects (
    subject_id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL REFERENCES projects (project_id),
    condition  TEXT NOT NULL,
    age        INTEGER CHECK (age >= 0),
    sex        TEXT NOT NULL CHECK (sex IN ('M', 'F')),
    treatment  TEXT NOT NULL,
    -- empty response for controls
    response   TEXT CHECK (response IN ('yes', 'no'))
);

CREATE TABLE samples (
    sample_id                 TEXT PRIMARY KEY,
    subject_id                TEXT NOT NULL REFERENCES subjects (subject_id),
    sample_type               TEXT NOT NULL,
    time_from_treatment_start INTEGER NOT NULL CHECK (time_from_treatment_start >= 0)
);

CREATE TABLE populations (
    population_id INTEGER PRIMARY KEY,
    name          TEXT NOT NULL UNIQUE
);

CREATE TABLE cell_counts (
    sample_id     TEXT NOT NULL REFERENCES samples (sample_id),
    population_id INTEGER NOT NULL REFERENCES populations (population_id),
    count         INTEGER NOT NULL CHECK (count >= 0),
    PRIMARY KEY (sample_id, population_id)
) WITHOUT ROWID;

CREATE INDEX idx_subjects_project             ON subjects (project_id);
CREATE INDEX idx_subjects_condition_treatment ON subjects (condition, treatment);
CREATE INDEX idx_samples_subject              ON samples (subject_id);
CREATE INDEX idx_samples_type_time            ON samples (sample_type, time_from_treatment_start);
CREATE INDEX idx_cell_counts_population       ON cell_counts (population_id);

-- counts with pct per sample
CREATE VIEW sample_population_frequencies AS
SELECT
    sa.sample_id,
    sa.subject_id,
    su.project_id,
    su.condition,
    su.age,
    su.sex,
    su.treatment,
    su.response,
    sa.sample_type,
    sa.time_from_treatment_start,
    p.population_id,
    p.name AS population,
    cc.count,
    SUM(cc.count) OVER (PARTITION BY cc.sample_id) AS total_count,
    100.0 * cc.count / SUM(cc.count) OVER (PARTITION BY cc.sample_id) AS percentage
FROM cell_counts AS cc
JOIN samples AS sa ON sa.sample_id = cc.sample_id
JOIN subjects AS su ON su.subject_id = sa.subject_id
JOIN populations AS p ON p.population_id = cc.population_id;
