import sqlite3
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DB_PATH = ROOT / "cell_counts.db"


def connect(db_path=DB_PATH):
    db_path = Path(db_path)
    if not db_path.exists():
        raise FileNotFoundError(f"{db_path} not found; run `python load_data.py` first")
    return sqlite3.connect(f"file:{db_path}?mode=ro", uri=True, check_same_thread=False)


def query(conn, sql, params=()):
    return pd.read_sql_query(sql, conn, params=params)
