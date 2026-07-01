"""Lightweight app-usage logging.

Records one row per Streamlit session (timestamp, client IP, user agent, and the
initial page) into a local SQLite database at ``data/usage.db``. Call
``log_access()`` from the main entry script; it is a no-op after the first run
of a given session.

Behind a reverse proxy (nginx, Apache, a load balancer, corporate gateway) the
direct socket peer is the proxy, so the real client IP is read from the
``X-Forwarded-For`` header when present. Configure your proxy to set that header.

Query the collected data with any SQLite client, e.g.::

    sqlite3 data/usage.db "SELECT ip, COUNT(*) FROM access_log GROUP BY ip;"
"""

from __future__ import annotations

import sqlite3
import threading
from datetime import datetime, timezone
from pathlib import Path

import streamlit as st

try:  # session id is best-effort; never fail logging because of it
    from streamlit.runtime.scriptrunner import get_script_run_ctx
except Exception:  # pragma: no cover - defensive
    get_script_run_ctx = None  # type: ignore[assignment]

_DB_PATH = Path(__file__).resolve().parent.parent / "data" / "usage.db"
_LOCK = threading.Lock()
_INITIALISED = False


def _init_db() -> None:
    global _INITIALISED
    if _INITIALISED:
        return
    _DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(_DB_PATH) as con:
        con.execute("PRAGMA journal_mode=WAL;")
        con.execute(
            """
            CREATE TABLE IF NOT EXISTS access_log (
                id             INTEGER PRIMARY KEY AUTOINCREMENT,
                ts_utc         TEXT NOT NULL,
                session_id     TEXT,
                client_ip      TEXT,
                forwarded_for  TEXT,
                user_agent     TEXT,
                page           TEXT
            )
            """
        )
        con.execute(
            "CREATE INDEX IF NOT EXISTS idx_access_ts ON access_log (ts_utc)"
        )
    _INITIALISED = True


def _session_id() -> str | None:
    if get_script_run_ctx is None:
        return None
    try:
        ctx = get_script_run_ctx()
        return ctx.session_id if ctx is not None else None
    except Exception:
        return None


def _client_ip() -> tuple[str | None, str | None]:
    """Return ``(best_guess_ip, raw_x_forwarded_for)``."""
    headers = {}
    try:
        headers = dict(st.context.headers or {})
    except Exception:
        headers = {}
    # Header lookup is case-insensitive in Streamlit, but be defensive.
    xff = headers.get("X-Forwarded-For") or headers.get("x-forwarded-for")
    if xff:
        # First entry is the original client; the rest are proxies.
        return xff.split(",")[0].strip(), xff
    try:
        return st.context.ip_address, None
    except Exception:
        return None, None


def log_access(page: str | None = None, *, once_per_session: bool = True) -> None:
    """Record an access event.

    Parameters
    ----------
    page:
        Optional label for the page/view being accessed.
    once_per_session:
        When ``True`` (default) only the first call within a browser session is
        written, giving a per-session access count. Set ``False`` to log every
        call (e.g. per page navigation).
    """
    if once_per_session and st.session_state.get("_usage_logged"):
        return

    try:
        _init_db()
        ip, xff = _client_ip()
        headers = {}
        try:
            headers = dict(st.context.headers or {})
        except Exception:
            headers = {}
        user_agent = headers.get("User-Agent") or headers.get("user-agent")
        row = (
            datetime.now(timezone.utc).isoformat(timespec="seconds"),
            _session_id(),
            ip,
            xff,
            user_agent,
            page,
        )
        with _LOCK, sqlite3.connect(_DB_PATH) as con:
            con.execute(
                "INSERT INTO access_log "
                "(ts_utc, session_id, client_ip, forwarded_for, user_agent, page) "
                "VALUES (?, ?, ?, ?, ?, ?)",
                row,
            )
        st.session_state["_usage_logged"] = True
    except Exception:
        # Usage logging must never break the app.
        pass


def fetch_access_log():
    """Return the full access log as a pandas ``DataFrame`` (newest first).

    ``ts_utc`` is parsed to timezone-aware datetimes. Returns an empty frame
    with the expected columns if no data has been recorded yet.
    """
    import pandas as pd

    columns = [
        "id", "ts_utc", "session_id", "client_ip",
        "forwarded_for", "user_agent", "page",
    ]
    try:
        _init_db()
        with sqlite3.connect(_DB_PATH) as con:
            df = pd.read_sql_query(
                "SELECT * FROM access_log ORDER BY ts_utc DESC", con
            )
        if not df.empty:
            df["ts_utc"] = pd.to_datetime(df["ts_utc"], utc=True, errors="coerce")
        return df
    except Exception:
        return pd.DataFrame(columns=columns)

