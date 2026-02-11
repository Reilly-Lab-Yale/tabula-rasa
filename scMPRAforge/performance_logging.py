# taskstream_log.py
from __future__ import annotations

import atexit
import json
import threading
import time
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Set, Tuple

from dask.distributed import Client


def _iter_taskstream_rows(payload: Iterable[Dict[str, Any]]) -> Iterable[Dict[str, Any]]:
    """
    Flatten Dask task-stream structures into one row per (key, worker, start, stop).
    """
    for ev in payload:
        key = ev.get("key")
        name = ev.get("name")
        for ss in ev.get("startstops", []) or []:
            yield {
                "key": key,
                "name": name,
                "worker": ss.get("worker"),
                "start": ss.get("start"),
                "stop": ss.get("stop"),
            }


class TaskStreamLogger:
    """
    Periodically polls client.get_task_stream() and appends compact JSONL rows.

    Robustness:
      - Writes incrementally + flushes every poll.
      - If the scheduler/cluster dies mid-run, you still keep everything already written.
      - Stop() does a final poll best-effort.
    """

    def __init__(
        self,
        client: Client,
        out_jsonl: str | Path,
        *,
        interval_s: float = 10.0,
        overlap_s: float = 2.0,
        include_header: bool = True,
    ):
        self.client = client
        self.out = Path(out_jsonl)
        self.out.parent.mkdir(parents=True, exist_ok=True)
        self.interval_s = float(interval_s)
        self.overlap_s = float(overlap_s)
        self._stop = threading.Event()
        self._thr: Optional[threading.Thread] = None
        self._seen: Set[Tuple[Any, Any, Any, Any]] = set()  # (key, worker, start, stop)
        self._last_poll_t = time.time()

        if include_header:
            header = {
                "type": "taskstream_header",
                "ts": time.time(),
                "dashboard_link": getattr(client, "dashboard_link", None),
            }
            with self.out.open("a", buffering=1) as f:
                f.write(json.dumps(header) + "\n")
                f.flush()

        atexit.register(self.stop)

    def start(self) -> "TaskStreamLogger":
        if self._thr is not None:
            return self
        self._thr = threading.Thread(target=self._loop, name="dask-taskstream-logger", daemon=True)
        self._thr.start()
        return self

    def stop(self) -> None:
        self._stop.set()
        if self._thr is not None:
            self._thr.join(timeout=2.0)
        # final best-effort poll
        self._poll_once(best_effort=True)

    def _loop(self) -> None:
        while not self._stop.is_set():
            self._poll_once(best_effort=True)
            self._stop.wait(self.interval_s)

    def _poll_once(self, *, best_effort: bool) -> None:
        now = time.time()
        # overlap to avoid missing tasks that straddle poll boundaries
        start = self._last_poll_t - self.overlap_s
        self._last_poll_t = now

        try:
            events = self.client.get_task_stream(start=start, stop=now)
        except Exception:
            if best_effort:
                return
            raise

        rows = []
        for r in _iter_taskstream_rows(events):
            ident = (r["key"], r["worker"], r["start"], r["stop"])
            if None in ident:
                continue
            if ident in self._seen:
                continue
            self._seen.add(ident)
            r["dur_s"] = float(r["stop"] - r["start"])
            r["ts"] = now  # when we recorded it
            rows.append(r)

        if not rows:
            return

        # append + flush so partial data survives failures
        with self.out.open("a", buffering=1) as f:
            for r in rows:
                f.write(json.dumps(r, sort_keys=True) + "\n")
            f.flush()


def start_taskstream_logger(
    client: Client,
    out_jsonl: str | Path,
    *,
    interval_s: float = 10.0,
    overlap_s: float = 2.0,
) -> TaskStreamLogger:
    """
    Convenience wrapper.
    """
    return TaskStreamLogger(
        client,
        out_jsonl,
        interval_s=interval_s,
        overlap_s=overlap_s,
        include_header=True,
    ).start()