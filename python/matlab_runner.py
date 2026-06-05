"""Shared helpers for batch-running the MATLAB synthesis examples (§4.4).

Reused by the experiment notebooks in ../scripts/ (e.g. run_matlab_timing.ipynb).
Keeping the orchestration here lets the notebooks stay focused on showing each
step's output while the reusable logic lives in one place.

The example scripts print `__TIMING__,<script>,<phase>,<seconds>` markers around
their two solve calls (design = reach_avoid_controller, sop_solve =
solvesop_bounded_control); `parse_timings` extracts them.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import time
from pathlib import Path

DEFAULT_EXAMPLES = [
    "example_dubins_car",
    "example_manipulator",
]

TIMING_RE = re.compile(r"__TIMING__,([^,]+),([^,]+),([0-9.eE+-]+)")


def repo_root() -> Path:
    """Repository root, resolved from this module's location (python/)."""
    return Path(__file__).resolve().parent.parent


def find_matlab(matlab: str = "matlab") -> str | None:
    """Resolve a MATLAB executable name/path to an absolute path, or None."""
    return shutil.which(matlab) or (matlab if os.path.exists(matlab) else None)


def build_matlab_command(matlab_dir: Path | str, example: str) -> str:
    """`matlab -batch` command: put matlab/ on the path and run the example.

    Wrapped in try/catch so a failure still flushes any timing markers already
    printed and reports the error message before exiting non-zero.
    """
    return (
        f"addpath('{matlab_dir}'); "
        f"try, {example}; "
        f"catch e, fprintf(2, '__ERROR__ %s\\n', e.message); exit(1); end; "
        f"exit(0);"
    )


def parse_timings(stdout: str) -> dict[str, float]:
    """Collect the last value seen for each phase from __TIMING__ markers."""
    phases: dict[str, float] = {}
    for m in TIMING_RE.finditer(stdout):
        phases[m.group(2)] = float(m.group(3))
    return phases


def measure_startup_baseline(matlab: str, timeout: float = 120.0) -> float | None:
    """Wall-clock of a trivial `matlab -batch` to estimate startup overhead."""
    cmd = [matlab, "-batch", "disp('__BASELINE_OK__'); exit(0);"]
    t0 = time.perf_counter()
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return None
    elapsed = time.perf_counter() - t0
    return elapsed if r.returncode == 0 else None


def run_one(matlab: str, matlab_dir: Path | str, example: str,
            timeout: float = 3600.0) -> dict:
    """Run one example via `matlab -batch`, returning timing + status."""
    cmd = [matlab, "-batch", build_matlab_command(matlab_dir, example)]
    t0 = time.perf_counter()
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        wall = time.perf_counter() - t0
        stdout, stderr, rc = r.stdout, r.stderr, r.returncode
    except subprocess.TimeoutExpired as e:
        wall = time.perf_counter() - t0
        stdout = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or "")
        stderr, rc = "TIMEOUT", -1

    phases = parse_timings(stdout)
    design, sos = phases.get("design"), phases.get("sop_solve")
    synth = (design + sos) if (design is not None and sos is not None) else None
    status = "ok" if rc == 0 else ("timeout" if stderr == "TIMEOUT" else "error")
    err_msg = ""
    if status == "error":
        m = re.search(r"__ERROR__ (.*)", stderr + "\n" + stdout)
        err_msg = m.group(1).strip() if m else (
            stderr.strip().splitlines()[-1] if stderr.strip() else "")

    return {
        "example": example,
        "status": status,
        "returncode": rc,
        "design_s": design,
        "sop_solve_s": sos,
        "synth_s": synth,
        "wall_s": wall,
        "stdout": stdout,
        "stderr": stderr,
        "error": err_msg,
    }


def write_csv(rows: list[dict], out_path: Path | str) -> Path:
    """Write timing rows to CSV (drops the bulky stdout/stderr columns)."""
    import csv

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["example", "repeat", "status", "returncode",
              "design_s", "sop_solve_s", "synth_s", "wall_s", "wall_net_s", "error"]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})
    return out_path
