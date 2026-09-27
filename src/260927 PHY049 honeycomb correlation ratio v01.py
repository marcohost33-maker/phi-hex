"""PHY049 - W4 v03 honeycomb normalization-free correlation ratio.

Pre-registered contract:
  spec/260927 PHI HEX w4 honeycomb preregistration v03 correlation-ratio.md

This module implements the measurement/preflight layer only. G4 synthetic FSS
recovery is deterministic and must pass at runtime before production, but a
repository-level validation claim additionally requires a committed gate log.
It MUST NOT emit an external T_BKT claim until real production evidence passes
G0/G5/G6. The fail-closed staging is intentional: measurement code can be
validated without looking at production physics data.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import numpy as np  # noqa: E402

_SRC = Path(__file__).resolve().parent
_ROOT = _SRC.parent
_RUNTIME_SOURCE_PATHS = (
    "src/260926 PHY045 helicity normalization O1 test v01.py",
    "src/260927 PHY049 honeycomb correlation ratio v01.py",
    "src/260927 PHY050 correlation ratio fss recovery v01.py",
    "spec/260927 PHI HEX w4 honeycomb preregistration v03 correlation-ratio.md",
    "spec/260927 PHI HEX w4 v03a fss estimator hardening.md",
)


def _load(name: str, filename: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _SRC / filename)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


_p45 = _load("phy045_normalization_for_phy049",
             "260926 PHY045 helicity normalization O1 test v01.py")

_p50 = _load("phy050_corr_ratio_fss_for_phy049",
             "260927 PHY050 correlation ratio fss recovery v01.py")

HAVE_NUMBA = _p45.HAVE_NUMBA
_py_wolff_sweep = _p45._py_wolff_sweep
_nb_wolff_sweep = getattr(_p45, "_nb_wolff_sweep", None)
W4V3_LADDER = (48, 72, 96, 144, 192)
W4V3_T_GRID = tuple(round(0.540 + 0.0025 * k, 4) for k in range(29))
W4V3_N_SEEDS = 12
W4V3_N_THERM = 1000
W4V3_N_MEAS = 4000
W4V3_SEED_BASE = 49_000_000
W4V3_MAX_WORKERS = 4
W4V3_WALL_BUDGET_H = 24.0
W4V3_RESULT_PATH = (_SRC.parent / "results" /
                    "260927 PHY049 honeycomb correlation-ratio production.json")
W4V3_MIN_DEN = 1e-6
W4V3_SPLAY_Z = 2.0
W4V3_SPLAY_MIN_POINTS = 3
# For ordered pairs L1<L2 above T_BKT, finite xi implies R_L2 < R_L1.
W4V3_SPLAY_EXPECTED_SIGN = -1.0
# Cached per interpreter/process. Numba is never production-selected solely
# because it is installed: the same-seed Python/Numba trajectory must first
# satisfy VAL-BIT in that process.
_VAL_BIT_OK: bool | None = None


def seed_for(L: int, t_idx: int, s: int) -> int:
    return W4V3_SEED_BASE + 1000 * L + 100 * t_idx + s


def _validate_L(L: int) -> None:
    if L <= 0 or L % 4:
        raise ValueError("PHY049 requires positive L divisible by 4")


def _py_corr(th: np.ndarray, L: int, r: int) -> float:
    """Same-sublattice correlation, averaged over a1/a2 and A/B."""
    if r <= 0 or r >= L:
        raise ValueError("r must satisfy 0 < r < L")
    acc = 0.0
    n_terms = 0
    for i in range(L):
        for j in range(L):
            for s in range(2):
                a = 2 * (i * L + j) + s
                b1 = 2 * (((i + r) % L) * L + j) + s
                b2 = 2 * (i * L + ((j + r) % L)) + s
                acc += math.cos(th[a] - th[b1])
                acc += math.cos(th[a] - th[b2])
                n_terms += 2
    return acc / n_terms


if HAVE_NUMBA:
    import numba

    _nb_corr = numba.njit(cache=False, fastmath=False)(_py_corr)

    @numba.njit(cache=False, fastmath=False)
    def _nb_run_corr(th, nbr, deg, beta, n_therm, n_meas, seed, L):
        np.random.seed(seed)
        n = th.shape[0]
        stack = np.zeros(2 * n, dtype=np.int64)
        incl = np.zeros(n, dtype=np.int64)
        for _ in range(n_therm):
            _nb_wolff_sweep(th, nbr, deg, beta, n, stack, incl)
        out = np.zeros((n_meas, 2))
        rq = L // 4
        rh = L // 2
        for k in range(n_meas):
            _nb_wolff_sweep(th, nbr, deg, beta, n, stack, incl)
            out[k, 0] = _nb_corr(th, L, rq)
            out[k, 1] = _nb_corr(th, L, rh)
        return out


def _py_run_corr(th, nbr, deg, beta, n_therm, n_meas, seed, L):
    np.random.seed(seed)
    n = th.shape[0]
    stack = np.zeros(2 * n, dtype=np.int64)
    incl = np.zeros(n, dtype=np.int64)
    for _ in range(n_therm):
        _py_wolff_sweep(th, nbr, deg, beta, n, stack, incl)
    out = np.zeros((n_meas, 2))
    for k in range(n_meas):
        _py_wolff_sweep(th, nbr, deg, beta, n, stack, incl)
        out[k, 0] = _py_corr(th, L, L // 4)
        out[k, 1] = _py_corr(th, L, L // 2)
    return out


def _job(args: tuple) -> dict:
    L, t_idx, T, s, n_therm, n_meas = args
    _validate_L(L)
    lat = _p45.build("honeycomb", L)
    if lat.n != 2 * L * L:
        raise RuntimeError("unexpected honeycomb indexing")
    nbr, deg = _p45._nbr_arrays(lat)
    th = np.zeros(lat.n, dtype=float)
    run = _production_run_backend()
    t0 = time.perf_counter()
    data = run(th, nbr, deg, 1.0 / T, n_therm, n_meas,
               seed_for(L, t_idx, s), L)
    return {
        "L": L, "t_idx": t_idx, "T": T, "s": s,
        "seed": seed_for(L, t_idx, s),
        "g_quarter": float(data[:, 0].mean()),
        "g_half": float(data[:, 1].mean()),
        "wall_s": time.perf_counter() - t0,
    }


def _campaign_contract(
    ladder: tuple[int, ...],
    t_grid: tuple[float, ...],
    n_seeds: int,
    n_therm: int,
    n_meas: int,
    max_workers: int,
    wall_budget_h: float,
) -> dict:
    """Immutable execution contract persisted in every checkpoint."""
    return {
        "ladder": list(ladder),
        "t_grid": list(t_grid),
        "n_seeds": n_seeds,
        "n_therm": n_therm,
        "n_meas": n_meas,
        "max_workers": max_workers,
        "wall_budget_h": wall_budget_h,
    }


def _atomic_write_json(path: Path, payload: dict) -> None:
    """Crash-durably replace a checkpoint; never expose partial JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, path)
    if os.name == "posix":
        dir_fd = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(dir_fd)
        finally:
            os.close(dir_fd)


def _strict_json_equal(value: object, expected: object) -> bool:
    """Compare JSON semantics without bool/int or int/float coercion."""
    try:
        lhs = json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        )
        rhs = json.dumps(
            expected, sort_keys=True, separators=(",", ":"), allow_nan=False
        )
    except (TypeError, ValueError):
        return False
    return lhs == rhs


def _json_number(value: object) -> bool:
    return type(value) in (int, float) and math.isfinite(float(value))


def _file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest().upper()


def _runtime_provenance() -> dict:
    """Immutable software/runtime fingerprint for one production campaign."""
    source_sha256 = {
        rel: _file_sha256(_ROOT / rel) for rel in _RUNTIME_SOURCE_PATHS
    }
    numba_version = None
    if HAVE_NUMBA:
        numba_version = str(getattr(_p45.numba, "__version__", "unknown"))
    return {
        "schema": "PHY049_RUNTIME_PROVENANCE_V1",
        "python_implementation": platform.python_implementation(),
        "python_version": platform.python_version(),
        "numpy_version": str(np.__version__),
        "numba_version": numba_version,
        "production_backend": "numba" if HAVE_NUMBA else "python",
        "sys_platform": sys.platform,
        "machine": platform.machine() or "unknown",
        "thread_env": {
            "OPENBLAS_NUM_THREADS": os.environ.get("OPENBLAS_NUM_THREADS"),
            "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS"),
            "MKL_NUM_THREADS": os.environ.get("MKL_NUM_THREADS"),
        },
        "source_sha256": source_sha256,
    }


def _completed_block_valid(
    rows: list[dict],
    L: int,
    t_grid: tuple[float, ...],
    n_seeds: int,
) -> bool:
    """Resume only an exact, finite, complete L block; partial blocks are toxic."""
    rr = [r for r in rows if isinstance(r, dict) and r.get("L") == L]
    expected_n = len(t_grid) * n_seeds
    if len(rr) != expected_n:
        return False
    seen = set()
    for r in rr:
        try:
            if type(r.get("L")) is not int or r["L"] != L:
                return False
            k = r["t_idx"]
            seed_idx = r["s"]
            if type(k) is not int or type(seed_idx) is not int:
                return False
            if not (0 <= k < len(t_grid) and 0 <= seed_idx < n_seeds):
                return False
            if (k, seed_idx) in seen:
                return False
            seen.add((k, seed_idx))
            if (
                type(r.get("seed")) is not int
                or r["seed"] != seed_for(L, k, seed_idx)
            ):
                return False
            if type(r.get("T")) not in (int, float):
                return False
            if float(r["T"]) != t_grid[k]:
                return False
            for key in ("g_quarter", "g_half", "wall_s"):
                if type(r.get(key)) not in (int, float):
                    return False
                if not math.isfinite(float(r[key])):
                    return False
            if not (
                -1.0 <= float(r["g_quarter"]) <= 1.0
                and -1.0 <= float(r["g_half"]) <= 1.0
                and float(r["wall_s"]) >= 0.0
            ):
                return False
        except (KeyError, TypeError, ValueError, OverflowError):
            return False
    return len(seen) == expected_n


def _checkpoint_payload(
    *,
    pre: dict,
    contract: dict,
    runtime_provenance: dict,
    rows: list[dict],
    unmeasured: list[dict],
    wall_s: float,
    status: str,
    inflight_L: int | None = None,
    inflight_started_epoch_s: float | None = None,
) -> dict:
    expected = (
        len(contract["ladder"])
        * len(contract["t_grid"])
        * contract["n_seeds"]
    )
    complete = len(rows) == expected and not unmeasured
    return {
        "module": "PHY049_honeycomb_correlation_ratio_v01",
        "spec": ("spec/260927 PHI HEX w4 honeycomb preregistration v03 "
                 "correlation-ratio.md"),
        "preflight_gates": dict(pre["gates"]),
        "runtime_provenance": runtime_provenance,
        **contract,
        "campaign_contract": dict(contract),
        "checkpoint_status": "COMPLETE" if complete else status,
        "wall_s": wall_s,
        "inflight_L": None if complete else inflight_L,
        "inflight_started_epoch_s": (
            None if complete else inflight_started_epoch_s
        ),
        "rows": rows,
        "unmeasured": unmeasured,
        "complete": complete,
    }


def produce(
    ladder=W4V3_LADDER,
    t_grid=W4V3_T_GRID,
    n_seeds=W4V3_N_SEEDS,
    n_therm=W4V3_N_THERM,
    n_meas=W4V3_N_MEAS,
    max_workers=W4V3_MAX_WORKERS,
    wall_budget_h=W4V3_WALL_BUDGET_H,
    *,
    checkpoint_path: str | Path | None = None,
    resume: bool = False,
) -> dict:
    """Produce raw seed means with fail-closed, whole-L transactional resume.

    Public production is locked behind VAL-BIT + G1-G4. Checkpoints are
    committed only after a complete lattice-size block; an interrupted
    in-flight L is recomputed rather than accepting partial evidence.
    """
    pre = preflight()
    if not pre["production_measurement_eligible"]:
        failed = [k for k, ok in pre["gates"].items() if not ok]
        raise RuntimeError(
            "PHY049 preflight failed; production is locked: "
            + ", ".join(failed)
        )

    try:
        raw_ladder = tuple(ladder)
        raw_t_grid = tuple(t_grid)
    except TypeError as exc:
        raise ValueError("ladder and t_grid must be finite sequences") from exc
    if (
        not raw_ladder
        or any(type(L) is not int for L in raw_ladder)
        or len(set(raw_ladder)) != len(raw_ladder)
    ):
        raise ValueError("ladder must contain unique literal integers")
    if (
        not raw_t_grid
        or any(type(T) not in (int, float) for T in raw_t_grid)
        or any(not math.isfinite(float(T)) or float(T) <= 0.0 for T in raw_t_grid)
        or len(set(float(T) for T in raw_t_grid)) != len(raw_t_grid)
    ):
        raise ValueError("t_grid must contain unique finite positive numbers")
    ladder = raw_ladder
    t_grid = tuple(float(T) for T in raw_t_grid)
    for L in ladder:
        _validate_L(L)
    if (
        type(n_seeds) is not int
        or type(n_therm) is not int
        or type(n_meas) is not int
        or type(max_workers) is not int
        or min(n_seeds, n_therm, n_meas, max_workers) <= 0
    ):
        raise ValueError("seed/sweep/worker counts must be positive integers")
    if (
        type(wall_budget_h) not in (int, float)
        or not math.isfinite(float(wall_budget_h))
        or float(wall_budget_h) < 0.0
    ):
        raise ValueError("wall_budget_h must be a finite non-negative number")

    wall_budget_h = float(wall_budget_h)
    runtime_provenance = _runtime_provenance()
    contract = _campaign_contract(
        ladder, t_grid, n_seeds, n_therm, n_meas, max_workers, wall_budget_h
    )
    checkpoint = None if checkpoint_path is None else Path(checkpoint_path)
    rows: list[dict] = []
    unmeasured: list[dict] = []
    prior_wall_s = 0.0
    completed: set[int] = set()

    if resume:
        if checkpoint is None or not checkpoint.exists():
            raise ValueError("resume requires an existing checkpoint_path")
        with checkpoint.open("r", encoding="utf-8") as handle:
            saved = json.load(handle)
        if not isinstance(saved, dict):
            raise RuntimeError("checkpoint root must be a JSON object")
        if saved.get("module") != "PHY049_honeycomb_correlation_ratio_v01":
            raise RuntimeError("checkpoint module identity is invalid")
        if saved.get("spec") != (
            "spec/260927 PHI HEX w4 honeycomb preregistration v03 "
            "correlation-ratio.md"
        ):
            raise RuntimeError("checkpoint spec identity is invalid")
        if not _strict_json_equal(saved.get("campaign_contract"), contract):
            raise RuntimeError("checkpoint contract does not match this campaign")
        if not _strict_json_equal(saved.get("preflight_gates"), pre["gates"]):
            raise RuntimeError("checkpoint preflight evidence does not match")
        if not _strict_json_equal(
            saved.get("runtime_provenance"), runtime_provenance
        ):
            raise RuntimeError(
                "checkpoint runtime/source provenance does not match"
            )

        status = saved.get("checkpoint_status")
        if status == "WALL_BUDGET_STOP":
            raise RuntimeError(
                "terminal WALL_BUDGET_STOP checkpoint cannot be resumed"
            )
        if status not in {"IN_PROGRESS", "BLOCK_IN_PROGRESS", "COMPLETE"}:
            raise RuntimeError("checkpoint_status is invalid")

        raw_unmeasured = saved.get("unmeasured")
        if not isinstance(raw_unmeasured, list):
            raise RuntimeError("checkpoint unmeasured must be a list")
        raw_rows = saved.get("rows")
        if not isinstance(raw_rows, list) or any(
            not isinstance(r, dict) for r in raw_rows
        ):
            raise RuntimeError("checkpoint rows must be a list of objects")

        wall_value = saved.get("wall_s")
        if not _json_number(wall_value) or float(wall_value) < 0.0:
            raise RuntimeError("checkpoint wall_s is invalid")
        prior_wall_s = float(wall_value)

        rows = list(raw_rows)
        if any(type(r.get("L")) is not int or r["L"] not in ladder for r in rows):
            raise RuntimeError("checkpoint contains foreign lattice rows")
        for L in ladder:
            count = sum(1 for r in rows if r.get("L") == L)
            if count == 0:
                continue
            if not _completed_block_valid(rows, L, t_grid, n_seeds):
                raise RuntimeError(f"checkpoint has partial/invalid L={L} block")
            completed.add(L)

        if status == "COMPLETE":
            if saved.get("complete") is not True or raw_unmeasured != []:
                raise RuntimeError("COMPLETE checkpoint envelope is incoherent")
            if (
                saved.get("inflight_L") is not None
                or saved.get("inflight_started_epoch_s") is not None
            ):
                raise RuntimeError("COMPLETE checkpoint cannot have in-flight state")
            if completed != set(ladder):
                raise RuntimeError("checkpoint claims complete with missing L")
            return saved

        if saved.get("complete") is not False or raw_unmeasured != []:
            raise RuntimeError("nonterminal checkpoint envelope is incoherent")

        if status == "IN_PROGRESS":
            if (
                saved.get("inflight_L") is not None
                or saved.get("inflight_started_epoch_s") is not None
            ):
                raise RuntimeError("IN_PROGRESS checkpoint has stale in-flight state")
        else:
            inflight_L = saved.get("inflight_L")
            started_epoch_s = saved.get("inflight_started_epoch_s")
            if (
                type(inflight_L) is not int
                or inflight_L not in ladder
                or inflight_L in completed
            ):
                raise RuntimeError("BLOCK_IN_PROGRESS lattice identity is invalid")
            if not _json_number(started_epoch_s) or float(started_epoch_s) < 0.0:
                raise RuntimeError("BLOCK_IN_PROGRESS start time is invalid")
            now_epoch_s = time.time()
            if (
                not math.isfinite(now_epoch_s)
                or now_epoch_s < float(started_epoch_s)
            ):
                raise RuntimeError("wall clock moved backwards across resume")
            prior_wall_s += now_epoch_s - float(started_epoch_s)
            if not math.isfinite(prior_wall_s):
                raise RuntimeError("checkpoint accumulated wall time is invalid")
    elif checkpoint is not None and checkpoint.exists():
        raise FileExistsError(
            "checkpoint exists; pass resume=True or choose a new path"
        )

    t0 = time.perf_counter()
    budget_s = 3600.0 * wall_budget_h

    for L in ladder:
        if L in completed:
            continue
        jobs = [
            (L, k, T, s, n_therm, n_meas)
            for k, T in enumerate(t_grid)
            for s in range(n_seeds)
        ]
        elapsed = prior_wall_s + (time.perf_counter() - t0)
        if elapsed >= budget_s:
            stop_index = ladder.index(L)
            for pending_L in ladder[stop_index:]:
                if pending_L in completed:
                    continue
                for k, T in enumerate(t_grid):
                    for seed_idx in range(n_seeds):
                        unmeasured.append({
                            "L": pending_L,
                            "t_idx": k,
                            "T": T,
                            "s": seed_idx,
                            "seed": seed_for(pending_L, k, seed_idx),
                            "reason": "WALL_BUDGET_STOP",
                        })
            break

        if checkpoint is not None:
            started_epoch_s = time.time()
            if not math.isfinite(started_epoch_s) or started_epoch_s < 0.0:
                raise RuntimeError("wall-clock anchor is invalid")
            payload = _checkpoint_payload(
                pre=pre,
                contract=contract,
                runtime_provenance=runtime_provenance,
                rows=rows,
                unmeasured=[],
                wall_s=elapsed,
                status="BLOCK_IN_PROGRESS",
                inflight_L=L,
                inflight_started_epoch_s=started_epoch_s,
            )
            _atomic_write_json(checkpoint, payload)

        if max_workers <= 1:
            block = [_job(j) for j in jobs]
        else:
            with ProcessPoolExecutor(max_workers=max_workers) as ex:
                block = list(ex.map(_job, jobs, chunksize=1))
        if not _completed_block_valid(block, L, t_grid, n_seeds):
            raise RuntimeError(f"internal production block validation failed L={L}")
        rows.extend(block)
        completed.add(L)

        if checkpoint is not None:
            elapsed = prior_wall_s + (time.perf_counter() - t0)
            payload = _checkpoint_payload(
                pre=pre,
                contract=contract,
                runtime_provenance=runtime_provenance,
                rows=rows,
                unmeasured=[],
                wall_s=elapsed,
                status="IN_PROGRESS",
            )
            _atomic_write_json(checkpoint, payload)

    elapsed = prior_wall_s + (time.perf_counter() - t0)
    status = "WALL_BUDGET_STOP" if unmeasured else "IN_PROGRESS"
    product = _checkpoint_payload(
        pre=pre,
        contract=contract,
        runtime_provenance=runtime_provenance,
        rows=rows,
        unmeasured=unmeasured,
        wall_s=elapsed,
        status=status,
    )
    if checkpoint is not None:
        _atomic_write_json(checkpoint, product)
    return product

def ratio_of_means(g_quarter: np.ndarray, g_half: np.ndarray) -> float | None:
    """R=<g(L/2)>/<g(L/4)>; never mean of per-seed ratios."""
    den = float(np.mean(g_quarter))
    if not math.isfinite(den) or abs(den) <= W4V3_MIN_DEN:
        return None
    num = float(np.mean(g_half))
    if not math.isfinite(num):
        return None
    return num / den


def jackknife_ratio(g_quarter: np.ndarray, g_half: np.ndarray) -> tuple[float | None, float | None]:
    """Ratio-of-means and delete-one-seed jackknife SEM."""
    q = np.asarray(g_quarter, dtype=float)
    h = np.asarray(g_half, dtype=float)
    if q.shape != h.shape or q.ndim != 1 or len(q) < 2:
        return None, None
    full = ratio_of_means(q, h)
    if full is None:
        return None, None
    vals = []
    for s in range(len(q)):
        v = ratio_of_means(np.delete(q, s), np.delete(h, s))
        if v is None:
            return full, None
        vals.append(v)
    vals = np.asarray(vals)
    m = float(vals.mean())
    se = math.sqrt((len(vals) - 1) / len(vals) * float(np.sum((vals - m) ** 2)))
    return full, se


def aggregate(prod: dict) -> dict:
    """Build R_L(T) curves with delete-one-seed errors; fail closed."""
    ladder = [int(x) for x in prod["ladder"]]
    t_grid = [float(x) for x in prod["t_grid"]]
    ns = int(prod["n_seeds"])
    curves = {}
    for L in ladder:
        vals, sems = [], []
        for k, _T in enumerate(t_grid):
            rr = [r for r in prod["rows"] if r.get("L") == L and r.get("t_idx") == k]
            rr.sort(key=lambda x: int(x.get("s", -1)))
            expected_s = list(range(ns))
            expected_T = t_grid[k]
            valid_rows = (
                len(rr) == ns
                and [int(r.get("s", -1)) for r in rr] == expected_s
                and all(float(r.get("T", math.nan)) == expected_T for r in rr)
                and all(
                    int(r.get("seed", -1)) == seed_for(L, k, int(r.get("s", -1)))
                    for r in rr
                )
            )
            if not valid_rows:
                vals.append(None)
                sems.append(None)
                continue
            q = np.array([r["g_quarter"] for r in rr], dtype=float)
            h = np.array([r["g_half"] for r in rr], dtype=float)
            if np.any(~np.isfinite(q)) or np.any(~np.isfinite(h)):
                vals.append(None)
                sems.append(None)
                continue
            v, e = jackknife_ratio(q, h)
            vals.append(v)
            sems.append(e)
        curves[str(L)] = {"R": vals, "SE": sems}
    return {"ladder": ladder, "t_grid": t_grid, "curves": curves}


def persistent_splay(t_grid, r1, e1, r2, e2,
                     z=W4V3_SPLAY_Z, min_points=W4V3_SPLAY_MIN_POINTS,
                     expected_sign=W4V3_SPLAY_EXPECTED_SIGN):
    """First persistent high-T separation for an ordered pair L1<L2.

    The preregistered sign is fixed before seeing production data. For the
    high-temperature phase g(r) decays with distance, so at equal fractional
    distances R_L2 < R_L1 and D=R_L2-R_L1 is negative.
    Missing, non-finite, zero-uncertainty, or opposite-sign points break
    persistence (fail closed).
    """
    n = len(t_grid)
    if n < min_points or expected_sign not in (-1.0, 1.0):
        return None
    for k in range(n - min_points + 1):
        ok = True
        for j in range(k, n):
            vals = (r1[j], e1[j], r2[j], e2[j])
            if any(v is None or not math.isfinite(float(v)) for v in vals):
                ok = False
                break
            d = expected_sign * (float(r2[j]) - float(r1[j]))
            sig = math.hypot(float(e1[j]), float(e2[j]))
            if not math.isfinite(sig) or sig <= 0.0 or d <= z * sig:
                ok = False
                break
        if ok:
            return float(t_grid[k])
    return None


def _geometry_invariants(L: int) -> dict:
    """Directly test the index identities required by preregistered G1."""
    _validate_L(L)
    wrap_a1 = True
    wrap_a2 = True
    same_sublattice = True
    for i in range(L):
        for j in range(L):
            for s in range(2):
                a = 2 * (i * L + j) + s
                a1_wrap = 2 * ((((i + L) % L) * L + j)) + s
                a2_wrap = 2 * ((i * L + ((j + L) % L))) + s
                wrap_a1 = wrap_a1 and a1_wrap == a
                wrap_a2 = wrap_a2 and a2_wrap == a
                for r in (L // 4, L // 2):
                    b1 = 2 * (((i + r) % L) * L + j) + s
                    b2 = 2 * (i * L + ((j + r) % L)) + s
                    same_sublattice = (
                        same_sublattice and b1 % 2 == s and b2 % 2 == s
                    )
    return {
        "translation_L_a1": wrap_a1,
        "translation_L_a2": wrap_a2,
        "same_sublattice_q_h": same_sublattice,
    }


def _backend_bit_identity() -> bool:
    """VAL-BIT: identical tiny trajectory before Numba is production-eligible."""
    if not HAVE_NUMBA:
        return True
    L = 8
    lat = _p45.build("honeycomb", L)
    nbr, deg = _p45._nbr_arrays(lat)
    args = (nbr, deg, 1.0 / 0.57, 2, 3, seed_for(L, 0, 0), L)
    py = _py_run_corr(np.zeros(lat.n), *args)
    nb = _nb_run_corr(np.zeros(lat.n), *args)
    return py.dtype == nb.dtype and py.shape == nb.shape and np.array_equal(py, nb)


def _production_run_backend():
    """Return the production kernel only after fail-closed VAL-BIT approval.

    The check is cached per interpreter. ProcessPool workers therefore verify
    their own runtime before their first Numba job; a direct _job() call cannot
    bypass VAL-BIT either.
    """
    global _VAL_BIT_OK
    if not HAVE_NUMBA:
        return _py_run_corr
    if _nb_run_corr is None:
        raise RuntimeError("Numba reported available but kernel is missing")
    if _VAL_BIT_OK is None:
        _VAL_BIT_OK = bool(_backend_bit_identity())
    if not _VAL_BIT_OK:
        raise RuntimeError(
            "VAL-BIT failed: Numba backend is not production-eligible"
        )
    return _nb_run_corr


def preflight() -> dict:
    """Cheap pre-production gates. Physics interpretation remains disabled."""
    seed_set = {seed_for(L, k, s) for L in W4V3_LADDER
                for k in range(len(W4V3_T_GRID))
                for s in range(W4V3_N_SEEDS)}
    expected = len(W4V3_LADDER) * len(W4V3_T_GRID) * W4V3_N_SEEDS
    geometry = {}
    for L in W4V3_LADDER:
        lat = _p45.build("honeycomb", L)
        th = np.zeros(lat.n)
        inv = _geometry_invariants(L)
        geometry[str(L)] = {
            "n_ok": lat.n == 2 * L * L,
            **inv,
            "R_aligned": _py_corr(th, L, L // 2) / _py_corr(th, L, L // 4),
        }
    val_bit = _backend_bit_identity()
    g4 = _p50.g4_synthetic_recovery()
    gates = {
        "VAL_BIT_numba": val_bit,
        "G1_geometry": all(
            v["n_ok"]
            and v["translation_L_a1"]
            and v["translation_L_a2"]
            and v["same_sublattice_q_h"]
            for v in geometry.values()
        ),
        "G2_aligned_limit": all(
            abs(v["R_aligned"] - 1.0) < 1e-15 for v in geometry.values()
        ),
        "G3_seed_unique": len(seed_set) == expected,
        "G4_fss_recovery": bool(g4["passed"]),
    }
    production_eligible = all(gates.values())
    return {
        "module": "PHY049_honeycomb_correlation_ratio_v01",
        "stage": (
            "PREPRODUCTION_G4_VALIDATED"
            if production_eligible
            else "PREPRODUCTION_GATE_FAILED"
        ),
        "geometry": geometry,
        "gates": gates,
        "g4": g4,
        "production_measurement_eligible": production_eligible,
        "overall_interpretation_enabled": False,
        "claim_ceiling": (
            "NO_PHYSICS_INTERPRETATION until real PHY049 production data "
            "satisfy G0 input, G5 power and G6 robustness with committed "
            "gate evidence"
        ),
    }

def _main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="PHY049 preflight or preregistered production campaign"
    )
    parser.add_argument(
        "--production",
        action="store_true",
        help="run the exact W4-v03 production campaign",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=W4V3_RESULT_PATH,
        help="atomic production checkpoint/result JSON",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="resume an exact whole-L checkpoint; never resumes budget-stop data",
    )
    args = parser.parse_args(argv)
    if not args.production:
        print(json.dumps(preflight(), indent=2, sort_keys=True))
        return 0
    product = produce(
        checkpoint_path=args.output,
        resume=args.resume,
    )
    summary = {
        "output": str(args.output),
        "checkpoint_status": product["checkpoint_status"],
        "complete": product["complete"],
        "rows": len(product["rows"]),
        "unmeasured": len(product["unmeasured"]),
        "wall_s": product["wall_s"],
        "claim_ceiling": (
            "NO_PHYSICS_INTERPRETATION: run PHY050 assess_production() only "
            "after an exact complete PHY049 product exists."
        ),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
