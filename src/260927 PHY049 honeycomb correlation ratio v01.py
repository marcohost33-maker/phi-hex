"""PHY049 - W4 v03 honeycomb normalization-free correlation ratio.

Pre-registered contract:
  spec/260927 PHI HEX w4 honeycomb preregistration v03 correlation-ratio.md

This module implements the measurement/preflight layer only. It MUST NOT emit
an external T_BKT claim until the preregistered nonlinear FSS recovery gate is
implemented and green. The fail-closed staging is intentional: measurement
code can be validated without looking at production physics data.
"""
from __future__ import annotations

import importlib.util
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: E402

_SRC = Path(__file__).resolve().parent


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
W4V3_MIN_DEN = 1e-6
W4V3_SPLAY_Z = 2.0
W4V3_SPLAY_MIN_POINTS = 3


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
    run = _nb_run_corr if HAVE_NUMBA else _py_run_corr
    t0 = time.perf_counter()
    data = run(th, nbr, deg, 1.0 / T, n_therm, n_meas,
               seed_for(L, t_idx, s), L)
    return {
        "L": L, "t_idx": t_idx, "T": T, "s": s,
        "g_quarter": float(data[:, 0].mean()),
        "g_half": float(data[:, 1].mean()),
        "wall_s": time.perf_counter() - t0,
    }


def produce(ladder=W4V3_LADDER, t_grid=W4V3_T_GRID,
            n_seeds=W4V3_N_SEEDS, n_therm=W4V3_N_THERM,
            n_meas=W4V3_N_MEAS, max_workers=W4V3_MAX_WORKERS) -> dict:
    """Produce raw seed-level means. No T_BKT estimator is run here."""
    for L in ladder:
        _validate_L(L)
    jobs = [(L, k, T, s, n_therm, n_meas)
            for L in ladder for k, T in enumerate(t_grid)
            for s in range(n_seeds)]
    jobs.sort(key=lambda x: (-x[0], x[1], x[3]))
    t0 = time.perf_counter()
    if max_workers <= 1:
        rows = [_job(j) for j in jobs]
    else:
        with ProcessPoolExecutor(max_workers=max_workers) as ex:
            rows = list(ex.map(_job, jobs, chunksize=1))
    return {
        "module": "PHY049_honeycomb_correlation_ratio_v01",
        "spec": ("spec/260927 PHI HEX w4 honeycomb preregistration v03 "
                 "correlation-ratio.md"),
        "ladder": list(ladder), "t_grid": list(t_grid),
        "n_seeds": n_seeds, "n_therm": n_therm, "n_meas": n_meas,
        "max_workers": max_workers, "wall_s": time.perf_counter() - t0,
        "rows": rows,
    }


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
            rr = [r for r in prod["rows"] if r["L"] == L and r["t_idx"] == k]
            if len(rr) != ns:
                vals.append(None); sems.append(None); continue
            rr.sort(key=lambda x: x["s"])
            q = np.array([r["g_quarter"] for r in rr])
            h = np.array([r["g_half"] for r in rr])
            v, e = jackknife_ratio(q, h)
            vals.append(v); sems.append(e)
        curves[str(L)] = {"R": vals, "SE": sems}
    return {"ladder": ladder, "t_grid": t_grid, "curves": curves}


def persistent_splay(t_grid, r1, e1, r2, e2,
                     z=W4V3_SPLAY_Z, min_points=W4V3_SPLAY_MIN_POINTS):
    """First persistent significant separation; sign learned from top-T tail.

    Requiring a stable tail sign avoids hard-coding lattice-ordering direction.
    Any missing value fails closed.
    """
    n = len(t_grid)
    if n < min_points:
        return None
    tail = []
    for j in range(n - min_points, n):
        if None in (r1[j], e1[j], r2[j], e2[j]):
            return None
        d = r2[j] - r1[j]
        sig = math.hypot(e1[j], e2[j])
        if sig <= 0 or abs(d) <= z * sig:
            return None
        tail.append(math.copysign(1.0, d))
    if not all(x == tail[0] for x in tail):
        return None
    sign = tail[0]
    for k in range(n - min_points + 1):
        ok = True
        for j in range(k, n):
            if None in (r1[j], e1[j], r2[j], e2[j]):
                ok = False; break
            d = sign * (r2[j] - r1[j])
            sig = math.hypot(e1[j], e2[j])
            if sig <= 0 or d <= z * sig:
                ok = False; break
        if ok:
            return float(t_grid[k])
    return None


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
        geometry[str(L)] = {
            "n_ok": lat.n == 2 * L * L,
            "R_aligned": _py_corr(th, L, L // 2) / _py_corr(th, L, L // 4),
        }
    g4 = _p50.g4_synthetic_recovery()
    gates = {
        "G1_geometry": all(v["n_ok"] for v in geometry.values()),
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

if __name__ == "__main__":
    import json
    print(json.dumps(preflight(), indent=2))
