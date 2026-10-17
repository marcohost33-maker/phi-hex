"""PHY044 — W4 honeycomb calibration before large-L production.

Purpose
-------
Pre-registered calibration for Issue #45.  This module does NOT estimate a
new T_BKT.  It measures whether the existing PHY041 Wang-Landau/1-t sampler
has enough convergence/coverage at a larger honeycomb lattice to justify a
subsequent >=3-walker W4 production ladder, and records a canonical Wolff
cross-check at the same temperatures.

The default L=64, two-WL-walker run is intentionally a CALIBRATION only:
two walkers can estimate cost and expose gross sampler disagreement, but
cannot satisfy the W4 production requirement of >=3 walkers per L.

Contract: spec/260927 PHI HEX w4 preregistration v01.md
Coworker Research / Coworkerz | 2026-09-27
"""
from __future__ import annotations

import importlib.util
import inspect
import json
import math
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

_SRC = Path(__file__).resolve().parent
_ROOT = _SRC.parent


def _load(name: str, filename: str):
    if name in __import__("sys").modules:
        return __import__("sys").modules[name]
    spec = importlib.util.spec_from_file_location(name, _SRC / filename)
    mod = importlib.util.module_from_spec(spec)
    __import__("sys").modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


_p41 = _load(
    "phy041_honeycomb_wl",
    "260702 PHY041 honeycomb wang-landau entropic helicity v01.py",
)

# Dedicated high stream namespace.  wl_entropic_lattice adds +L internally.
_STREAM_CAL_WL_BASE = 2_000_000

DEFAULT_TEMPS = (0.565, 0.575, 0.585)
SPREAD_LIMIT = 0.04
MASS_LIMIT = 1e-3
LNF_FINAL = 1e-5


def _stream_for_calibration(walker: int) -> int:
    if walker < 0:
        raise ValueError("walker muss >=0 sein")
    return _STREAM_CAL_WL_BASE + 1000 * walker


def _prod_sweeps_for(nbins: int) -> int:
    """Match the PHY042 production-density floor without importing PHY042."""
    if nbins <= 0:
        raise ValueError("nbins muss >0 sein")
    return max(30000, 60 * nbins)


def _validate_inputs(L: int, n_wl_walkers: int, temps) -> tuple[float, ...]:
    if isinstance(L, bool) or not isinstance(L, int) or L < 8:
        raise ValueError("L muss eine ganze Zahl >=8 sein")
    if isinstance(n_wl_walkers, bool) or not isinstance(n_wl_walkers, int) or n_wl_walkers < 2:
        raise ValueError("Kalibrierung braucht mindestens 2 WL-Walker (ganze Zahl)")
    ts = tuple(float(t) for t in temps)
    if not ts or any((not math.isfinite(t) or t <= 0.0) for t in ts):
        raise ValueError("temps muessen endlich und >0 sein")
    if tuple(sorted(ts)) != ts or len(set(ts)) != len(ts):
        raise ValueError("temps muessen streng aufsteigend sein")
    return ts


def _git_state() -> tuple[str | None, bool | None]:
    """Commit and dirty flag of the producing tree; (None, None) if unknown."""
    try:
        head = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=_ROOT, capture_output=True,
            text=True, check=True).stdout.strip()
        status = subprocess.run(
            ["git", "status", "--porcelain"], cwd=_ROOT, capture_output=True,
            text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return None, None
    return head, bool(status.strip())


def _provenance(L: int, n_wl_walkers: int, wolff_n_seeds: int) -> dict:
    """Spec §7 reproducibility record: code, environment, lnf_final, streams.

    Stream IDs are derived from the constants the samplers actually use
    (wl_entropic_lattice adds +L; Wolff anchor/reference use
    _STREAM_ANCHOR/_STREAM_REF + s + 1000*L), so the record cannot drift
    from the run (PR #51 Codex round 2026-10-16).
    """
    commit, dirty = _git_state()
    anchor_seeds = inspect.signature(_p41.wolff_anchor).parameters["n_seeds"].default
    return {
        "git_commit": commit,
        "git_dirty": dirty,
        "python_version": sys.version.split()[0],
        "numpy_version": np.__version__,
        "platform": platform.platform(),
        "lnf_final": LNF_FINAL,
        "streams": {
            "wl": [_stream_for_calibration(w) + L for w in range(n_wl_walkers)],
            "wolff_anchor": [_p41._STREAM_ANCHOR + s + 1000 * L
                             for s in range(anchor_seeds)],
            "wolff_reference": [_p41._STREAM_REF + s + 1000 * L
                                for s in range(wolff_n_seeds)],
        },
    }


def _calibration_verdict(*, all_1t: bool, leak_max: float,
                         uncovered_max: float, max_y2_spread: float) -> str:
    """Pure, pre-registered stop/continue gate; never yields a T_BKT claim."""
    if not all_1t:
        return "CALIBRATION_FAIL_1T_NOT_ENGAGED"
    # Each accepting branch is bound to a TRUE comparison, so NaN fails closed
    # (``NaN >= limit`` is False and used to fall through to PASS; 007acc2 R7).
    if not (leak_max < MASS_LIMIT and uncovered_max < MASS_LIMIT):
        return "CALIBRATION_FAIL_COVERAGE"
    if not max_y2_spread < SPREAD_LIMIT:
        return "CALIBRATION_NEEDS_MORE_WL"
    return "CALIBRATION_PASS_FOR_COST_MODEL"


CALIBRATION_L = 64  # the preregistered W4-G1 calibration lattice


def _scale_decision(verdict: str, L: int) -> str:
    """Production planning needs the L=64 contract run; smaller lattices are
    smoke runs and never eligible (PR #51 Codex 2026-10-17)."""
    if L != CALIBRATION_L:
        return "NOT_A_CONTRACT_RUN_SMOKE_ONLY"
    if verdict == "CALIBRATION_PASS_FOR_COST_MODEL":
        # Spec §6: scaling also needs the Wolff comparison and information per
        # runtime; no preregistered threshold exists for them, so a passing
        # calibration leaves the scale decision UNASSESSED (PR #51 Codex 2026-10-17).
        return "SCALE_UNASSESSED_PENDING_WOLFF_AND_COST_REVIEW"
    return "DO_NOT_SCALE_UNTIL_CALIBRATION_FAILURE_IS_RESOLVED"


def _run_one_wl(L: int, walker: int, e_lo: float, e_hi: float,
                temps: tuple[float, ...], master_seed: int) -> dict:
    _, nbr_list, ei, ej, ax, ay, n = _p41.honeycomb_arrays(L)
    nbins = int(round((e_hi - e_lo) * n))
    prod_sweeps = _prod_sweeps_for(nbins)
    t0 = time.perf_counter()
    res = _p41.wl_entropic_lattice(
        nbr_list, ei, ej, ax, ay, n, L, e_lo, e_hi,
        lnf_final=LNF_FINAL,
        prod_sweeps=prod_sweeps,
        seed=master_seed,
        stream=_stream_for_calibration(walker),
        verbose=False,
    )
    wall = time.perf_counter() - t0
    curves = _p41.upsilon_curves(res, np.asarray(temps))
    leaks = [_p41.canonical_edge_leak(res, t) for t in temps]
    uncovered = [_p41.uncovered_canonical_mass(res, t) for t in temps]
    return {
        "walker": walker,
        "stream_base": _stream_for_calibration(walker),
        "wall_s": wall,
        "wl_sweeps": int(res.wl_sweeps),
        "prod_sweeps": int(res.prod_sweeps),
        "one_over_t_engaged": bool(res.one_over_t_engaged),
        "sweeps_at_1t": int(res.sweeps_at_1t),
        "covered_bins": [int(res.mask.sum()), int(len(res.mask))],
        "leak": [float(x) for x in leaks],
        "uncovered_mass": [float(x) for x in uncovered],
        "y2": [float(x) for x in curves["y2"]],
        "energy_per_site": [float(x) / res.n for x in curves["E"]],
    }


def _run_wolff_crosscheck(L: int, temps: tuple[float, ...],
                          master_seed: int,
                          n_seeds: int = 6,
                          n_measure: int = 300,
                          n_burn: int = 200) -> dict:
    """Canonical comparator.  It is evidence, not an automatic sampler winner."""
    rows = []
    t0 = time.perf_counter()
    for T in temps:
        a = time.perf_counter()
        ref = _p41.wolff_reference(
            L, T, n_seeds=n_seeds, n_measure=n_measure, n_burn=n_burn,
            master_seed=master_seed,
        )
        rows.append({
            "T": T,
            "wall_s": time.perf_counter() - a,
            "energy_per_site": float(ref["E_ps"]),
            "y2": float(ref["y2"]),
            "y2_sem": float(ref["y2_sem"]),
        })
    return {
        "n_seeds": n_seeds,
        "n_measure": n_measure,
        "n_burn": n_burn,
        "wall_s_total": time.perf_counter() - t0,
        "rows": rows,
    }


def run_phy044_calibration(
    L: int = 64,
    n_wl_walkers: int = 2,
    temps=DEFAULT_TEMPS,
    master_seed: int = 42,
    wolff_n_seeds: int = 6,
    wolff_n_measure: int = 300,
    wolff_n_burn: int = 200,
) -> dict:
    """Run the W4 preflight calibration.

    This is deliberately sequential: per-walker wall times remain directly
    interpretable for the cost model.  A later production scheduler may run
    independent walkers in parallel after this calibration is accepted.
    """
    temps = _validate_inputs(L, n_wl_walkers, temps)
    if isinstance(wolff_n_seeds, bool) or not isinstance(wolff_n_seeds, int) or wolff_n_seeds < 2:
        # The Wolff SEM uses ddof=1; one seed gives NaN, which would let the
        # verdict pass without the mandatory Wolff SEM. Refuse before the
        # expensive run (PR #51 Codex round 2026-10-16).
        raise ValueError("Wolff-Quercheck braucht mindestens 2 Seeds (SEM mit ddof=1)")
    # An empty measurement budget averages empty arrays to NaN; the Wolff values
    # never enter the verdict, so the mandatory cross-check would be serialized as
    # null under a passing verdict (PR #51 Codex round 2026-10-17).
    if isinstance(wolff_n_measure, bool) or not isinstance(wolff_n_measure, int) or wolff_n_measure < 1:
        raise ValueError("wolff_n_measure muss eine ganze Zahl >= 1 sein")
    if isinstance(wolff_n_burn, bool) or not isinstance(wolff_n_burn, int) or wolff_n_burn < 0:
        raise ValueError("wolff_n_burn muss eine ganze Zahl >= 0 sein")
    total_t0 = time.perf_counter()

    # Shared anchor/window: same binning for every walker.
    lo_mean, lo_std = _p41.wolff_anchor(L, 0.50, master_seed=master_seed)
    hi_mean, hi_std = _p41.wolff_anchor(L, 0.70, master_seed=master_seed)
    e_lo, e_hi = _p41.window_from_anchors(lo_mean, lo_std, hi_mean, hi_std)
    # Bind bookkeeping to the actual lattice builder, not to an assumed
    # closed-form site count.
    _, _, _, _, _, _, n_sites = _p41.honeycomb_arrays(L)
    nbins = int(round((e_hi - e_lo) * n_sites))

    walkers = [
        _run_one_wl(L, w, e_lo, e_hi, temps, master_seed)
        for w in range(n_wl_walkers)
    ]

    y2 = np.asarray([w["y2"] for w in walkers], dtype=float)
    y2_spread = y2.max(axis=0) - y2.min(axis=0)
    max_spread = float(y2_spread.max())
    # np.max propagates NaN; the builtin max() drops it depending on the order.
    leak_max = float(np.max(np.concatenate([w["leak"] for w in walkers])))
    uncovered_max = float(np.max(np.concatenate([w["uncovered_mass"] for w in walkers])))
    all_1t = all(w["one_over_t_engaged"] and w["sweeps_at_1t"] > 0
                 for w in walkers)

    wolff = _run_wolff_crosscheck(
        L, temps, master_seed,
        n_seeds=wolff_n_seeds,
        n_measure=wolff_n_measure,
        n_burn=wolff_n_burn,
    )

    # Apples-to-apples discrepancy table at fixed T.  We report, but do not
    # collapse, WL walker spread and Wolff SEM into one pseudo-uncertainty.
    comparison = []
    wl_mean = y2.mean(axis=0)
    for k, T in enumerate(temps):
        wr = wolff["rows"][k]
        comparison.append({
            "T": T,
            "wl_y2_mean": float(wl_mean[k]),
            "wl_y2_spread": float(y2_spread[k]),
            "wolff_y2": wr["y2"],
            "wolff_y2_sem": wr["y2_sem"],
            "abs_delta_y2": abs(float(wl_mean[k]) - wr["y2"]),
        })

    verdict = _calibration_verdict(
        all_1t=all_1t,
        leak_max=leak_max,
        uncovered_max=uncovered_max,
        max_y2_spread=max_spread,
    )
    return {
        "module": "PHY044_honeycomb_w4_calibration_v01",
        "date": "2026-09-27",
        "claim_ceiling": "CALIBRATION_ONLY_NOT_T_BKT",
        "spec": "spec/260927 PHI HEX w4 preregistration v01.md",
        "L": L,
        "n_sites": n_sites,
        "n_wl_walkers": n_wl_walkers,
        "production_walker_requirement_met": n_wl_walkers >= 3,
        "calibration_is_not_production_validity": True,
        "master_seed": master_seed,
        "provenance": _provenance(L, n_wl_walkers, wolff_n_seeds),
        "temperatures": list(temps),
        "window": {
            "anchor_lo": {"T": 0.50, "mean_ps": lo_mean, "std_ps": lo_std},
            "anchor_hi": {"T": 0.70, "mean_ps": hi_mean, "std_ps": hi_std},
            "window_ps": [e_lo, e_hi],
            "nbins": nbins,
        },
        "wl": {
            "walkers": walkers,
            "y2_spread_by_T": [float(x) for x in y2_spread],
            "max_y2_spread": max_spread,
            "leak_max": float(leak_max),
            "uncovered_mass_max": float(uncovered_max),
            "all_one_over_t_engaged": bool(all_1t),
            "wall_s_sum": float(sum(w["wall_s"] for w in walkers)),
        },
        "wolff_crosscheck": wolff,
        "comparison": comparison,
        "verdict": verdict,
        "scale_decision": _scale_decision(verdict, L),
        "runtime_s_total": time.perf_counter() - total_t0,
    }


def _clean(obj):
    if isinstance(obj, dict):
        return {str(k): _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        obj = float(obj)
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj


if __name__ == "__main__":
    report = run_phy044_calibration()
    out = _ROOT / "results" / "260927 PHY044 honeycomb w4 calibration report.json"
    out.write_text(
        json.dumps(_clean(report), indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(_clean(report), indent=2, allow_nan=False))
    print(f"\nReport geschrieben: {out}")
