"""PHY050 - deterministic W4-v03a BKT correlation-ratio FSS and G4 recovery.

This module hardens the pre-data W4-v03 analysis contract. It deliberately avoids
an unconstrained spline fit: candidate (T_BKT, c) values are scored by symmetric
cross-size interpolation in the common scaling variable

    X = L * exp(-c / sqrt(T - T_BKT)),  T > T_BKT.

The score is a weighted reduced mismatch, not a formal chi-square p-value because
interpolated comparisons are correlated. A profile-width gate rejects data sets
whose T_BKT is not identifiable even if their raw collapse score is small.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

W4V3_LADDER = (48, 72, 96, 144, 192)
W4V3_T_GRID = tuple(round(0.540 + 0.0025 * k, 4) for k in range(29))
FSS_T_WINDOW = (0.585, 0.610)
FSS_T0_BOUNDS = (0.550, 0.582)
FSS_C_BOUNDS = (0.15, 2.50)
FSS_MIN_PAIR_POINTS = 4
FSS_MIN_COMPARISONS = 60
FSS_PROFILE_DELTA = 0.20
FSS_MAX_PROFILE_WIDTH = 0.012
FSS_MAX_SCORE = 2.50
G4_MAX_ABS_ERROR = 0.003
BOOTSTRAP_SEED = 50_049
MODEL_SIGMA_FLOOR = 0.002
W4V3_SEED_BASE = 49_000_000
W4V3_N_SEEDS = 12
W4V3_N_THERM = 1000
W4V3_N_MEAS = 4000
W4V3_N_BOOT = 1000
DECISION_LOW_EDGE = 0.570
DECISION_HIGH_EDGE = 0.573
MAX_SIGMA_TOT = 0.010
MAX_ROBUST_DELTA = 0.008


@dataclass(frozen=True)
class CollapseFit:
    tbkt: float
    c: float
    score: float
    n_comparisons: int
    n_size_pairs: int
    profile_width: float
    identifiable: bool
    boundary_hit: bool

    @property
    def quotable(self) -> bool:
        return (
            self.identifiable
            and not self.boundary_hit
            and math.isfinite(self.score)
            and self.score <= FSS_MAX_SCORE
        )


def _curve(agg: dict, L: int) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    try:
        t = np.asarray(agg["t_grid"], dtype=float)
        block = agg["curves"][str(L)]
        r = np.asarray(block["R"], dtype=float)
        e = np.asarray(block["SE"], dtype=float)
    except (KeyError, TypeError, ValueError):
        return None
    if t.shape != r.shape or r.shape != e.shape or t.ndim != 1:
        return None
    return t, r, e


def validate_aggregate(agg: dict, require_complete: bool = True) -> bool:
    try:
        ladder = tuple(int(x) for x in agg["ladder"])
        t_grid = tuple(float(x) for x in agg["t_grid"])
    except (KeyError, TypeError, ValueError):
        return False
    if require_complete and ladder != W4V3_LADDER:
        return False
    if require_complete and not np.allclose(
        t_grid, W4V3_T_GRID, atol=0.0, rtol=0.0
    ):
        return False
    if len(ladder) < 3 or len(t_grid) < 6:
        return False
    for L in ladder:
        got = _curve(agg, L)
        if got is None:
            return False
        t, r, e = got
        if len(t) != len(t_grid):
            return False
        if np.any(~np.isfinite(r)) or np.any(~np.isfinite(e)) or np.any(e <= 0.0):
            return False
    return True


def _interp_with_se(
    x: float, xp: np.ndarray, yp: np.ndarray, ep: np.ndarray
) -> tuple[float, float] | None:
    if x < xp[0] or x > xp[-1]:
        return None
    j = int(np.searchsorted(xp, x))
    if j == 0:
        return float(yp[0]), float(ep[0])
    if j == len(xp):
        return float(yp[-1]), float(ep[-1])
    x0, x1 = float(xp[j - 1]), float(xp[j])
    if not x1 > x0:
        return None
    w = (x - x0) / (x1 - x0)
    y = (1.0 - w) * yp[j - 1] + w * yp[j]
    var = (1.0 - w) ** 2 * ep[j - 1] ** 2 + w**2 * ep[j] ** 2
    if not math.isfinite(float(var)) or var <= 0.0:
        return None
    return float(y), math.sqrt(float(var))


def collapse_score(
    agg: dict,
    tbkt: float,
    c: float,
    *,
    ladder: tuple[int, ...] | None = None,
    t_window: tuple[float, float] = FSS_T_WINDOW,
) -> dict | None:
    """Symmetric pairwise collapse mismatch for one (T_BKT, c)."""
    if ladder is None:
        ladder = tuple(int(x) for x in agg.get("ladder", ()))
    lo_t, hi_t = t_window
    if not (FSS_T0_BOUNDS[0] <= tbkt <= FSS_T0_BOUNDS[1]):
        return None
    if not (FSS_C_BOUNDS[0] <= c <= FSS_C_BOUNDS[1]):
        return None
    if not tbkt < lo_t < hi_t:
        return None

    series: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for L in ladder:
        got = _curve(agg, L)
        if got is None:
            return None
        t, r, e = got
        mask = (t >= lo_t) & (t <= hi_t) & (t > tbkt)
        t, r, e = t[mask], r[mask], e[mask]
        if (
            len(t) < 6
            or np.any(~np.isfinite(r))
            or np.any(~np.isfinite(e))
            or np.any(e <= 0.0)
        ):
            return None
        x = L * np.exp(-c / np.sqrt(t - tbkt))
        order = np.argsort(x)
        x, r, e = x[order], r[order], e[order]
        if np.any(np.diff(x) <= 0.0):
            return None
        series[L] = (x, r, e)

    total = 0.0
    n_total = 0
    n_pairs = 0
    for ia, la in enumerate(ladder):
        xa, ya, ea = series[la]
        for lb in ladder[ia + 1 :]:
            xb, yb, eb = series[lb]
            pair_n = 0
            pair_total = 0.0
            for x1, y1, e1, x2, y2, e2 in (
                (xa, ya, ea, xb, yb, eb),
                (xb, yb, eb, xa, ya, ea),
            ):
                for x, y, err in zip(x1, y1, e1, strict=True):
                    pred = _interp_with_se(float(x), x2, y2, e2)
                    if pred is None:
                        continue
                    yhat, ehat = pred
                    var = float(err * err + ehat * ehat)
                    if var <= 0.0 or not math.isfinite(var):
                        return None
                    pair_total += (float(y) - yhat) ** 2 / var
                    pair_n += 1
            if pair_n >= FSS_MIN_PAIR_POINTS:
                total += pair_total
                n_total += pair_n
                n_pairs += 1

    if n_total < FSS_MIN_COMPARISONS or n_pairs < 6:
        return None
    return {
        "score": total / n_total,
        "n_comparisons": n_total,
        "n_size_pairs": n_pairs,
    }


def fit_collapse(
    agg: dict,
    *,
    ladder: tuple[int, ...] | None = None,
    t_window: tuple[float, float] = FSS_T_WINDOW,
    require_complete: bool = True,
) -> CollapseFit | None:
    """Deterministic two-stage grid search with an identifiability profile."""
    if not validate_aggregate(agg, require_complete=require_complete):
        return None
    if ladder is None:
        ladder = tuple(int(x) for x in agg["ladder"])
    if len(ladder) < 3 or any(L not in agg["ladder"] for L in ladder):
        return None

    t0_grid = np.arange(FSS_T0_BOUNDS[0], FSS_T0_BOUNDS[1] + 1e-12, 0.001)
    c_grid = np.geomspace(FSS_C_BOUNDS[0], FSS_C_BOUNDS[1], 80)
    best: tuple[float, float, float, int, int] | None = None
    profile: list[tuple[float, float]] = []

    for t0 in t0_grid:
        t_best = math.inf
        for c in c_grid:
            got = collapse_score(
                agg, float(t0), float(c), ladder=ladder, t_window=t_window
            )
            if got is None:
                continue
            value = float(got["score"])
            t_best = min(t_best, value)
            if best is None or value < best[0]:
                best = (
                    value,
                    float(t0),
                    float(c),
                    int(got["n_comparisons"]),
                    int(got["n_size_pairs"]),
                )
        if math.isfinite(t_best):
            profile.append((float(t0), t_best))

    if best is None:
        return None

    # Profile and profile minimum must come from the same objective grid.
    coarse_best_score = best[0]
    _, bt, bc, _, _ = best
    t_ref = np.arange(
        max(FSS_T0_BOUNDS[0], bt - 0.0025),
        min(FSS_T0_BOUNDS[1], bt + 0.0025) + 1e-12,
        0.0002,
    )
    c_ref = np.linspace(
        max(FSS_C_BOUNDS[0], bc * 0.70),
        min(FSS_C_BOUNDS[1], bc * 1.30),
        31,
    )
    for t0 in t_ref:
        for c in c_ref:
            got = collapse_score(
                agg, float(t0), float(c), ladder=ladder, t_window=t_window
            )
            if got is None:
                continue
            value = float(got["score"])
            if value < best[0]:
                best = (
                    value,
                    float(t0),
                    float(c),
                    int(got["n_comparisons"]),
                    int(got["n_size_pairs"]),
                )

    score_best, bt, bc, nc, npairs = best
    prof_t = [
        t
        for t, value in profile
        if value <= coarse_best_score + FSS_PROFILE_DELTA
    ]
    profile_width = math.inf if not prof_t else max(prof_t) - min(prof_t)
    identifiable = (
        math.isfinite(profile_width) and profile_width <= FSS_MAX_PROFILE_WIDTH
    )
    boundary_hit = (
        bt <= FSS_T0_BOUNDS[0] + 0.0002
        or bt >= FSS_T0_BOUNDS[1] - 0.0002
        or bc <= FSS_C_BOUNDS[0] * 1.02
        or bc >= FSS_C_BOUNDS[1] * 0.98
    )
    return CollapseFit(
        tbkt=bt,
        c=bc,
        score=score_best,
        n_comparisons=nc,
        n_size_pairs=npairs,
        profile_width=profile_width,
        identifiable=identifiable,
        boundary_hit=boundary_hit,
    )


def synthetic_aggregate(
    tbkt: float, c: float, *, noise: float = 0.001, seed: int = 1
) -> dict:
    """Generate BKT-like curves from a known collapse function for G4 only."""
    rng = np.random.default_rng(seed)
    t = np.asarray(W4V3_T_GRID, dtype=float)
    curves = {}
    for L in W4V3_LADDER:
        vals = []
        for temp in t:
            if temp > tbkt:
                x = L * math.exp(-c / math.sqrt(temp - tbkt))
                mean = 0.90 - 0.18 * x / (1.0 + x)
            else:
                mean = 0.90
            vals.append(mean + float(rng.normal(0.0, noise)))
        curves[str(L)] = {"R": vals, "SE": [noise] * len(t)}
    return {
        "ladder": list(W4V3_LADDER),
        "t_grid": list(W4V3_T_GRID),
        "curves": curves,
    }


def flat_null_aggregate(*, noise: float = 0.001, seed: int = 9) -> dict:
    """Size-independent null: low score is possible but T_BKT is unidentifiable."""
    rng = np.random.default_rng(seed)
    t = np.asarray(W4V3_T_GRID, dtype=float)
    curves = {}
    for L in W4V3_LADDER:
        vals = 0.80 + rng.normal(0.0, noise, len(t))
        curves[str(L)] = {"R": vals.tolist(), "SE": [noise] * len(t)}
    return {
        "ladder": list(W4V3_LADDER),
        "t_grid": list(W4V3_T_GRID),
        "curves": curves,
    }


def noncollapse_null_aggregate(*, noise: float = 0.001, seed: int = 17) -> dict:
    """Size-dependent curves without a shared BKT scaling function."""
    rng = np.random.default_rng(seed)
    t = np.asarray(W4V3_T_GRID, dtype=float)
    curves = {}
    for idx, L in enumerate(W4V3_LADDER):
        z = t - 0.575
        vals = (
            0.82
            - (0.5 + 0.2 * idx) * z
            - (2.0 - 0.25 * idx) * z * z
            + rng.normal(0.0, noise, len(t))
        )
        curves[str(L)] = {"R": vals.tolist(), "SE": [noise] * len(t)}
    return {
        "ladder": list(W4V3_LADDER),
        "t_grid": list(W4V3_T_GRID),
        "curves": curves,
    }


def g4_synthetic_recovery() -> dict:
    cases = (
        (0.565, 0.90, 1_049),
        (0.573, 0.90, 2_049),
        (0.579, 0.75, 3_049),
    )
    recovered = []
    pass_recovery = True
    for truth, c, seed in cases:
        fit = fit_collapse(synthetic_aggregate(truth, c, seed=seed))
        err = math.inf if fit is None else abs(fit.tbkt - truth)
        ok = fit is not None and fit.quotable and err <= G4_MAX_ABS_ERROR
        pass_recovery = pass_recovery and ok
        recovered.append(
            {
                "truth": truth,
                "c_truth": c,
                "fit": None if fit is None else fit.tbkt,
                "score": None if fit is None else fit.score,
                "abs_error": err,
                "quotable": False if fit is None else fit.quotable,
                "pass": ok,
            }
        )

    missing = synthetic_aggregate(0.573, 0.90, seed=4_049)
    missing["ladder"] = missing["ladder"][:-1]
    missing["curves"].pop("192")
    missing_fail_closed = fit_collapse(missing) is None

    flat_fit = fit_collapse(flat_null_aggregate())
    flat_rejected = (
        flat_fit is not None
        and not flat_fit.identifiable
        and not flat_fit.quotable
    )
    noncollapse_fit = fit_collapse(noncollapse_null_aggregate())
    noncollapse_rejected = (
        noncollapse_fit is not None
        and not noncollapse_fit.identifiable
        and not noncollapse_fit.quotable
    )

    passed = (
        pass_recovery
        and missing_fail_closed
        and flat_rejected
        and noncollapse_rejected
    )
    return {
        "gate": "G4_SYNTHETIC_RECOVERY",
        "passed": passed,
        "recovery_cases": recovered,
        "missing_input_fail_closed": missing_fail_closed,
        "flat_null_rejected_as_unidentifiable": flat_rejected,
        "noncollapse_null_rejected_as_unidentifiable": noncollapse_rejected,
        "claim_ceiling": (
            "G4 validates estimator mechanics only; no PHI-Hex physics claim without "
            "production G0/G5/G6 and a committed gate log."
        ),
    }


def _product_groups(prod: dict) -> dict[tuple[int, int], list[dict]] | None:
    """Validate the exact preregistered raw-production contract, fail closed."""
    try:
        ladder = tuple(int(x) for x in prod["ladder"])
        t_grid = tuple(float(x) for x in prod["t_grid"])
        n_seeds = int(prod["n_seeds"])
        n_therm = int(prod["n_therm"])
        n_meas = int(prod["n_meas"])
        rows = list(prod["rows"])
    except (KeyError, TypeError, ValueError):
        return None

    if ladder != W4V3_LADDER or not np.allclose(
        t_grid, W4V3_T_GRID, atol=0.0, rtol=0.0
    ):
        return None
    if (
        n_seeds != W4V3_N_SEEDS
        or n_therm != W4V3_N_THERM
        or n_meas != W4V3_N_MEAS
    ):
        return None
    if prod.get("complete") is not True or prod.get("unmeasured") != []:
        return None

    expected_rows = len(W4V3_LADDER) * len(W4V3_T_GRID) * W4V3_N_SEEDS
    if len(rows) != expected_rows:
        return None

    groups: dict[tuple[int, int], list[dict]] = {}
    try:
        for L in ladder:
            for k in range(len(t_grid)):
                rr = [
                    r
                    for r in rows
                    if int(r.get("L", -1)) == L
                    and int(r.get("t_idx", -1)) == k
                ]
                if len(rr) != W4V3_N_SEEDS:
                    return None
                rr.sort(key=lambda r: int(r.get("s", -1)))
                if [int(r.get("s", -1)) for r in rr] != list(
                    range(W4V3_N_SEEDS)
                ):
                    return None
                expected_t = t_grid[k]
                for r in rr:
                    seed_idx = int(r.get("s", -1))
                    expected_seed = (
                        W4V3_SEED_BASE + 1000 * L + 100 * k + seed_idx
                    )
                    if float(r.get("T", math.nan)) != expected_t:
                        return None
                    if int(r.get("seed", -1)) != expected_seed:
                        return None
                    q = float(r["g_quarter"])
                    h = float(r["g_half"])
                    if not math.isfinite(q) or not math.isfinite(h):
                        return None
                groups[(L, k)] = rr
    except (KeyError, TypeError, ValueError, OverflowError):
        return None
    return groups

def product_to_aggregate(
    prod: dict, *, rng: np.random.Generator | None = None
) -> dict | None:
    """Build R curves from seed rows; optional paired resampling per (L,T)."""
    groups = _product_groups(prod)
    if groups is None:
        return None
    ns = int(prod["n_seeds"])
    curves = {}
    for L in W4V3_LADDER:
        vals, sems = [], []
        for k in range(len(W4V3_T_GRID)):
            rr = groups[(L, k)]
            try:
                q0 = np.asarray(
                    [float(r["g_quarter"]) for r in rr], dtype=float
                )
                h0 = np.asarray(
                    [float(r["g_half"]) for r in rr], dtype=float
                )
            except (KeyError, TypeError, ValueError, OverflowError):
                return None
            if np.any(~np.isfinite(q0)) or np.any(~np.isfinite(h0)):
                return None
            idx = np.arange(ns) if rng is None else rng.integers(0, ns, size=ns)
            q, h = q0[idx], h0[idx]
            den = float(q.mean())
            if abs(den) <= 1e-6 or not math.isfinite(den):
                return None
            vals.append(float(h.mean()) / den)

            # Fixed original-sample jackknife SE avoids nested bootstrap variance.
            jk = []
            for s in range(ns):
                qq = np.delete(q0, s)
                hh = np.delete(h0, s)
                d = float(qq.mean())
                if abs(d) <= 1e-6:
                    return None
                jk.append(float(hh.mean()) / d)
            jk = np.asarray(jk)
            jm = float(jk.mean())
            se = math.sqrt(
                (ns - 1) / ns * float(np.sum((jk - jm) ** 2))
            )
            if not math.isfinite(se) or se <= 0.0:
                return None
            sems.append(se)
        curves[str(L)] = {"R": vals, "SE": sems}
    return {
        "ladder": list(W4V3_LADDER),
        "t_grid": list(W4V3_T_GRID),
        "curves": curves,
    }


def bootstrap_tbkt(
    prod: dict, *, n_boot: int = W4V3_N_BOOT, seed: int = BOOTSTRAP_SEED
) -> dict | None:
    """Paired-seed bootstrap; invalid/non-identifiable replicas fail closed."""
    if n_boot < 20:
        raise ValueError("n_boot must be >= 20")
    base = product_to_aggregate(prod)
    if base is None:
        return None
    primary = fit_collapse(base)
    if primary is None or not primary.quotable:
        return None
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        agg = product_to_aggregate(prod, rng=rng)
        if agg is None:
            continue
        fit = fit_collapse(agg)
        if fit is not None and fit.quotable:
            vals.append(fit.tbkt)
    min_valid = math.ceil(0.90 * n_boot)
    if len(vals) < min_valid:
        return None
    a = np.asarray(vals, dtype=float)
    return {
        "tbkt": primary.tbkt,
        "c": primary.c,
        "collapse_score": primary.score,
        "sigma_boot": float(np.std(a, ddof=1)),
        "bootstrap_valid": len(vals),
        "bootstrap_requested": n_boot,
        "q025": float(np.quantile(a, 0.025)),
        "q975": float(np.quantile(a, 0.975)),
    }


def robustness_variants(agg: dict, primary: CollapseFit) -> dict:
    variants = {}
    specs = {
        "window_hi_trim": (tuple(W4V3_LADDER), (0.585, 0.6075), True),
        "window_lo_raise": (tuple(W4V3_LADDER), (0.5875, 0.610), True),
    }
    for ldrop in W4V3_LADDER:
        ladder = tuple(L for L in W4V3_LADDER if L != ldrop)
        specs[f"leave_L{ldrop}_out"] = (ladder, FSS_T_WINDOW, False)

    diffs = []
    for name, (ladder, window, require_complete) in specs.items():
        fit = fit_collapse(
            agg,
            ladder=ladder,
            t_window=window,
            require_complete=require_complete,
        )
        ok = fit is not None and fit.quotable
        delta = None if not ok else abs(fit.tbkt - primary.tbkt)
        if delta is not None:
            diffs.append(delta)
        variants[name] = {
            "tbkt": None if fit is None else fit.tbkt,
            "score": None if fit is None else fit.score,
            "quotable": ok,
            "delta_primary": delta,
        }
    max_delta = (
        math.inf if len(diffs) != len(specs) else max(diffs, default=0.0)
    )
    return {
        "variants": variants,
        "max_delta": max_delta,
        "passed": max_delta <= MAX_ROBUST_DELTA,
    }



def decision_label(
    tbkt: float,
    sigma_tot: float,
    *,
    gates_passed: bool = True,
) -> str:
    """Internal discrimination label with the full hypothesis overlap respected."""
    if (
        not gates_passed
        or not math.isfinite(tbkt)
        or not math.isfinite(sigma_tot)
        or sigma_tot < 0.0
        or sigma_tot > MAX_SIGMA_TOT
    ):
        return "INCONCLUSIVE"
    lo = tbkt - 2.0 * sigma_tot
    hi = tbkt + 2.0 * sigma_tot
    if hi < DECISION_LOW_EDGE:
        return "SUPPORTED_LOW"
    if lo > DECISION_HIGH_EDGE:
        return "SUPPORTED_LITERATURE"
    # Any 95% interval intersecting the preregistered H_A/H_B overlap
    # [0.570, 0.573] stays explicitly non-discriminating.
    if hi >= DECISION_LOW_EDGE and lo <= DECISION_HIGH_EDGE:
        return "OVERLAP"
    return "INCONCLUSIVE"


def assess_production(prod: dict, *, n_boot: int = W4V3_N_BOOT) -> dict:
    """Single fail-closed G0/G4/G5/G6 adjudicator for real PHY049 data."""
    g4 = g4_synthetic_recovery()
    agg = product_to_aggregate(prod)
    g0 = agg is not None
    if not g0:
        return {
            "gates": {"G0_input": False, "G4_synthetic_recovery": bool(g4["passed"]),
                      "G5_power": False, "G6_robustness": False},
            "decision": "INCONCLUSIVE",
            "physics_interpretation_enabled": False,
            "claim_ceiling": "G0 failed: incomplete or invalid production evidence.",
        }

    if n_boot != W4V3_N_BOOT:
        return {
            "gates": {
                "G0_input": True,
                "G4_synthetic_recovery": bool(g4["passed"]),
                "G5_power": False,
                "G6_robustness": False,
            },
            "decision": "INCONCLUSIVE",
            "physics_interpretation_enabled": False,
            "bootstrap_requested": n_boot,
            "bootstrap_required": W4V3_N_BOOT,
            "claim_ceiling": (
                "Non-production bootstrap count: production adjudication "
                f"requires exactly {W4V3_N_BOOT} replicas."
            ),
        }

    primary = fit_collapse(agg)
    if primary is None or not primary.quotable:
        return {
            "gates": {"G0_input": True, "G4_synthetic_recovery": bool(g4["passed"]),
                      "G5_power": False, "G6_robustness": False},
            "decision": "INCONCLUSIVE",
            "physics_interpretation_enabled": False,
            "claim_ceiling": "Primary collapse is non-quotable.",
        }

    boot = bootstrap_tbkt(prod, n_boot=n_boot)
    robust = robustness_variants(agg, primary)
    if boot is None:
        sigma_boot = math.inf
    else:
        sigma_boot = float(boot["sigma_boot"])

    variant_values = [primary.tbkt]
    for item in robust["variants"].values():
        if item["quotable"] and item["tbkt"] is not None:
            variant_values.append(float(item["tbkt"]))
    all_variants_quotable = len(variant_values) == 1 + len(robust["variants"])
    sigma_model = (
        max(MODEL_SIGMA_FLOOR,
            0.5 * (max(variant_values) - min(variant_values)))
        if all_variants_quotable
        else math.inf
    )
    sigma_tot = math.hypot(sigma_boot, sigma_model)
    g5 = bool(math.isfinite(sigma_tot) and sigma_tot <= MAX_SIGMA_TOT)
    g6 = bool(robust["passed"] and all_variants_quotable)
    gates = {
        "G0_input": True,
        "G4_synthetic_recovery": bool(g4["passed"]),
        "G5_power": g5,
        "G6_robustness": g6,
    }
    enabled = all(gates.values())
    return {
        "gates": gates,
        "primary": {
            "tbkt": primary.tbkt,
            "c": primary.c,
            "score": primary.score,
            "profile_width": primary.profile_width,
        },
        "bootstrap": boot,
        "robustness": robust,
        "sigma_boot": sigma_boot,
        "sigma_model": sigma_model,
        "sigma_tot": sigma_tot,
        "decision": decision_label(primary.tbkt, sigma_tot, gates_passed=enabled),
        "physics_interpretation_enabled": enabled,
        "claim_ceiling": (
            "INTERNAL_W4_V03_DISCRIMINATION_ONLY; not an external best-value claim."
            if enabled
            else "NO_PHYSICS_INTERPRETATION until G0/G4/G5/G6 are green."
        ),
    }


def preflight() -> dict:
    g4 = g4_synthetic_recovery()
    return {
        "module": "PHY050_correlation_ratio_fss_recovery_v01",
        "gates": {"G4_synthetic_recovery": bool(g4["passed"])},
        "g4": g4,
        "production_interpretation_enabled": False,
        "claim_ceiling": (
            "NO_PHYSICS_INTERPRETATION until PHY049 production data satisfy "
            "G0 input, G5 power and G6 robustness with committed gate evidence."
        ),
    }


if __name__ == "__main__":
    import json

    print(json.dumps(preflight(), indent=2))
