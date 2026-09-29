"""PHY051 - triangular + kagome T_BKT mit Upsilon PRO FLAECHE (Wolff-Numba,
Paar-Crossings, gewichtetes Paar-Mittel) als O1-Diskriminator gegen
normierungsfreie Literaturwerte - exakt nach Vorregistrierung.

Coworker Research / Coworkerz, 28. September 2026
de-CH Konventionen, ASCII-Quotes, ss statt Eszett, keine Personennamen.

Vertrag: `spec/260928 PHI HEX phy051 triangular kagome per-area
preregistration v01.md` (VOR der Datennahme committet). Normierung:
`spec/260926 PHI HEX O1 helicity normalization decision v01.md`. Kernel und
Gitter-Builder: PHY045 (validiert: T->0-Orakel, Quercheck gegen den
Python-Wolff von PHY031). Schaetzer-Bausteine 1:1 aus PHY046 (HKS-Fit,
WM-Fit mit freiem C) und PHY040 (C-eliminiertes Paar-Crossing).

Ablauf je Gitter: je (L, T, Seed) ein Wolff-Lauf -> per-Seed-Upsilon_A und
Upsilon_site; Seed-Mittelkurven; Paar-Crossings (L1, L2) je Kanal;
Primaer = gewichtetes Mittel der grossen Paare (Lehre PHY048); Varianten
v1..v5; Jackknife ueber Seeds; Verdikt je Kanal gegen das vorregistrierte
Band; O1-Wahrheitstafel (Spec Abschnitt 6). Kein Parameter dieses Moduls darf nach
Sicht der PHY051-Daten geaendert werden.

EHRLICHKEIT: FINDING relativ zu vorregistrierten Baendern, kein Bestwert;
Pipeline-Gates pruefen Integritaet, nicht Physik-Wahrheit.
"""
from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import datetime  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import importlib.util  # noqa: E402
import multiprocessing  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

_SRC = Path(__file__).resolve().parent
_ROOT = _SRC.parent


def _load(name: str, filename: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _SRC / filename)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


_p45 = _load("phy045_normalization",
             "260926 PHY045 helicity normalization O1 test v01.py")
_p46 = _load("phy046_w4_v02",
             "260926 PHY046 honeycomb w4 wolff area-helicity v01.py")

tbkt_pair_from_curves = _p46.tbkt_pair_from_curves
fit_hks = _p46.fit_hks
wm_free_c_fit = _p46.wm_free_c_fit
AREA_PER_SITE = _p45.AREA_PER_SITE
UPSILON0_SITE = _p45.UPSILON0_SITE

# ---- vorregistrierte Konstanten (Spec v01; Test-gebunden) -----------------
PHY051_NORMALIZATION = "per_area"
PHY051_N_SEEDS = 8
PHY051_N_THERM = 500
PHY051_N_MEAS = 1500
PHY051_SEED_BASE = 51_000_000
PHY051_LAT_ID = {"triangular": 0, "kagome": 1}
PHY051_SIGMA_FSS_FLOOR_REL = 0.005     # 0.5 % von T_P
PHY051_SIGMA_MAX_REL = 0.010           # S2: sigma_tot > 1 % von T_P
PHY051_MIN_PAIRS = 3                   # S1
PHY051_MIN_PRIMARY_PAIRS = 2           # S1b
PHY051_HKS_B_RANGE = (0.1, 10.0)
PHY051_WALL_BUDGET_H = 3.0             # S0, je Gitter, 4 Prozesse


def _grid(t0: float, t1: float, step: float = 0.005) -> tuple:
    n = int(round((t1 - t0) / step)) + 1
    return tuple(round(t0 + step * k, 3) for k in range(n))


PHY051_PLAN = {
    "triangular": {
        "ladder": (33, 49, 65, 97, 129, 193, 257),
        "pairs": ((33, 65), (49, 97), (65, 129), (97, 193), (129, 257)),
        "t_grid": _grid(1.380, 1.510),
        "band": (1.450, 1.480),
        "L_min_primary": 65,
        "ref_free": {"label": "HT-Reihe Butera & Pernici arXiv:0806.1496, "
                              "J_c = 0.6824(8)",
                     "T": 1.4654, "sigma": 0.0017},
        "ref_free_alt": {"label": "MC arXiv:1010.3075, J_c = 0.6833(6)",
                         "T": 1.4635, "sigma": 0.0013},
        "ref_helicity": {"label": "Sorokin (zitiert in arXiv:2305.00651), "
                                  "Helicity, Normierung unbelegt",
                         "T": 1.418, "sigma": 0.002},
    },
    "kagome": {
        "ladder": (32, 48, 64, 96, 128, 192),
        "pairs": ((32, 64), (48, 96), (64, 128), (96, 192)),
        "t_grid": _grid(0.780, 0.860),
        "band": (0.808, 0.842),
        "L_min_primary": 64,
        "ref_free": {"label": "Okabe & Otsuka arXiv:2501.07388 Tab. 1, "
                              "rough estimate (kein Fehlerbalken)",
                     "T": 0.825, "sigma": None},
        "ref_free_alt": None,
        "ref_helicity": None,
    },
}


def seed_for(lattice: str, L: int, t_idx: int, s: int) -> int:
    """Seed-Vertrag v01: 51_000_000 + 1_000_000*lat_id + 1000*L + 10*t_idx + s."""
    return (PHY051_SEED_BASE + 1_000_000 * PHY051_LAT_ID[lattice]
            + 1000 * L + 10 * t_idx + s)


# ============================================================================
# 1. Produktion (L-Bloecke absteigend, Budget-Stop S0 zwischen Bloecken)
# ============================================================================

def _job(args: tuple) -> dict:
    lattice, L, t_idx, T, s, n_therm, n_meas = args
    lat = _p45.build(lattice, L)
    nbr, deg = _p45._nbr_arrays(lat)
    th = np.zeros(lat.n)
    run = _p45._nb_run if _p45.HAVE_NUMBA else _p45._py_run
    t0 = time.perf_counter()
    d = run(th, nbr, deg, lat.ei, lat.ej,
            np.ascontiguousarray(lat.disp[:, 0]),
            np.ascontiguousarray(lat.disp[:, 1]),
            1.0 / T, n_therm, n_meas, seed_for(lattice, L, t_idx, s))
    beta = 1.0 / T
    num = 0.5 * ((d[:, 1].mean() - beta * np.mean(d[:, 2] ** 2))
                 + (d[:, 3].mean() - beta * np.mean(d[:, 4] ** 2)))
    return {"L": L, "t_idx": t_idx, "s": s,
            "ups_area": float(num / lat.torus_area),
            "ups_site": float(num / lat.n),
            "wall_s": time.perf_counter() - t0}


def _pool(max_workers: int):
    """Prozess-Pool NUR mit fork-Kontext: dieses Modul ist unter einem
    Kunstnamen per importlib geladen (Dateiname mit Leerzeichen) und damit
    fuer spawn/forkserver-Kinder (Windows; Python-3.14-Default auf Linux)
    NICHT per Modulname importierbar. Ohne fork wird sequenziell gerechnet -
    die Ergebnisse sind per (Gitter, L, T, Seed) deterministisch und vom
    Prozessmodell unabhaengig; nur die Wall-Zeit aendert sich."""
    if max_workers <= 1:
        return None
    try:
        ctx = multiprocessing.get_context("fork")
    except ValueError:
        print("PHY051: kein fork-Kontext verfuegbar -> sequenziell (Wall "
              "laenger, Ergebnisse identisch)", flush=True)
        return None
    return ProcessPoolExecutor(max_workers=max_workers, mp_context=ctx)


def produce(lattice: str, ladder=None, t_grid=None, n_seeds=PHY051_N_SEEDS,
            n_therm=PHY051_N_THERM, n_meas=PHY051_N_MEAS, max_workers=4,
            wall_budget_h=PHY051_WALL_BUDGET_H) -> dict:
    plan = PHY051_PLAN[lattice]
    ladder = tuple(plan["ladder"] if ladder is None else ladder)
    t_grid = tuple(plan["t_grid"] if t_grid is None else t_grid)
    area = {L: np.full((len(t_grid), n_seeds), np.nan) for L in ladder}
    site = {L: np.full((len(t_grid), n_seeds), np.nan) for L in ladder}
    unmeasured, cpu = [], 0.0
    run_utc_start = datetime.datetime.now(datetime.timezone.utc).isoformat(
        timespec="seconds")
    t_start = time.perf_counter()
    ex = _pool(max_workers)
    try:
        done = 0
        for L in sorted(ladder, reverse=True):
            # S0: Budget wird ZWISCHEN L-Bloecken geprueft; der erste Block
            # laeuft immer (sonst gaebe es keine Evidenz fuer den Report).
            if done and (time.perf_counter() - t_start) / 3600.0 > wall_budget_h:
                unmeasured.append(L)
                continue
            jobs = [(lattice, L, k, T, s, n_therm, n_meas)
                    for k, T in enumerate(t_grid) for s in range(n_seeds)]
            res = (list(ex.map(_job, jobs, chunksize=1)) if ex
                   else [_job(j) for j in jobs])
            for r in res:
                area[L][r["t_idx"], r["s"]] = r["ups_area"]
                site[L][r["t_idx"], r["s"]] = r["ups_site"]
                cpu += r["wall_s"]
            done += 1
    finally:
        if ex:
            ex.shutdown()
    return {"lattice": lattice, "t_grid": list(t_grid),
            "ladder": list(ladder), "n_seeds": n_seeds,
            "n_therm": n_therm, "n_meas": n_meas,
            "ups_area": {str(L): area[L].tolist() for L in ladder},
            "ups_site": {str(L): site[L].tolist() for L in ladder},
            "unmeasured_L": sorted(unmeasured),
            "wall_s": time.perf_counter() - t_start, "cpu_s": cpu,
            "max_workers": max_workers, "have_numba": bool(_p45.HAVE_NUMBA),
            # Lauf-Zeitstempel (Reality-Anchor); Laeufe vor Einfuehrung des
            # Feldes tragen None (Zeitpunkt dann aus Commit-Historie).
            "run_utc_start": run_utc_start,
            "run_utc_end": datetime.datetime.now(
                datetime.timezone.utc).isoformat(timespec="seconds")}


# ============================================================================
# 2. Auswertung (rein, deterministisch aus den Produktionsdaten)
# ============================================================================

def pair_crossings(curves: dict, t_grid, pairs) -> dict:
    t = np.asarray(t_grid, dtype=float)
    out = {}
    for a, b in pairs:
        if a in curves and b in curves:
            out[(a, b)] = tbkt_pair_from_curves(t, curves[a], a, curves[b], b)
        else:
            out[(a, b)] = None
    return out


def weighted_pair_mean(pc: dict, pair_sig: dict | None, L_min: int) -> dict:
    """Inverse-varianz-gewichtetes Mittel der Paar-Crossings mit L1 >= L_min.
    Fehlt eine SE (oder ist sie <= 0), wird ungewichtet gemittelt."""
    pts = [((a, b), T) for (a, b), T in pc.items()
           if T is not None and a >= L_min]
    if not pts:
        return {"T": None, "n_pairs": 0, "weighted": False}
    sig = ([pair_sig.get(p) for p, _ in pts] if pair_sig is not None
           else [None] * len(pts))
    Ts = np.array([T for _, T in pts], dtype=float)
    if any(x is None or not (x > 0.0) for x in sig):
        return {"T": float(Ts.mean()), "n_pairs": len(pts),
                "weighted": False}
    w = 1.0 / np.asarray(sig, dtype=float) ** 2
    return {"T": float(np.sum(w * Ts) / np.sum(w)), "n_pairs": len(pts),
            "weighted": True}


def estimate(curves: dict, curves_se: dict, t_grid, plan: dict,
             pair_sig: dict | None = None) -> dict:
    """Primaer (gewichtetes Mittel grosser Paare) + Varianten v1..v5."""
    pairs = plan["pairs"]
    ladder = [L for L in plan["ladder"] if L in curves]
    pc = pair_crossings(curves, t_grid, pairs)
    pts = [(a, T) for (a, b), T in pc.items() if T is not None]
    Ls = [a for a, _ in pts]
    Ts = [T for _, T in pts]
    sig = None
    if pair_sig is not None:
        key_of = {a: (a, b) for (a, b) in pc}       # L1 ist je Paar eindeutig
        sig = [pair_sig.get(key_of[a]) for a in Ls]
        # fehlt fuer ein Paar die Jackknife-SE (oder ist sie 0), wird
        # ungewichtet gefittet - kein Mischen von echten und Ersatz-Gewichten
        if any(x is None or x <= 0.0 for x in sig):
            sig = None
    primary = weighted_pair_mean(pc, pair_sig, plan["L_min_primary"])
    v1 = weighted_pair_mean(pc, pair_sig, 0)
    v2 = pc.get(tuple(pairs[-1]))        # groesstes Paar allein (Spec 5, v2)
    v3 = fit_hks(Ls, Ts, sig, b_range=PHY051_HKS_B_RANGE)
    v4 = fit_hks(Ls, Ts, sig, b_fixed=1.0)
    v5 = wm_free_c_fit(curves, curves_se, t_grid, ladder) if ladder else None
    return {"pairs": {f"{a}_{b}": T for (a, b), T in pc.items()},
            "n_pairs_with_crossing": len(pts),
            "primary_weighted_large_pairs": primary,
            "v1_weighted_all_pairs": v1, "v2_largest_pair": v2,
            "v3_hks_3par": v3, "v4_hks_b1": v4, "v5_wm_free_c": v5}


def _jk_se(vals: list) -> float | None:
    v = [x for x in vals if x is not None]
    n = len(v)
    if n < 2:
        return None
    m = float(np.mean(v))
    return float(math.sqrt((n - 1) / n * sum((x - m) ** 2 for x in v)))


def verdict(t_p: float | None, sigma_sampler: float | None,
            sigma_fss: float | None, n_pairs: int, n_primary_pairs: int,
            band: tuple, ladder_complete: bool = True) -> dict:
    """Regeln Spec Abschnitt 6 in fester Reihenfolge: S0, S1, S1b, S2, dann C/I."""
    if not ladder_complete:
        return {"verdict": "NEGATIVE_RESULT", "rule": "S0", "sigma_tot": None}
    if n_pairs < PHY051_MIN_PAIRS:
        return {"verdict": "NEGATIVE_RESULT", "rule": "S1", "sigma_tot": None}
    if n_primary_pairs < PHY051_MIN_PRIMARY_PAIRS or t_p is None:
        return {"verdict": "NEGATIVE_RESULT", "rule": "S1b", "sigma_tot": None}
    if sigma_sampler is None or sigma_fss is None:
        return {"verdict": "NEGATIVE_RESULT", "rule": "S2", "sigma_tot": None}
    s_fss = max(float(sigma_fss), PHY051_SIGMA_FSS_FLOOR_REL * t_p)
    s_tot = math.hypot(float(sigma_sampler), s_fss)
    if s_tot > PHY051_SIGMA_MAX_REL * t_p:
        return {"verdict": "NEGATIVE_RESULT", "rule": "S2", "sigma_tot": s_tot}
    lo, hi = t_p - 2.0 * s_tot, t_p + 2.0 * s_tot
    ok = hi >= band[0] and lo <= band[1]
    return {"verdict": "CONSISTENT" if ok else "INCONSISTENT",
            "rule": "C" if ok else "I", "sigma_tot": s_tot,
            "sigma_fss_floored": s_fss, "interval": [lo, hi],
            "band": list(band)}


def o1_label(area_verdict: str, site_verdict: str) -> str:
    """Wahrheitstafel Spec Abschnitt 6 (fails-closed)."""
    if "NEGATIVE_RESULT" in (area_verdict, site_verdict):
        return "NEGATIVE_RESULT"
    table = {("CONSISTENT", "INCONSISTENT"): "O1_PER_AREA_CORROBORATED",
             ("INCONSISTENT", "CONSISTENT"): "O1_CHALLENGED",
             ("CONSISTENT", "CONSISTENT"): "NON_DISCRIMINATING",
             ("INCONSISTENT", "INCONSISTENT"): "TENSION_BOTH_CHANNELS"}
    return table.get((area_verdict, site_verdict), "NEGATIVE_RESULT")


def _z(t: float | None, ref: dict | None, sigma_tot: float | None):
    if t is None or ref is None or sigma_tot is None:
        return None
    s = math.hypot(sigma_tot, ref["sigma"] or 0.0)
    return float((t - ref["T"]) / s) if s > 0 else None


def analyse(prod: dict, key: str = "ups_area") -> dict:
    """Vollstaendige vorregistrierte Auswertung eines Kanals."""
    lattice = prod["lattice"]
    plan = PHY051_PLAN[lattice]
    t_grid = prod["t_grid"]
    unmeasured = set(int(x) for x in prod.get("unmeasured_L", []))
    ladder = [int(x) for x in prod["ladder"] if int(x) not in unmeasured]
    ns = prod["n_seeds"]
    data = {L: np.asarray(prod[key][str(L)], dtype=float) for L in ladder}
    ladder = [L for L in ladder if np.all(np.isfinite(data[L]))]
    data = {L: data[L] for L in ladder}
    complete = (set(ladder) == set(int(x) for x in plan["ladder"]))
    mean = {L: data[L].mean(axis=1) for L in ladder}
    se = {L: data[L].std(axis=1, ddof=1) / math.sqrt(ns) for L in ladder}
    jk_curves = [{L: np.delete(data[L], s, axis=1).mean(axis=1)
                  for L in ladder} for s in range(ns)]
    pairs = plan["pairs"]
    pc_full = pair_crossings(mean, t_grid, pairs)
    pc_jk = [pair_crossings(c, t_grid, pairs) for c in jk_curves]
    pair_sig = {p: _jk_se([j[p] for j in pc_jk]) for p in pc_full}
    full = estimate(mean, se, t_grid, plan, pair_sig)
    jk_T = [estimate(c, se, t_grid, plan, pair_sig)
            ["primary_weighted_large_pairs"]["T"] for c in jk_curves]
    prim = full["primary_weighted_large_pairs"]
    t_p = prim["T"]
    variants = [t_p, full["v1_weighted_all_pairs"]["T"],
                full["v2_largest_pair"],
                None if full["v3_hks_3par"] is None
                else full["v3_hks_3par"]["T_c"],
                None if full["v4_hks_b1"] is None
                else full["v4_hks_b1"]["T_c"],
                None if full["v5_wm_free_c"] is None
                else full["v5_wm_free_c"]["T"]]
    vv = [v for v in variants if v is not None]
    sigma_fss = 0.5 * (max(vv) - min(vv)) if len(vv) >= 2 else None
    # Diagnostik (Nachtrag v01a, vor Sicht der Daten; NICHT entscheidend):
    # halbe Spannweite OHNE die 3-Parameter-HKS-Variante v3, die PHY048 auf
    # square als rauschverstaerkend ausgewiesen hat. Zeigt, ob ein S2-Stop
    # allein von v3 getragen wuerde. Das Verdikt bleibt bei sigma_fss (v01).
    vv_no3 = [v for k, v in enumerate(variants) if v is not None and k != 3]
    sigma_fss_no_hks3 = (0.5 * (max(vv_no3) - min(vv_no3))
                         if len(vv_no3) >= 2 else None)
    sigma_sampler = _jk_se(jk_T)
    vd = verdict(t_p, sigma_sampler, sigma_fss, full["n_pairs_with_crossing"],
                 prim["n_pairs"], plan["band"], ladder_complete=complete)
    vd_no3 = verdict(t_p, sigma_sampler, sigma_fss_no_hks3,
                     full["n_pairs_with_crossing"], prim["n_pairs"],
                     plan["band"], ladder_complete=complete)
    return {"channel": key, "ladder_used": ladder,
            "ladder_complete": complete, "estimates": full,
            "pair_jackknife_se": {f"{a}_{b}": v
                                  for (a, b), v in pair_sig.items()},
            "T_P": t_p, "sigma_sampler": sigma_sampler,
            "sigma_fss_raw": sigma_fss, "variants": variants,
            "verdict": vd,
            "diagnostic_without_hks3par": {"sigma_fss_raw": sigma_fss_no_hks3,
                                           "verdict": vd_no3,
                                           "decisive": False},
            "z_vs_ref_free": _z(t_p, plan["ref_free"], vd.get("sigma_tot")),
            "z_vs_ref_free_alt": _z(t_p, plan["ref_free_alt"],
                                    vd.get("sigma_tot")),
            "z_vs_ref_helicity": _z(t_p, plan["ref_helicity"],
                                    vd.get("sigma_tot")),
            "curve_mean": {str(L): mean[L].tolist() for L in ladder},
            "curve_se": {str(L): se[L].tolist() for L in ladder},
            "curve_se_max": {str(L): float(np.max(se[L])) for L in ladder}}


# ============================================================================
# 3. Gates (Pipeline-Integritaet, Spec Abschnitt 7) + Report
# ============================================================================

def geometry_oracle(lattice: str, ladder) -> dict:
    out = {}
    for L in ladder:
        lat = _p45.build(lattice, int(L))
        y0 = _p45.upsilon0_site_from_geometry(lat)
        out[str(L)] = {"upsilon0_site": y0,
                       "upsilon0_area": y0 / lat.area_per_site,
                       "area_per_site": lat.area_per_site, "n": lat.n}
    return out


def gates(prod: dict, area: dict, oracle: dict) -> dict:
    lattice = prod["lattice"]
    plan = PHY051_PLAN[lattice]
    t_grid = prod["t_grid"]
    ns = prod["n_seeds"]
    ladder = [int(x) for x in prod["ladder"]]
    shapes_ok = all(np.asarray(prod["ups_area"][str(L)]).shape
                    == (len(t_grid), ns) for L in ladder)
    finite = all(np.all(np.isfinite(np.asarray(prod[k][str(L)], dtype=float)))
                 for L in ladder for k in ("ups_area", "ups_site"))
    y0 = UPSILON0_SITE[lattice]
    oracle_ok = all(abs(v["upsilon0_site"] - y0) < 1e-12
                    and abs(v["area_per_site"] - AREA_PER_SITE[lattice]) < 1e-12
                    for v in oracle.values())
    used = area["ladder_used"]
    cm = area["curve_mean"]
    if len(used) >= 2:
        lo = np.array([cm[str(L)][0] for L in used])
        hi = np.array([cm[str(L)][-1] for L in used])
        merge = float((lo.max() - lo.min()) / max(abs(lo.mean()), 1e-12))
        splay = bool(hi[used.index(max(used))] < hi[used.index(min(used))])
    else:
        merge, splay = math.inf, False
    se_max = max(area["curve_se_max"].values()) if used else math.inf
    se_bound = 0.02 * 2.0 * t_grid[-1] / math.pi
    g = {"PASS_INPUT_COMPLETE": bool(shapes_ok and finite
                                     and not prod.get("unmeasured_L")),
         "PASS_T0_GEOMETRY_ORACLE": bool(oracle_ok),
         "PASS_LOW_T_MERGE": bool(merge < 0.08),
         "PASS_HIGH_T_SPLAY": splay,
         "PASS_CURVE_SE_BOUNDED": bool(se_max < se_bound),
         "PASS_MIN_PAIRS": bool(area["estimates"]["n_pairs_with_crossing"]
                                >= PHY051_MIN_PAIRS)}
    return {"gates": g, "overall_pass": all(g.values()),
            "diagnostics": {"low_t_max_rel_spread": merge,
                            "curve_se_max": se_max, "curve_se_bound": se_bound,
                            "band": list(plan["band"])}}


def report_stem(lattice: str) -> str:
    return f"260928 PHY051 {lattice} area-helicity report"


def build_report(prod: dict) -> dict:
    lattice = prod["lattice"]
    plan = PHY051_PLAN[lattice]
    area = analyse(prod, "ups_area")
    site = analyse(prod, "ups_site")
    oracle = geometry_oracle(lattice, prod["ladder"])
    g = gates(prod, area, oracle)
    label = o1_label(area["verdict"]["verdict"], site["verdict"]["verdict"])
    return {"module": "PHY051_area_helicity_crosscheck",
            "attribution": "Coworker Research / Coworkerz",
            "date": "2026-09-28", "lattice": lattice,
            "spec": "spec/260928 PHI HEX phy051 triangular kagome per-area "
                    "preregistration v01.md",
            "normalization_primary": PHY051_NORMALIZATION,
            "band": list(plan["band"]),
            "references": {k: plan[k] for k in
                           ("ref_free", "ref_free_alt", "ref_helicity")},
            "production": {**{k: prod[k] for k in (
                "t_grid", "ladder", "n_seeds", "n_therm", "n_meas",
                "unmeasured_L", "wall_s", "cpu_s", "max_workers",
                "have_numba")},
                "run_utc_start": prod.get("run_utc_start"),
                "run_utc_end": prod.get("run_utc_end")},
            "geometry_oracle": oracle,
            "primary_per_area": area,
            "convention_crosscheck_per_site": site,
            "o1_label": label,
            "pass_gates": g["gates"], "overall_pass": g["overall_pass"],
            "gate_diagnostics": g["diagnostics"],
            "claim_ceiling": "FINDING relative to preregistered band; no "
                             "best value; cross-family review before any "
                             "external statement",
            "data": {"ups_area": prod["ups_area"],
                     "ups_site": prod["ups_site"]}}


def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(x) for x in o]
    if isinstance(o, (np.floating, float)):
        v = float(o)
        return v if math.isfinite(v) else None
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def _fmt(x, f=".4f"):
    return "null" if x is None else format(x, f)


def write_text_report(rep: dict, path: Path) -> None:
    p = rep["production"]
    lines = [f"PHY051 - {rep['lattice']}: T_BKT mit Upsilon pro Flaeche "
             "(Wolff, Paar-Crossings) - O1-Diskriminator",
             "Coworker Research / Coworkerz, 2026-09-28", "=" * 70,
             f"Spec: {rep['spec']}",
             f"Leiter {p['ladder']}, T-Gitter {p['t_grid'][0]}..{p['t_grid'][-1]}"
             f" ({len(p['t_grid'])} Punkte), {p['n_seeds']} Seeds x "
             f"({p['n_therm']}+{p['n_meas']}) Sweeps, numba={p['have_numba']}",
             f"Wall {p['wall_s']:.0f} s auf {p['max_workers']} Prozessen, "
             f"CPU {p['cpu_s'] / 3600:.2f} h; unmeasured_L={p['unmeasured_L']}",
             f"Lauf (UTC): {p.get('run_utc_start')} .. {p.get('run_utc_end')}",
             ""]
    for title, ch in (("PRIMAER: Upsilon pro Flaeche", rep["primary_per_area"]),
                      ("QUERCHECK: per Site (Konvention, NICHT entscheidend)",
                       rep["convention_crosscheck_per_site"])):
        e = ch["estimates"]
        lines += [f"--- {title} ---",
                  "  Paar-Crossings T*(L1,L2) +- Jackknife-SE:"]
        for k, T in e["pairs"].items():
            lines.append(f"    ({k.replace('_', ',')}): {_fmt(T)} "
                         f"+- {_fmt(ch['pair_jackknife_se'].get(k))}")
        pr = e["primary_weighted_large_pairs"]
        lines += [f"  PRIMAER gewichtetes Mittel grosser Paare: {_fmt(ch['T_P'])}"
                  f" (n={pr['n_pairs']}, weighted={pr['weighted']})",
                  f"  v1 alle Paare gewichtet: {_fmt(e['v1_weighted_all_pairs']['T'])}",
                  f"  v2 groesstes Paar: {_fmt(e['v2_largest_pair'])}",
                  f"  v3 HKS 3-Par: {_fmt((e['v3_hks_3par'] or {}).get('T_c'))}"
                  f" (b={_fmt((e['v3_hks_3par'] or {}).get('b'), '.2f')})",
                  f"  v4 HKS b=1: {_fmt((e['v4_hks_b1'] or {}).get('T_c'))}",
                  f"  v5 WM freies C: {_fmt((e['v5_wm_free_c'] or {}).get('T'))}",
                  f"  sigma_sampler (Jackknife) = {_fmt(ch['sigma_sampler'])}",
                  f"  sigma_FSS (halbe Spannweite, vor Floor) = "
                  f"{_fmt(ch['sigma_fss_raw'])}",
                  f"  Verdikt vs Band {rep['band']}: {ch['verdict']}",
                  f"  Diagnostik ohne v3 (NICHT entscheidend): sigma_FSS = "
                  f"{_fmt(ch['diagnostic_without_hks3par']['sigma_fss_raw'])}"
                  f" -> {ch['diagnostic_without_hks3par']['verdict'].get('verdict')}"
                  f" ({ch['diagnostic_without_hks3par']['verdict'].get('rule')})",
                  f"  z vs ref_free = {_fmt(ch['z_vs_ref_free'], '.2f')}; "
                  f"z vs ref_free_alt = {_fmt(ch['z_vs_ref_free_alt'], '.2f')}; "
                  f"z vs ref_helicity = {_fmt(ch['z_vs_ref_helicity'], '.2f')}",
                  ""]
    lines += ["--- O1-Diskriminator (Spec Abschnitt 6) ---",
              f"  {rep['o1_label']}", "",
              "--- Gates (Pipeline-Integritaet, keine Physik-Wahrheit) ---"]
    for g, v in rep["pass_gates"].items():
        lines.append(f"  [{'PASS' if v else 'FAIL'}] {g}")
    lines += [f"  OVERALL: {'PASS' if rep['overall_pass'] else 'FAIL'}", "",
              f"Referenzen: {json.dumps(rep['references'])}", "",
              "Grenzen: FINDING relativ zum vorregistrierten Band, kein "
              "Bestwert;", "externe Aussage erst nach Cross-Family-Review."]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_lattice(lattice: str, max_workers: int = 4,
                out_dir: Path | None = None) -> dict:
    prod = produce(lattice, max_workers=max_workers)
    rep = _clean(build_report(prod))
    res = out_dir or (_ROOT / "results")
    stem = report_stem(lattice)
    (res / f"{stem}.json").write_text(
        json.dumps(rep, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    write_text_report(rep, res / f"{stem}.txt")
    print(f"Report: results/{stem}.{{json,txt}}; O1={rep['o1_label']}; "
          f"per-area {rep['primary_per_area']['verdict']['verdict']}; "
          f"gates {'PASS' if rep['overall_pass'] else 'FAIL'}", flush=True)
    return rep


def prod_from_report(rep: dict) -> dict:
    """Produktions-Dict (Rohdaten + Metadaten) aus einem Report-JSON."""
    pr = rep["production"]
    return {"lattice": rep["lattice"], "t_grid": pr["t_grid"],
            "ladder": pr["ladder"], "n_seeds": pr["n_seeds"],
            "n_therm": pr["n_therm"], "n_meas": pr["n_meas"],
            "unmeasured_L": pr["unmeasured_L"], "wall_s": pr["wall_s"],
            "cpu_s": pr["cpu_s"], "max_workers": pr["max_workers"],
            "have_numba": pr["have_numba"],
            "run_utc_start": pr.get("run_utc_start"),
            "run_utc_end": pr.get("run_utc_end"),
            "ups_area": rep["data"]["ups_area"],
            "ups_site": rep["data"]["ups_site"]}


def regenerate(lattice: str, out_dir: Path | None = None) -> dict:
    """Report MC-frei aus den committeten Rohdaten neu ableiten (Auswertung
    ist eine deterministische Funktion der Rohdaten; Rohdaten unveraendert)."""
    res = out_dir or (_ROOT / "results")
    stem = report_stem(lattice)
    old = json.loads((res / f"{stem}.json").read_text(encoding="utf-8"))
    rep = _clean(build_report(prod_from_report(old)))
    (res / f"{stem}.json").write_text(
        json.dumps(rep, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    write_text_report(rep, res / f"{stem}.txt")
    print(f"Regenerated: results/{stem}.{{json,txt}}; O1={rep['o1_label']}",
          flush=True)
    return rep


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which == "regenerate":
        for name in (sys.argv[2:] or ("triangular", "kagome")):
            regenerate(name)
        sys.exit(0)
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    names = ("triangular", "kagome") if which == "all" else (which,)
    for name in names:
        run_lattice(name, max_workers=workers)
