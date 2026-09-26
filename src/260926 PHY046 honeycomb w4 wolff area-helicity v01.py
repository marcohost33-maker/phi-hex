"""PHY046 - W4 v02: honeycomb T_BKT aus Upsilon PRO FLAECHE, kanonischer
Wolff (Numba), HKS-Paar-Extrapolation - exakt nach Vorregistrierung.

Coworker Research / Coworkerz, 26. September 2026
de-CH Konventionen, ASCII-Quotes, ss statt Eszett, keine Personennamen.

Vertrag: `spec/260926 PHI HEX w4 honeycomb preregistration v02.md` (vor jeder
W4-Datennahme committet). Normierung: `spec/260926 PHI HEX O1 helicity
normalization decision v01.md`. Kernel: PHY045 (validiert: T->0-Orakel,
Quercheck gegen PHY031-Python-Wolff).

Ablauf: je (L, T, Seed) ein Wolff-Lauf -> per-Seed-Upsilon_A; Seed-Mittel-
kurven je L; C-eliminierte WM-Paare (L, 2L); HKS-Extrapolation
T*(L) = T_c + a / ln^2(b L); FSS-Varianten v1..v4; Jackknife ueber Seeds;
Entscheid via `w4_verdict` (PHY044). Kein Parameter dieses Moduls darf nach
Sicht der W4-Daten geaendert werden (sonst: v03 als post-hoc-Nachtrag).

EHRLICHKEIT: das Ergebnis ist ein FINDING relativ zum vorregistrierten Band,
kein Bestwert; externe Aussagen erst nach Cross-Family-Review.
"""
from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import json  # noqa: E402
import math  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import importlib.util  # noqa: E402
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
_p44 = _load("phy044_wl_calibration",
             "260926 PHY044 honeycomb wl calibration v01.py")
_p40 = _load("phy040_wang_landau",
             "260616 PHY040 wang-landau entropic helicity v01.py")

tbkt_pair_from_curves = _p40.tbkt_pair_from_curves
w4_verdict = _p44.w4_verdict
W4_BAND = _p44.W4_BAND

# ---- vorregistrierte Konstanten (Spec v02; Test-gebunden) -----------------
W4V2_NORMALIZATION = "per_area"
W4V2_LADDER = (32, 48, 64, 96, 128, 192, 256)
W4V2_T_GRID = tuple(round(0.545 + 0.005 * k, 3) for k in range(16))
W4V2_N_SEEDS = 8
W4V2_N_THERM = 500
W4V2_N_MEAS = 1500
W4V2_SEED_BASE = 47_000_000
W4V2_HKS_B_RANGE = (0.1, 10.0)
W4V2_WALL_BUDGET_H = 24.0

AREA_PER_SITE = _p45.AREA_PER_SITE["honeycomb"]


def seed_for(L: int, t_idx: int, s: int) -> int:
    """Seed-Vertrag v02: 47_000_000 + 1000*L + 10*t_idx + s."""
    return W4V2_SEED_BASE + 1000 * L + 10 * t_idx + s


def hks_pairs(ladder) -> list:
    """(L, 2L)-Paare innerhalb der Leiter."""
    s = set(ladder)
    return [(L, 2 * L) for L in sorted(ladder) if 2 * L in s]


# ============================================================================
# 1. Produktion
# ============================================================================

def _job(args: tuple) -> dict:
    L, t_idx, T, s, n_therm, n_meas = args
    lat = _p45.build("honeycomb", L)
    nbr, deg = _p45._nbr_arrays(lat)
    th = np.zeros(lat.n)
    run = _p45._nb_run if _p45.HAVE_NUMBA else _p45._py_run
    t0 = time.perf_counter()
    d = run(th, nbr, deg, lat.ei, lat.ej,
            np.ascontiguousarray(lat.disp[:, 0]),
            np.ascontiguousarray(lat.disp[:, 1]),
            1.0 / T, n_therm, n_meas, seed_for(L, t_idx, s))
    beta = 1.0 / T
    num = 0.5 * ((d[:, 1].mean() - beta * np.mean(d[:, 2] ** 2))
                 + (d[:, 3].mean() - beta * np.mean(d[:, 4] ** 2)))
    return {"L": L, "t_idx": t_idx, "s": s,
            "ups_area": float(num / lat.torus_area),
            "ups_site": float(num / lat.n),
            "wall_s": time.perf_counter() - t0}


def produce(ladder=W4V2_LADDER, t_grid=W4V2_T_GRID, n_seeds=W4V2_N_SEEDS,
            n_therm=W4V2_N_THERM, n_meas=W4V2_N_MEAS,
            max_workers=4) -> dict:
    jobs = [(L, k, T, s, n_therm, n_meas) for L in ladder
            for k, T in enumerate(t_grid) for s in range(n_seeds)]
    jobs.sort(key=lambda j: (-j[0], j[1], j[3]))
    t0 = time.perf_counter()
    if max_workers <= 1:
        res = [_job(j) for j in jobs]
    else:
        with ProcessPoolExecutor(max_workers=max_workers) as ex:
            res = list(ex.map(_job, jobs, chunksize=1))
    area = {L: np.zeros((len(t_grid), n_seeds)) for L in ladder}
    site = {L: np.zeros((len(t_grid), n_seeds)) for L in ladder}
    cpu = 0.0
    for r in res:
        area[r["L"]][r["t_idx"], r["s"]] = r["ups_area"]
        site[r["L"]][r["t_idx"], r["s"]] = r["ups_site"]
        cpu += r["wall_s"]
    return {"t_grid": list(t_grid), "ladder": list(ladder),
            "n_seeds": n_seeds, "n_therm": n_therm, "n_meas": n_meas,
            "ups_area": {str(L): area[L].tolist() for L in ladder},
            "ups_site": {str(L): site[L].tolist() for L in ladder},
            "wall_s": time.perf_counter() - t0, "cpu_s": cpu,
            "max_workers": max_workers}


# ============================================================================
# 2. Auswertung (rein, deterministisch aus den Produktionsdaten)
# ============================================================================

def pair_crossings(curves: dict, t_grid, ladder) -> dict:
    t = np.asarray(t_grid, dtype=float)
    return {(a, b): tbkt_pair_from_curves(t, curves[a], a, curves[b], b)
            for a, b in hks_pairs(ladder)}


def fit_hks(Ls, Ts, sig=None, b_fixed=None, b_range=W4V2_HKS_B_RANGE):
    """T*(L) = T_c + a / ln^2(b L); b-Scan (log) + gewichtete lineare LSQ.
    Rueckgabe dict oder None (zu wenige Punkte)."""
    Ls = np.asarray(Ls, dtype=float)
    Ts = np.asarray(Ts, dtype=float)
    n_par = 2 if b_fixed is not None else 3
    if len(Ls) < n_par:
        return None
    w = (np.ones_like(Ts) if sig is None
         else 1.0 / np.maximum(np.asarray(sig, dtype=float), 1e-6) ** 2)
    bs = ([b_fixed] if b_fixed is not None
          else np.exp(np.linspace(math.log(b_range[0]),
                                  math.log(b_range[1]), 400)))
    best = None
    for b in bs:
        x = 1.0 / np.log(b * Ls) ** 2
        if np.any(~np.isfinite(x)) or np.any(b * Ls <= 1.0):
            continue
        X = np.column_stack([np.ones_like(x), x])
        W = np.diag(w)
        try:
            cov = np.linalg.inv(X.T @ W @ X)
        except np.linalg.LinAlgError:
            continue
        p = cov @ X.T @ W @ Ts
        chi2 = float((Ts - X @ p) @ W @ (Ts - X @ p))
        if best is None or chi2 < best["chi2"]:
            best = {"T_c": float(p[0]), "a": float(p[1]), "b": float(b),
                    "chi2": chi2, "dof": len(Ls) - n_par}
    return best


def wm_free_c_fit(curves_mean, curves_se, t_grid, ladder) -> dict | None:
    """v4: je T WM-Fit mit freiem C ueber alle L; T am chi^2-Minimum
    (parabolisch verfeinert). Upsilon = (2T/pi)(1 + 1/(2 ln L + C))."""
    t = np.asarray(t_grid, dtype=float)
    Cs = np.linspace(-20.0, 60.0, 1601)
    chi = []
    for k, T in enumerate(t):
        y = np.array([curves_mean[L][k] for L in ladder])
        e = np.array([max(curves_se[L][k], 1e-6) for L in ladder])
        lnL = np.log(np.asarray(ladder, dtype=float))
        best = math.inf
        for C in Cs:
            den = 2.0 * lnL + C
            if np.any(den <= 0.05):
                continue
            m = (2.0 * T / math.pi) * (1.0 + 1.0 / den)
            best = min(best, float(np.sum(((y - m) / e) ** 2)))
        chi.append(best)
    chi = np.asarray(chi)
    if not np.any(np.isfinite(chi)):
        return None
    k = int(np.nanargmin(chi))
    T_min = float(t[k])
    if 0 < k < len(t) - 1 and all(np.isfinite(chi[k - 1:k + 2])):
        c0, c1, c2 = chi[k - 1], chi[k], chi[k + 1]
        den = c0 - 2 * c1 + c2
        if den > 0:
            T_min = float(t[k] + 0.5 * (t[1] - t[0]) * (c0 - c2) / den)
    return {"T": T_min, "chi2_min": float(chi[k]),
            "at_grid_edge": bool(k in (0, len(t) - 1))}


def estimate(curves: dict, curves_se: dict, t_grid, ladder,
             pair_sig=None) -> dict:
    """Primaer + Varianten aus EINEM Satz Kurven."""
    pc = pair_crossings(curves, t_grid, ladder)
    pts = [(a, T) for (a, b), T in pc.items() if T is not None]
    Ls = [a for a, _ in pts]
    Ts = [T for _, T in pts]
    sig = None
    if pair_sig is not None:
        sig = [pair_sig.get((a, 2 * a)) for a in Ls]
        # fehlt fuer ein Paar die Jackknife-SE (oder ist sie 0), wird
        # ungewichtet gefittet - kein Mischen von echten und Ersatz-Gewichten
        if any(x is None or x <= 0.0 for x in sig):
            sig = None
    prim = fit_hks(Ls, Ts, sig)
    v1 = fit_hks(Ls, Ts, sig, b_fixed=1.0)
    v2 = fit_hks(Ls[1:], Ts[1:], None if sig is None else sig[1:]) \
        if len(Ls) >= 4 else None
    v3 = Ts[-1] if Ts else None
    v4 = wm_free_c_fit(curves, curves_se, t_grid, ladder)
    return {"pairs": {f"{a}_{b}": T for (a, b), T in pc.items()},
            "n_pairs_with_crossing": len(pts),
            "primary_hks": prim, "v1_hks_b1": v1, "v2_hks_drop_smallest": v2,
            "v3_largest_pair": v3, "v4_wm_free_c": v4}


def _jk_se(vals: list) -> float | None:
    v = [x for x in vals if x is not None]
    n = len(v)
    if n < 2:
        return None
    m = float(np.mean(v))
    return float(math.sqrt((n - 1) / n * sum((x - m) ** 2 for x in v)))


def analyse(prod: dict, key: str = "ups_area") -> dict:
    """Vollstaendige vorregistrierte Auswertung inkl. Jackknife + Verdikt."""
    t_grid = prod["t_grid"]
    ladder = [int(x) for x in prod["ladder"]]
    data = {L: np.asarray(prod[key][str(L)]) for L in ladder}
    ns = prod["n_seeds"]
    mean = {L: data[L].mean(axis=1) for L in ladder}
    se = {L: data[L].std(axis=1, ddof=1) / math.sqrt(ns) for L in ladder}
    jk_curves = [{L: np.delete(data[L], s, axis=1).mean(axis=1)
                  for L in ladder} for s in range(ns)]
    # Paar-Jackknife-SE (Gewichte fuer HKS)
    pc_full = pair_crossings(mean, t_grid, ladder)
    pc_jk = [pair_crossings(c, t_grid, ladder) for c in jk_curves]
    pair_sig = {p: _jk_se([j[p] for j in pc_jk]) for p in pc_full}
    full = estimate(mean, se, t_grid, ladder, pair_sig)
    # Jackknife des Primaer-Schaetzers (Gewichte fest aus dem Vollsatz)
    jk_T = []
    for c in jk_curves:
        e = estimate(c, se, t_grid, ladder, pair_sig)
        jk_T.append(None if e["primary_hks"] is None
                    else e["primary_hks"]["T_c"])
    T_w4 = None if full["primary_hks"] is None else full["primary_hks"]["T_c"]
    variants = [T_w4,
                None if full["v1_hks_b1"] is None else full["v1_hks_b1"]["T_c"],
                None if full["v2_hks_drop_smallest"] is None
                else full["v2_hks_drop_smallest"]["T_c"],
                full["v3_largest_pair"],
                None if full["v4_wm_free_c"] is None
                else full["v4_wm_free_c"]["T"]]
    vv = [v for v in variants if v is not None]
    sigma_fss = 0.5 * (max(vv) - min(vv)) if len(vv) >= 2 else None
    sigma_sampler = _jk_se(jk_T)
    verdict = w4_verdict(T_w4, sigma_sampler, sigma_fss,
                         full["n_pairs_with_crossing"])
    return {"normalization": key, "estimates": full,
            "pair_jackknife_se": {f"{a}_{b}": v
                                  for (a, b), v in pair_sig.items()},
            "T_w4": T_w4, "sigma_sampler": sigma_sampler,
            "sigma_fss_raw": sigma_fss, "variants": variants,
            "verdict": verdict,
            "curve_se_max": {str(L): float(np.max(se[L])) for L in ladder}}


REPORT_STEM = "260926 PHY046 honeycomb w4 wolff area-helicity report"


def build_report(prod: dict) -> dict:
    area = analyse(prod, "ups_area")
    site = analyse(prod, "ups_site")   # nur Konventions-Quercheck
    return {"module": "PHY046_honeycomb_w4_v02",
            "attribution": "Coworker Research / Coworkerz",
            "date": "2026-09-26",
            "spec": "spec/260926 PHI HEX w4 honeycomb preregistration v02.md",
            "band_B": list(W4_BAND),
            "production": {k: prod[k] for k in (
                "t_grid", "ladder", "n_seeds", "n_therm", "n_meas",
                "wall_s", "cpu_s", "max_workers")},
            "primary_per_area": area,
            "convention_crosscheck_per_site": site,
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
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def _fmt(x, f=".4f"):
    return "null" if x is None else format(x, f)


def write_text_report(rep: dict, path: Path) -> None:
    a = rep["primary_per_area"]
    s = rep["convention_crosscheck_per_site"]
    p = rep["production"]
    L_ = ["PHY046 - W4 v02: honeycomb T_BKT (Upsilon pro Flaeche, Wolff, HKS)",
          "Coworker Research / Coworkerz, 2026-09-26", "=" * 70,
          f"Spec: {rep['spec']}",
          f"Leiter {p['ladder']}, T-Gitter {p['t_grid'][0]}..{p['t_grid'][-1]} "
          f"({len(p['t_grid'])} Punkte), {p['n_seeds']} Seeds x "
          f"({p['n_therm']}+{p['n_meas']}) Sweeps",
          f"Wall {p['wall_s']:.0f} s auf {p['max_workers']} Prozessen, "
          f"CPU {p['cpu_s'] / 3600:.2f} h", "",
          "--- Paar-Crossings T*(L,2L) (Upsilon pro Flaeche) ---"]
    for k, T in a["estimates"]["pairs"].items():
        L_.append(f"  ({k.replace('_', ',')}): {_fmt(T)} "
                  f"+- {_fmt(a['pair_jackknife_se'].get(k))}")
    e = a["estimates"]
    L_ += ["", "--- Schaetzer ---",
           f"  PRIMAER HKS (T_c, a, b frei): {_fmt(a['T_w4'])} "
           f"(chi2 {_fmt((e['primary_hks'] or {}).get('chi2'), '.2f')})",
           f"  v1 HKS b=1: {_fmt((e['v1_hks_b1'] or {}).get('T_c'))}",
           f"  v2 HKS ohne kleinstes Paar: "
           f"{_fmt((e['v2_hks_drop_smallest'] or {}).get('T_c'))}",
           f"  v3 groesstes Paar: {_fmt(e['v3_largest_pair'])}",
           f"  v4 WM freies C: {_fmt((e['v4_wm_free_c'] or {}).get('T'))}",
           f"  sigma_sampler (Jackknife) = {_fmt(a['sigma_sampler'])}",
           f"  sigma_FSS (halbe Spannweite, vor Floor) = "
           f"{_fmt(a['sigma_fss_raw'])}", "",
           "--- Vorregistriertes Verdikt (w4_verdict) ---",
           f"  Band B = {rep['band_B']}",
           f"  -> {a['verdict']}", "",
           "--- Konventions-Quercheck (per Site; NICHT entscheidend) ---",
           f"  T_W4(per Site) = {_fmt(s['T_w4'])}; Paare "
           f"{ {k: _fmt(v) for k, v in s['estimates']['pairs'].items()} }", "",
           "Grenzen: FINDING relativ zum vorregistrierten Band, kein Bestwert;",
           "externe Aussage erst nach Cross-Family-Review."]
    path.write_text("\n".join(L_) + "\n", encoding="utf-8")


if __name__ == "__main__":
    out_json = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    prod = produce()
    rep = _clean(build_report(prod))
    res = _ROOT / "results"
    (res / f"{REPORT_STEM}.json").write_text(
        json.dumps(rep, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    write_text_report(rep, res / f"{REPORT_STEM}.txt")
    if out_json:
        out_json.write_text(json.dumps(_clean(prod)), encoding="utf-8")
    print(f"Report: results/{REPORT_STEM}.{{json,txt}}; "
          f"Verdikt {rep['primary_per_area']['verdict']}")
