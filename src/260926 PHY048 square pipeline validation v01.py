"""PHY048 - POST-HOC V&V: dieselbe W4-v02-Pipeline (Wolff, Upsilon pro
Flaeche, HKS-Paare, Varianten, Jackknife) auf dem QUADRATGITTER, dessen
T_BKT hochpraezise bekannt ist (0.8929-0.8935; HKS arXiv:1302.2900 0.8935(1),
TRG/Tensor-Netz ~0.89294). Zweck: den Methoden-Bias der Pipeline bei gleicher
Leiter (32..256) messen, um das W4-Ergebnis einordnen zu koennen.

Coworker Research / Coworkerz, 26. September 2026
de-CH Konventionen, ASCII-Quotes, ss statt Eszett, keine Personennamen.

POST-HOC-Kennzeichnung: dieser Lauf war NICHT vorregistriert; er wurde nach
Sicht des W4-Ergebnisses (PHY046) aufgesetzt. Er aendert das vorregistrierte
W4-Verdikt nicht und wird getrennt berichtet. Auf dem square-Gitter sind
per Site und per Flaeche identisch (a_s = 1).
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
_p46 = _load("phy046_w4_v02",
             "260926 PHY046 honeycomb w4 wolff area-helicity v01.py")

LADDER = _p46.W4V2_LADDER                      # identische Leiter
T_GRID = tuple(round(0.860 + 0.005 * k, 3) for k in range(16))
N_SEEDS, N_THERM, N_MEAS = (_p46.W4V2_N_SEEDS, _p46.W4V2_N_THERM,
                            _p46.W4V2_N_MEAS)
SEED_BASE = 48_000_000
T_REF = {"HKS_2013": 0.8935, "TRG_2020": 0.89290, "level_spectroscopy":
         0.892943}


def _job(args: tuple) -> dict:
    L, t_idx, T, s = args
    lat = _p45.build("square", L)
    nbr, deg = _p45._nbr_arrays(lat)
    th = np.zeros(lat.n)
    run = _p45._nb_run if _p45.HAVE_NUMBA else _p45._py_run
    t0 = time.perf_counter()
    d = run(th, nbr, deg, lat.ei, lat.ej,
            np.ascontiguousarray(lat.disp[:, 0]),
            np.ascontiguousarray(lat.disp[:, 1]),
            1.0 / T, N_THERM, N_MEAS, SEED_BASE + 1000 * L + 10 * t_idx + s)
    beta = 1.0 / T
    num = 0.5 * ((d[:, 1].mean() - beta * np.mean(d[:, 2] ** 2))
                 + (d[:, 3].mean() - beta * np.mean(d[:, 4] ** 2)))
    return {"L": L, "t_idx": t_idx, "s": s,
            "ups_area": float(num / lat.torus_area),
            "wall_s": time.perf_counter() - t0}


def produce(max_workers: int = 4) -> dict:
    jobs = [(L, k, T, s) for L in LADDER for k, T in enumerate(T_GRID)
            for s in range(N_SEEDS)]
    jobs.sort(key=lambda j: (-j[0], j[1], j[3]))
    t0 = time.perf_counter()
    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        res = list(ex.map(_job, jobs, chunksize=1))
    area = {L: np.zeros((len(T_GRID), N_SEEDS)) for L in LADDER}
    cpu = 0.0
    for r in res:
        area[r["L"]][r["t_idx"], r["s"]] = r["ups_area"]
        cpu += r["wall_s"]
    return {"t_grid": list(T_GRID), "ladder": list(LADDER),
            "n_seeds": N_SEEDS, "n_therm": N_THERM, "n_meas": N_MEAS,
            "ups_area": {str(L): area[L].tolist() for L in LADDER},
            "wall_s": time.perf_counter() - t0, "cpu_s": cpu,
            "max_workers": max_workers}


STEM = "260926 PHY048 square pipeline validation report"


def build_report(prod: dict) -> dict:
    a = _p46.analyse(prod, "ups_area")
    bias = {k: (a["T_w4"] - v) / v if a["T_w4"] is not None else None
            for k, v in T_REF.items()}
    pair_bias = {k: (T - T_REF["TRG_2020"]) / T_REF["TRG_2020"]
                 if T is not None else None
                 for k, T in a["estimates"]["pairs"].items()}
    return {"module": "PHY048_square_pipeline_validation",
            "attribution": "Coworker Research / Coworkerz",
            "date": "2026-09-26", "post_hoc": True,
            "purpose": "method bias of the W4-v02 pipeline on the square "
                       "lattice at the same ladder",
            "references": T_REF,
            "production": {k: prod[k] for k in (
                "t_grid", "ladder", "n_seeds", "n_therm", "n_meas",
                "wall_s", "cpu_s", "max_workers")},
            "analysis": a, "relative_bias_T_w4": bias,
            "relative_bias_pairs_vs_TRG": pair_bias,
            "data": {"ups_area": prod["ups_area"]}}


def write_text(rep: dict, path: Path) -> None:
    a = rep["analysis"]
    f = _p46._fmt
    L_ = ["PHY048 - POST-HOC V&V der W4-v02-Pipeline auf dem Quadratgitter",
          "Coworker Research / Coworkerz, 2026-09-26", "=" * 70,
          "NICHT vorregistriert; aendert das W4-Verdikt nicht.",
          f"Referenzen: {rep['references']}",
          f"Leiter {rep['production']['ladder']}, T {rep['production']['t_grid'][0]}"
          f"..{rep['production']['t_grid'][-1]}, Wall "
          f"{rep['production']['wall_s']:.0f} s", "",
          "--- Paar-Crossings T*(L,2L) und Abweichung vs 0.89290 ---"]
    for k, T in a["estimates"]["pairs"].items():
        b = rep["relative_bias_pairs_vs_TRG"][k]
        L_.append(f"  ({k.replace('_', ',')}): {f(T)} +- "
                  f"{f(a['pair_jackknife_se'].get(k))}  "
                  f"({'n/a' if b is None else f'{100 * b:+.2f} %'})")
    L_ += ["", f"  HKS-Primaer T_c = {f(a['T_w4'])}; Varianten "
           f"{[None if v is None else round(v, 4) for v in a['variants']]}",
           f"  sigma_sampler = {f(a['sigma_sampler'])}, sigma_FSS roh = "
           f"{f(a['sigma_fss_raw'])}",
           "  relative Abweichung T_c: " + ", ".join(
               f"{k} {100 * v:+.2f} %" for k, v in
               rep["relative_bias_T_w4"].items() if v is not None)]
    path.write_text("\n".join(L_) + "\n", encoding="utf-8")


if __name__ == "__main__":
    rep = _p46._clean(build_report(produce()))
    res = _ROOT / "results"
    (res / f"{STEM}.json").write_text(
        json.dumps(rep, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    write_text(rep, res / f"{STEM}.txt")
    print((res / f"{STEM}.txt").read_text(encoding="utf-8"))
    sys.exit(0 if math.isfinite(rep["analysis"]["T_w4"] or float("nan"))
             else 1)
