"""PHY047 - Neuauswertung COMMITTETER Helicity-Daten mit Upsilon pro Flaeche
(Folge von Audit O1; Lineage: die Original-Reports bleiben unveraendert).

Coworker Research / Coworkerz, 26. September 2026
de-CH Konventionen, ASCII-Quotes, ss statt Eszett, keine Personennamen.

Reihenfolge-Vertrag: ausgefuehrt NACH dem Commit der W4-v02-Vorregistrierung
(82e2738) - die Neuauswertung alter Daten kann das W4-Protokoll nicht mehr
beeinflussen.

Quellen (committet, SHA-gepinnt):
- honeycomb PHY042 (WL, L=24/32/48, Walker-Mittelkurven, JSON)
- honeycomb PHY032 (Wolff, L=12/24/48, Tabelle Upsilon(T,L) +- SEM)
- kagome   PHY033 (Wolff, L=12/24/36, Tabelle)
- triangular PHY030 v02 (Wolff, L=9/13/19, T <= 1.46) - pro Flaeche liegt das
  NK-Crossing oberhalb des Gitters: keine Schaetzung moeglich (ehrlich
  ausgewiesen).

Methoden identisch zu den Originalen (C-eliminierte WM-Paare nach PHY040,
per-L-NK-Crossing, WM-Fit mit freiem C nach PHY046 v4); einziger Unterschied
ist Upsilon_area = Upsilon_site / a_s. Beide Normierungen werden nebeneinander
berichtet; die per-Site-Spalte muss die Original-Reports reproduzieren
(Drift-Guard).

EHRLICHKEIT: kleine L (<= 48), keine HKS-Extrapolation -> FINDING zur Groesse
und Richtung des Normierungseffekts, kein Bestwert.
"""
from __future__ import annotations

import json
import math
import re
import sys
import importlib.util
from pathlib import Path

import numpy as np

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


_p40 = _load("phy040_wang_landau",
             "260616 PHY040 wang-landau entropic helicity v01.py")
_p45 = _load("phy045_normalization",
             "260926 PHY045 helicity normalization O1 test v01.py")
_p46 = _load("phy046_w4_v02",
             "260926 PHY046 honeycomb w4 wolff area-helicity v01.py")

tbkt_pair = _p40.tbkt_pair_from_curves
A_S = _p45.AREA_PER_SITE

RES = _ROOT / "results"
SOURCES = {
    "honeycomb_phy042": RES / "260707 PHY042 honeycomb wl-fss L24-32-48 gate report.json",
    "honeycomb_phy032": RES / "260607 PHY032 honeycomb wm-logfit bootstrap report.txt",
    "kagome_phy033": RES / "260609 PHY033 kagome tbkt per-site report.txt",
    "triangular_phy030v02": RES / "260604 PHY030 v02 wm-logfit report.txt",
}

_ROW = re.compile(r"^\s*L=\s*(\d+):\s*(.*)$")
_CELL = re.compile(r"T=([0-9.]+):([0-9.]+)\+-([0-9.]+)")


def parse_table(path: Path) -> tuple[list, dict, dict]:
    """Tabelle 'L=..: T=..:u+-sem ...' aus dem Abschnitt '--- Messung'."""
    text = path.read_text(encoding="utf-8")
    after = text.split("--- Messung", 1)[1].split("\n", 1)[1]
    sec = after.split("\n---", 1)[0]
    data: dict[int, dict[float, tuple[float, float]]] = {}
    for line in sec.splitlines():
        m = _ROW.match(line)
        if not m:
            continue
        L = int(m.group(1))
        data[L] = {float(t): (float(u), float(e))
                   for t, u, e in _CELL.findall(m.group(2))}
    Ls = sorted(data)
    temps = sorted(set.intersection(*(set(d) for d in data.values())))
    mean = {L: np.array([data[L][t][0] for t in temps]) for L in Ls}
    sem = {L: np.array([data[L][t][1] for t in temps]) for L in Ls}
    return temps, mean, sem


def nk_crossing(temps, y) -> float | None:
    """Erstes Vorzeichenwechsel-Crossing Upsilon(T) = 2T/pi (benachbart)."""
    d = [yy - 2.0 * t / math.pi for t, yy in zip(temps, y)]
    for i in range(len(d) - 1):
        if d[i] == 0:
            return float(temps[i])
        if d[i] > 0 > d[i + 1]:
            return float(temps[i] + (temps[i + 1] - temps[i])
                         * d[i] / (d[i] - d[i + 1]))
    return None


def analyse_set(temps, mean_site, sem_site, a_s: float) -> dict:
    """Beide Normierungen: per-L-NK, alle Paare, WM-Fit freies C."""
    out = {}
    Ls = sorted(mean_site)
    for norm, fac in (("per_site", 1.0), ("per_area", 1.0 / a_s)):
        m = {L: mean_site[L] * fac for L in Ls}
        e = {L: sem_site[L] * fac for L in Ls}
        pairs = {f"{a}_{b}": tbkt_pair(np.asarray(temps), m[a], a, m[b], b)
                 for i, a in enumerate(Ls) for b in Ls[i + 1:]}
        wm = _p46.wm_free_c_fit(m, e, temps, Ls)
        out[norm] = {"nk_per_L": {str(L): nk_crossing(temps, m[L])
                                  for L in Ls},
                     "pairs": pairs, "wm_free_c": wm,
                     "grid": [temps[0], temps[-1]]}
    return out


def run() -> dict:
    rep = {"module": "PHY047_per_area_reanalysis",
           "attribution": "Coworker Research / Coworkerz",
           "date": "2026-09-26",
           "decision": "spec/260926 PHI HEX O1 helicity normalization "
                       "decision v01.md",
           "order_contract": "executed after W4 v02 pre-registration "
                             "commit 82e2738",
           "sets": {}}
    # honeycomb PHY042 (WL-Mittelkurven, keine SEM -> Walker-Spread/2)
    j = json.loads(SOURCES["honeycomb_phy042"].read_text(encoding="utf-8"))
    temps = j["main_curves"]["24"]["T"]
    mean = {int(L): np.asarray(c["y2"]) for L, c in j["main_curves"].items()}
    sem = {}
    for L in mean:
        nw = int(j["n_walkers"][str(L)])
        ys = np.array([j["curves"][f"{L}_{w}"]["y2"] for w in range(nw)])
        sem[L] = (ys.std(axis=0, ddof=1) / math.sqrt(nw) if nw > 1
                  else np.full(len(temps), 0.01))
    rep["sets"]["honeycomb_phy042_wl"] = {
        "a_s": A_S["honeycomb"], "committed_pairs_per_site":
            j["pair_tbkt_mean_curves"],
        **analyse_set(temps, mean, sem, A_S["honeycomb"])}
    for key, lat in (("honeycomb_phy032", "honeycomb"),
                     ("kagome_phy033", "kagome"),
                     ("triangular_phy030v02", "triangular")):
        t, m, e = parse_table(SOURCES[key])
        rep["sets"][key] = {"a_s": A_S[lat], **analyse_set(t, m, e, A_S[lat])}
    return rep


def _fmt(x):
    return "ausserhalb Gitter" if x is None else f"{x:.4f}"


def write_text(rep: dict, path: Path) -> None:
    L_ = ["PHY047 - Neuauswertung committeter Daten mit Upsilon pro Flaeche",
          "Coworker Research / Coworkerz, 2026-09-26", "=" * 70,
          f"Entscheid: {rep['decision']}",
          f"Reihenfolge: {rep['order_contract']}",
          "Methoden identisch zu den Originalen; einzig Upsilon_area = "
          "Upsilon_site / a_s.", ""]
    for name, s in rep["sets"].items():
        L_.append(f"--- {name} (a_s = {s['a_s']:.4f}, T-Gitter "
                  f"{s['per_site']['grid'][0]}..{s['per_site']['grid'][1]}) ---")
        for norm in ("per_site", "per_area"):
            r = s[norm]
            nk = ", ".join(f"L{L}={_fmt(v)}" for L, v in r["nk_per_L"].items())
            pr = ", ".join(f"({k.replace('_', ',')})={_fmt(v)}"
                           for k, v in r["pairs"].items())
            wm = r["wm_free_c"]
            wm_t = "n/a" if wm is None else (
                f"{wm['T']:.4f}" + (" [Gitterrand]" if wm["at_grid_edge"]
                                    else ""))
            L_.append(f"  {norm:<9} NK: {nk}")
            L_.append(f"  {'':<9} Paare: {pr}")
            L_.append(f"  {'':<9} WM freies C: {wm_t}")
        L_.append("")
    L_.append("Grenzen: kleine L (<= 48), keine HKS-Extrapolation; FINDING zur "
              "Groesse und Richtung des Normierungseffekts, kein Bestwert.")
    path.write_text("\n".join(L_) + "\n", encoding="utf-8")


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


STEM = "260926 PHY047 per-area reanalysis committed report"

if __name__ == "__main__":
    rep = _clean(run())
    (RES / f"{STEM}.json").write_text(
        json.dumps(rep, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    write_text(rep, RES / f"{STEM}.txt")
    print((RES / f"{STEM}.txt").read_text(encoding="utf-8"))
