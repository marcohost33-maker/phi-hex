"""PHY037 - HKS-L->inf-Extrapolation GITTER-UEBERGREIFEND (honeycomb / kagome /
triangular) auf den jeweils committed Paar-Schaetzern.

Generalisiert den HKS-Apparat (arXiv:1302.2900; PHY034) auf alle drei vermessenen
Gitter und zieht die EHRLICHE Methoden-Bilanz: wo konvergiert die C-eliminierte
Paar-Extrapolation sauber zur Referenz, wo nicht?

DATENQUELLE (committed Reports, zur LAUFZEIT geparst - single source of truth,
kein eingebettetes Daten-Duplikat):
  honeycomb : results/260607 PHY032 ... report.txt   (L=12,24,48)
  kagome    : results/260609 PHY033 ... report.txt   (L=12,24,36)
  triangular: results/260604 PHY030 v02 ... report.txt (L=9,13,19)
Geparst werden die "T_BKT(L=a,b) = X"-Paar-Zeilen (Methode B, C-frei).

BEFUND (ehrlich, NICHT getunt):
  - kagome:    Paare fallen monoton -> Extrap ~0.824 (~ref 0.825). SAUBER.
  - honeycomb: Paare fallen (NUR mit L=48!) -> Extrap ~0.574 (~ref 0.573/0.576).
               SAUBER - aber fragil: ohne L=48 (PHY035, L<=32) kippt es (~0.605).
  - triangular: Paare STEIGEN (von unten) bei sehr kleinem L (<=19) -> Extrap
               ~1.457, OBERHALB ref 1.418. UEBERSCHIESSEND/unzuverlaessig.
  => Die Klein-L-Paar-Extrapolation ist gitter-abhaengig und empfindlich auf das
  groesste enthaltene L; sie ist KEIN universelles Wundermittel. Das deckt sich
  mit HKS (L=48..192 noetig) und mit der bestehenden Repo-Ehrlichkeit.

EVIDENZ: deterministischer Gate-Report nach results/ (keine RNG, rein analytisch).
"""
from __future__ import annotations

import importlib.util
import math
import re
import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parent
_RESULTS = _SRC.parent / "results"


def _load(name: str, filename: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _SRC / filename)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


phy034 = _load("phy034_hks_extrapolation",
               "260611 PHY034 honeycomb hks extrapolation v01.py")
ols_intercept = phy034.ols_intercept
char_size = phy034.char_size
log_variable = phy034.log_variable

# Gitter -> (Report-Datei, Referenz(en), erwartete Paare fuer Drift-Guard).
LATTICES = {
    "honeycomb": {
        "report": "260607 PHY032 honeycomb wm-logfit bootstrap report.txt",
        "refs": (0.573, 0.576),
        "expected": {(12, 24): 0.6014, (12, 48): 0.5968, (24, 48): 0.5917},
    },
    "kagome": {
        "report": "260609 PHY033 kagome tbkt per-site report.txt",
        "refs": (0.825,),
        "expected": {(12, 24): 0.8437, (12, 36): 0.8410, (24, 36): 0.8377},
    },
    "triangular": {
        "report": "260604 PHY030 v02 wm-logfit report.txt",
        "refs": (1.418,),
        "expected": {(9, 13): 1.3663, (9, 19): 1.3852, (13, 19): 1.3888},
    },
}

_PAIR_RE = re.compile(r"T_BKT\(L=(\d+),(\d+)\)\s*=\s*([0-9.]+)")


def parse_pairs(report_name: str):
    """Parse die committed 'T_BKT(L=a,b) = X'-Paar-Zeilen eines Reports."""
    txt = (_RESULTS / report_name).read_text(encoding="utf-8")
    out = {}
    for m in _PAIR_RE.finditer(txt):
        l1, l2, val = int(m.group(1)), int(m.group(2)), float(m.group(3))
        out[(l1, l2)] = val
    return out


def extrapolate(pairs, power=2, size_mode="geom"):
    """L->inf-Intercept der Paar-Schaetzer in u=1/(ln Lc)^power (OLS)."""
    us, ts = [], []
    for (l1, l2), t in pairs.items():
        us.append(log_variable(char_size(l1, l2, size_mode), power))
        ts.append(t)
    if len({round(u, 12) for u in us}) < 2:
        return None
    intercept, _ = ols_intercept(us, ts)
    return intercept if math.isfinite(intercept) else None


def systematic_band(pairs):
    vals = []
    for mode in ("small", "geom", "large"):
        for p in (1, 2, 3):
            v = extrapolate(pairs, power=p, size_mode=mode)
            if v is not None:
                vals.append(v)
    return (min(vals), max(vals)) if vals else (None, None)


def bare_largest_pair(pairs):
    """Paar mit der groessten charakteristischen Groesse (konservative Schranke)."""
    key = max(pairs, key=lambda lp: char_size(lp[0], lp[1], "geom"))
    return key, pairs[key]


def _closest_ref_dev(value, refs):
    """Naechste Referenz + relative Abweichung in %."""
    ref = min(refs, key=lambda r: abs(value - r))
    return ref, (value - ref) / ref * 100.0


def analyse_lattice(name):
    cfg = LATTICES[name]
    pairs = parse_pairs(cfg["report"])
    refs = cfg["refs"]
    primary = extrapolate(pairs)               # p=2, geom (HKS-Standard)
    band = systematic_band(pairs)
    bkey, bval = bare_largest_pair(pairs)
    ref_e, dev_e = _closest_ref_dev(primary, refs)
    ref_b, dev_b = _closest_ref_dev(bval, refs)
    # "sauber" = Extrapolation naeher an ref als das blanke groesste Paar
    # UND innerhalb 1% einer Referenz.
    clean = abs(dev_e) < 1.0 and abs(dev_e) <= abs(dev_b) + 1e-9
    return {
        "pairs": pairs, "refs": refs, "primary": primary, "band": band,
        "bare_key": bkey, "bare_val": bval,
        "dev_extrap": dev_e, "ref_extrap": ref_e,
        "dev_bare": dev_b, "ref_bare": ref_b, "clean": clean,
    }


def run_phy037():
    res = {name: analyse_lattice(name) for name in LATTICES}

    # Drift-Guard: geparste Paare == committed Erwartung.
    drift_ok = True
    for name, cfg in LATTICES.items():
        got = res[name]["pairs"]
        for k, exp in cfg["expected"].items():
            if k not in got or abs(got[k] - exp) > 1e-4:
                drift_ok = False

    gates = {
        "PASS_DRIFT_GUARD_PAIRS_MATCH_REPORTS": drift_ok,
        "PASS_THREE_PAIRS_PER_LATTICE": all(
            len(res[n]["pairs"]) >= 3 for n in LATTICES),
        "PASS_ALL_EXTRAPOLATIONS_FINITE": all(
            res[n]["primary"] is not None for n in LATTICES),
        # Robuster Positiv-Befund: kagome konvergiert sauber zur Referenz.
        "PASS_KAGOME_CONVERGES_TO_REF": abs(res["kagome"]["dev_extrap"]) < 1.0,
    }
    overall = all(gates.values())
    return {"res": res, "gates": gates, "overall": overall}


def write_report(rep, path):
    res = rep["res"]
    L = []
    L.append("PHY037 v01 - HKS-L->inf-Extrapolation GITTER-UEBERGREIFEND")
    L.append("            (honeycomb / kagome / triangular)")
    L.append("Selbst-Audit / Coworkerz, 2026-06-11")
    L.append("=" * 72)
    L.append("")
    L.append("Methode: C-eliminierte Paar-Schaetzer (committed Reports) linear")
    L.append("  in u=1/(ln Lc)^2 (geom. Mittel) zum Intercept L->inf (HKS,")
    L.append("  arXiv:1302.2900). Rein analytisch, kein neuer MC-Lauf.")
    L.append("")
    hdr = (f"  {'Gitter':11} {'Extrap':>8} {'Abw':>8} {'bestes Paar':>12} "
           f"{'Abw':>8} {'Systematik-Band':>20} {'Verdikt':>10}")
    L.append(hdr)
    for name in ("honeycomb", "kagome", "triangular"):
        r = res[name]
        blo, bhi = r["band"]
        band = f"[{blo:.4f},{bhi:.4f}]" if blo is not None else "[n/a]"
        verdikt = "SAUBER" if r["clean"] else "UEBERSCHIESST"
        L.append(f"  {name:11} {r['primary']:8.4f} {r['dev_extrap']:+7.2f}% "
                 f"{r['bare_val']:12.4f} {r['dev_bare']:+7.2f}% "
                 f"{band:>20} {verdikt:>10}")
    L.append("")
    for name in ("honeycomb", "kagome", "triangular"):
        r = res[name]
        L.append(f"--- {name} (refs {', '.join(str(x) for x in r['refs'])}) ---")
        for (l1, l2), t in sorted(
                r["pairs"].items(),
                key=lambda kv: char_size(kv[0][0], kv[0][1], "geom")):
            lc = char_size(l1, l2, "geom")
            L.append(f"  Paar({l1},{l2}) = {t:.4f}  [Lc(geom)={lc:.2f}]")
        L.append(f"  -> L->inf (p=2,geom) = {r['primary']:.4f}  "
                 f"({r['dev_extrap']:+.2f}% vs {r['ref_extrap']})")
        L.append("")
    L.append("--- BEFUND (ehrliche Methoden-Bilanz) ---")
    L.append("  kagome:    Paare fallen monoton -> Extrap ~ref. SAUBER.")
    L.append("  honeycomb: Paare fallen (NUR mit L=48) -> Extrap ~ref. SAUBER,")
    L.append("             aber fragil (ohne L=48 kippt es, s. PHY035).")
    L.append("  triangular: Paare STEIGEN bei L<=19 -> Extrap UEBERSCHIESST ueber")
    L.append("             1.418; das blanke groesste Paar ist hier naeher an ref.")
    L.append("  => Klein-L-Paar-Extrapolation ist gitter-abhaengig + empfindlich")
    L.append("  auf das groesste enthaltene L - KEIN universelles Wundermittel")
    L.append("  (deckt sich mit HKS: L=48..192 noetig). Das blanke groesste Paar")
    L.append("  bleibt eine konkurrenzfaehige konservative Schranke.")
    L.append("")
    L.append("--- PASS-Gates (Analyse-Integritaet + robuster kagome-Positiv) ---")
    for k, v in rep["gates"].items():
        L.append(f"  [{'PASS' if v else 'FAIL'}] {k}")
    L.append(f"  OVERALL: {'PASS' if rep['overall'] else 'FAIL'}")
    L.append("")
    L.append("--- Caveats ---")
    L.append("  - triangular/honeycomb-Extrapolationen aus <=3 kleinen L sind")
    L.append("    sub-leading-log-limitiert (ehrlich; vgl. PHY030/032-Caveats).")
    L.append("  - 'SAUBER/UEBERSCHIESST' ist eine qualitative Methoden-Aussage,")
    L.append("    kein neuer Praezisions-Punktschaetzer fuer T_BKT.")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> int:
    rep = run_phy037()
    out = _RESULTS / "260611 PHY037 hks multilattice extrapolation report.txt"
    write_report(rep, out)
    print(f"PHY037 Report -> {out}")
    for name in ("honeycomb", "kagome", "triangular"):
        r = rep["res"][name]
        print(f"  {name:11} extrap={r['primary']:.4f} "
              f"({r['dev_extrap']:+.2f}%)  "
              f"{'SAUBER' if r['clean'] else 'UEBERSCHIESST'}")
    print(f"  OVERALL: {'PASS' if rep['overall'] else 'FAIL'}")
    return 0 if rep["overall"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
