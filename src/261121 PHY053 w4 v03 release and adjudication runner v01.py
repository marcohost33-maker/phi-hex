"""PHY053 - W4 v03: Freigabe, Produktion und Adjudikation (reiner Runner).

Coworker Research / Coworkerz, 21. November 2026
de-CH Konventionen, ASCII-Quotes, ss statt Eszett, keine Personennamen.

Dieses Modul enthaelt KEINE eigene Physik, keine Schwellen und keine
Schaetzer. Es verbindet nur die vorregistrierten Bausteine in der
vorgeschriebenen Reihenfolge:

  1. Freigabe (Spec PHY052 Abschnitt 5): committeter PHY052-Pilot-Report
     `MIXING_PILOT_PASS` UND nativer PHY049-Preflight (5 Gates literal
     true) -> `ELIGIBLE_TO_START_PHY049_MEASUREMENT`, sonst Abbruch.
  2. Produktion: `PHY049.produce` mit den vorregistrierten Defaults
     (5x29x12, 1000+4000 Sweeps, 4 Prozesse, 24 h) und Checkpoint/Resume.
  3. Adjudikation: ausschliesslich `PHY050.assess_production` (v03a
     Abschnitt 9: "Single adjudicator"), 1000 Bootstrap-Replikate.

Warum ein eigenes Modul: PHY049, PHY050 und PHY045 stehen im
Provenienz-Fingerprint des Produkts (v03a Abschnitt 12). Jede Aenderung
dort - auch ein CLI-Zusatz - machte das Produkt unadjudizierbar. PHY053
steht bewusst NICHT im Fingerprint.

Claim-Decke: das Adjudikations-Label ist interne Diskriminationssemantik
(v03 Abschnitt 6/8), kein Bestwert; externe Aussage erst nach
Cross-Family-Review.
"""
from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse  # noqa: E402
import hashlib  # noqa: E402
import importlib.util  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

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


_p52 = _load("phy052_mixing_pilot",
             "261121 PHY052 honeycomb correlation-ratio mixing pilot v01.py")
_p49 = _p52._p49
_p50 = _p49._p50

PILOT_REPORT = _ROOT / "results" / f"{_p52.REPORT_STEM}.json"
PRODUCT_PATH = Path(_p49.W4V3_RESULT_PATH)
ADJ_STEM = "261121 PHY053 w4 v03 production adjudication report"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ============================================================================
# 1. Freigabe
# ============================================================================

def release_decision(pilot_report: dict | None, preflight: dict) -> dict:
    """Freigabe nur mit committetem, vorregistriertem PASS-Pilot und
    nativem PHY049-Preflight. Der Pilot-Entscheid wird aus den Zellen neu
    berechnet (nicht dem gespeicherten Label geglaubt)."""
    if not isinstance(pilot_report, dict) or "pilot" not in pilot_report:
        return {"label": "DO_NOT_START_PHY049_MEASUREMENT",
                "reason": "no PHY052 pilot report"}
    recomputed = _p52.decision(pilot_report["pilot"])
    if recomputed != pilot_report.get("decision"):
        return {"label": "DO_NOT_START_PHY049_MEASUREMENT",
                "reason": "stored pilot decision disagrees with its cells"}
    label = _p52.production_release(recomputed, preflight)
    return {"label": label, "pilot_decision": recomputed["label"],
            "preflight_gates": preflight.get("gates"),
            "reason": None if label.startswith("ELIGIBLE") else
            "pilot not PASS or PHY049 preflight not green"}


def load_pilot_report(path: Path = PILOT_REPORT) -> dict | None:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


# ============================================================================
# 2. Produktion (vorregistrierte Defaults, nichts ueberschreibbar)
# ============================================================================

def run_production(*, resume: bool = False,
                   product_path: Path = PRODUCT_PATH) -> dict:
    pilot = load_pilot_report()
    pre = _p49.preflight()
    rel = release_decision(pilot, pre)
    if rel["label"] != "ELIGIBLE_TO_START_PHY049_MEASUREMENT":
        raise RuntimeError(f"production not released: {rel}")
    t0 = time.perf_counter()
    prod = _p49.produce(checkpoint_path=product_path, resume=resume)
    return {"release": rel, "wall_s_this_call": time.perf_counter() - t0,
            "checkpoint_status": prod.get("checkpoint_status"),
            "complete": prod.get("complete"), "rows": len(prod["rows"])}


# ============================================================================
# 3. Adjudikation
# ============================================================================

def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(x) for x in o]
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    if hasattr(o, "item") and not isinstance(o, (str, bytes)):
        return _clean(o.item())
    return o


def adjudicate(prod: dict, *, product_sha256: str | None = None,
               n_boot: int = _p50.W4V3_N_BOOT) -> dict:
    """PHY050 entscheidet; PHY053 legt nur Kontext daneben (Kurven,
    Splay-Diagnostik), der in kein Label eingeht."""
    t0 = time.perf_counter()
    assessment = _p50.assess_production(prod, n_boot=n_boot)
    agg = _p49.aggregate(prod)
    splay = {}
    ladder = agg["ladder"]
    for i, L1 in enumerate(ladder):
        for L2 in ladder[i + 1:]:
            c1, c2 = agg["curves"][str(L1)], agg["curves"][str(L2)]
            splay[f"{L1}_{L2}"] = _p49.persistent_splay(
                agg["t_grid"], c1["R"], c1["SE"], c2["R"], c2["SE"])
    return {
        "module": "PHY053_w4_v03_release_and_adjudication_v01",
        "attribution": "Coworker Research / Coworkerz",
        "adjudicator": "PHY050.assess_production (single adjudicator, "
                       "v03a section 9)",
        "product": {"path": str(PRODUCT_PATH.relative_to(_ROOT)),
                    "sha256": product_sha256,
                    "checkpoint_status": prod.get("checkpoint_status"),
                    "rows": len(prod.get("rows", [])),
                    "wall_s": prod.get("wall_s"),
                    "runtime_provenance": prod.get("runtime_provenance")},
        "pilot_report_sha256": (_sha256(PILOT_REPORT)
                                if PILOT_REPORT.exists() else None),
        "assessment": assessment,
        "context_not_decisive": {
            "curves": agg["curves"], "t_grid": agg["t_grid"],
            "persistent_splay_T_by_pair": splay},
        "n_boot": n_boot,
        "adjudication_wall_s": time.perf_counter() - t0,
        "claim_ceiling": assessment.get("claim_ceiling"),
    }


def write_text_report(rep: dict, path: Path) -> None:
    a = rep["assessment"]

    def f(x, fmt=".4f"):
        return "null" if x is None else format(x, fmt)
    pr = a.get("primary") or {}
    boot = a.get("bootstrap") or {}
    lines = ["PHY053 - W4 v03 Adjudikation (PHY050.assess_production)",
             "Coworker Research / Coworkerz, 2026-11-21", "=" * 70,
             f"Produkt: {rep['product']['path']} (sha256 "
             f"{rep['product']['sha256']}), Status "
             f"{rep['product']['checkpoint_status']}, "
             f"{rep['product']['rows']} Zeilen, Wall "
             f"{f(rep['product']['wall_s'], '.0f')} s",
             f"Pilot-Report sha256 {rep['pilot_report_sha256']}", "",
             f"Gates: {a.get('gates')}",
             f"Primaer (Collapse): T_BKT = {f(pr.get('tbkt'))}, c = "
             f"{f(pr.get('c'))}, Score {f(pr.get('score'), '.3f')}, "
             f"Profilbreite {f(pr.get('profile_width'))}",
             f"Bootstrap: sigma_boot = {f(a.get('sigma_boot'))} "
             f"({boot.get('bootstrap_valid')}/{boot.get('bootstrap_requested')}"
             f" gueltig), 95%-Quantile [{f(boot.get('q025'))}, "
             f"{f(boot.get('q975'))}]",
             f"sigma_model = {f(a.get('sigma_model'))}, sigma_tot = "
             f"{f(a.get('sigma_tot'))}", ""]
    rob = (a.get("robustness") or {}).get("variants", {})
    for name, v in rob.items():
        lines.append(f"  Variante {name}: T = {f(v.get('tbkt'))}, "
                     f"quotierbar {v.get('quotable')}, |Delta| = "
                     f"{f(v.get('delta_primary'))}")
    lines += ["", f"ENTSCHEID: {a.get('decision')} "
              f"(Physik-Interpretation freigegeben: "
              f"{a.get('physics_interpretation_enabled')})",
              f"Claim-Decke: {a.get('claim_ceiling')}", "",
              "Kontext (nicht entscheidend): persistenter 2-sigma-Splay je "
              "Paar (L1,L2):"]
    for k, v in rep["context_not_decisive"][
            "persistent_splay_T_by_pair"].items():
        lines.append(f"  ({k.replace('_', ',')}): {f(v)}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("release", help="nur Freigabe pruefen")
    p = sub.add_parser("produce", help="Freigabe pruefen, dann PHY049")
    p.add_argument("--resume", action="store_true")
    sub.add_parser("adjudicate", help="PHY050 auf das Produkt anwenden")
    ns = ap.parse_args(argv)
    if ns.cmd == "release":
        rel = release_decision(load_pilot_report(), _p49.preflight())
        print(json.dumps(_clean(rel), indent=1))
        return 0 if rel["label"].startswith("ELIGIBLE") else 3
    if ns.cmd == "produce":
        print(json.dumps(_clean(run_production(resume=ns.resume)), indent=1))
        return 0
    prod = json.loads(PRODUCT_PATH.read_text(encoding="utf-8"))
    rep = _clean(adjudicate(prod, product_sha256=_sha256(PRODUCT_PATH)))
    res = _ROOT / "results"
    (res / f"{ADJ_STEM}.json").write_text(
        json.dumps(rep, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    write_text_report(rep, res / f"{ADJ_STEM}.txt")
    print(f"PHY053: {rep['assessment'].get('decision')}; results/"
          f"{ADJ_STEM}.{{json,txt}}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
