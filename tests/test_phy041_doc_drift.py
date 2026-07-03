"""Doku-Drift-Gate fuer PHY041-Honeycomb-Referenzen.

Der Review-Vertrag seit 2026-07-03: Honeycomb-Referenzen werden als Band
mit Quelle, Observable und beta/T-Konvention gefuehrt. Eine Rueckkehr zur
alten Einzelanker-Sprache ist nicht mergefaehig.
"""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

DOC_PATHS = [
    ROOT / "README.md",
    ROOT / "spec" / "260702 PHI HEX phy041 honeycomb entropic tbkt method v01.md",
    ROOT / "spec" / "260703 PHI HEX honeycomb reference conventions audit v01.md",
    ROOT / "results" / "260702 PHY041 honeycomb wang-landau entropic helicity report.txt",
]


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_no_stale_single_anchor_attribution() -> None:
    offenders: list[str] = []
    stale_patterns = [
        re.compile(r"0\.576\(3\).{0,40}2406\.12076"),
        re.compile(r"2406\.12076.{0,40}0\.576\(3\)"),
        re.compile(r"vs\s+0\.576(?:\(3\))?"),
    ]
    for path in DOC_PATHS:
        text = _read(path)
        for pat in stale_patterns:
            if pat.search(text):
                offenders.append(str(path.relative_to(ROOT)))
                break
    assert not offenders, "stale honeycomb reference wording in: " + ", ".join(offenders)


def test_reference_band_contract_is_present() -> None:
    combined = "\n".join(_read(path).lower() for path in DOC_PATHS)
    for token in ["reference", "beta", "t = 1 / beta", "0.5928", "0.6116", "0.5800"]:
        assert token in combined
