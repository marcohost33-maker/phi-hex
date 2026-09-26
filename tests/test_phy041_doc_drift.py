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
    """Kein Einzelanker-Vergleich ("vs 0.576") in den Vertrags-Dokumenten.

    Nachtrag 2026-09-26 (Issue #45 §3, spec/260926 provenance addendum §3):
    bis dahin verbot dieses Gate zusaetzlich jede Naehe von "0.576(3)" und
    "2406.12076" - auf der Annahme, 0.576(3) sei dieser Quelle faelschlich
    zugeschrieben. Das Abstract (v4 / Physica Scripta) berichtet 0.576(3)
    aber direkt als WL-T-Wert. Das Attributions-Verbot ist deshalb durch
    die Positiv-Pruefung unten ersetzt (0.576(3) nur als direkter T-Kanal
    IM Band); das Einzelanker-Verbot bleibt unveraendert.
    """
    offenders: list[str] = []
    stale_patterns = [
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


def test_0576_is_attributed_only_as_direct_t_channel_in_band() -> None:
    """Nachtrag 2026-09-26: README fuehrt 0.576(3) als eigene Bandzeile
    (direkter T-Kanal, WL) - nie als Umrechnung eines beta-Kanals."""
    readme = _read(ROOT / "README.md")
    row = re.search(r"^\| arXiv:2406\.12076[^|]*\| T_BKT = 0\.576\(3\) \|"
                    r"[^|]*\|[^|]*direkt T[^|]*\|$", readme, re.MULTILINE)
    assert row, "0.576(3)-Bandzeile (direkter T-Kanal) fehlt im README"
    assert not re.search(r"beta_BKT = [0-9.()]+ \| 0\.576", readme)
