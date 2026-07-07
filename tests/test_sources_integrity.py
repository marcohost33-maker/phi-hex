"""SHA-First-Integritaets-Gate fuer SOURCES.md (Audit 2026-07-07, Issue #21).

Vertrag: SOURCES.md ist die Integritaets-Wahrheit (AGENTS.md §Working
agreements, "SHA-First"). Dieses Gate erzwingt das maschinell:

1. Jede committete Datei unter src/ tests/ spec/ results/ archive/ muss
   mindestens eine SOURCES.md-Zeile haben, deren 16-Zeichen-SHA-256-Praefix
   die aktuellen Datei-Bytes matcht. Fuer Alt-Eintraege aus dem Windows-
   Staging gelten CRLF/LF(+BOM)-Varianten als Match.
2. Der PHY042-Gate-Report ist zusaetzlich mit voller SHA-256 gepinnt
   (Issue #21: Provenance-Reparatur, keine Ergebnis-Aenderung).

Wird eine gedeckte Datei geaendert, schlaegt (1) fehl, bis eine neue Zeile
angehaengt ist (append-only; Alt-Zeilen bleiben als Lineage stehen).
"""
from __future__ import annotations

import hashlib
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ROOT / "SOURCES.md"
SCOPE = ("src", "tests", "spec", "results", "archive")

PHY042_REPORT = (
    ROOT / "results" / "260707 PHY042 honeycomb wl-fss L24-32-48 gate report.json"
)
PHY042_SHA256_FULL = (
    "19a9ce3c799401dbb55e519c09b390bee9dfbe23792f99b2e7c650c3eeefa3cf"
)


def _row_hash_prefixes() -> set[str]:
    text = SOURCES.read_text(encoding="utf-8")
    return set(re.findall(r"\| ([0-9A-F]{16}) \|", text))


def _byte_variants(raw: bytes) -> list[bytes]:
    lf = raw.replace(b"\r\n", b"\n")
    crlf = lf.replace(b"\n", b"\r\n")
    bom = b"\xef\xbb\xbf"
    return [raw, lf, crlf, bom + lf, bom + crlf]


def _scoped_files() -> list[Path]:
    out = subprocess.run(
        ["git", "ls-files", "-z", "--", *SCOPE],
        cwd=ROOT,
        capture_output=True,
        check=True,
    ).stdout
    files = [ROOT / p.decode("utf-8") for p in out.split(b"\0") if p]
    assert files, "git ls-files lieferte keine Dateien im Scope"
    return files


def _is_covered(path: Path, prefixes: set[str]) -> bool:
    raw = path.read_bytes()
    return any(
        hashlib.sha256(v).hexdigest().upper()[:16] in prefixes
        for v in _byte_variants(raw)
    )


def test_every_scoped_file_has_matching_sources_row() -> None:
    prefixes = _row_hash_prefixes()
    uncovered = [
        str(p.relative_to(ROOT))
        for p in _scoped_files()
        if not _is_covered(p, prefixes)
    ]
    assert not uncovered, (
        "Dateien ohne byte-genau passende SOURCES.md-Zeile (neue Zeile "
        "anhaengen, append-only): " + ", ".join(sorted(uncovered))
    )


def test_phy042_gate_report_pinned_full_sha256() -> None:
    actual = hashlib.sha256(PHY042_REPORT.read_bytes()).hexdigest()
    assert actual == PHY042_SHA256_FULL, (
        "PHY042-Gate-Report weicht vom in SOURCES.md gepinnten Stand ab: "
        f"{actual} != {PHY042_SHA256_FULL}"
    )
    text = SOURCES.read_text(encoding="utf-8")
    assert PHY042_SHA256_FULL in text, (
        "volle PHY042-SHA-256 fehlt in SOURCES.md"
    )
    assert actual.upper()[:16] in _row_hash_prefixes(), (
        "16-Zeichen-Praefix des PHY042-Reports fehlt als SOURCES.md-Zeile"
    )
