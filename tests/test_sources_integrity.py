"""SHA-First-Integritaets-Gate fuer SOURCES.md (Audit 2026-07-07, Issue #21;
gehaertet 2026-07-10, Code-Audit).

Vertrag: SOURCES.md ist die Integritaets-Wahrheit (AGENTS.md §Working
agreements, "SHA-First"). Dieses Gate erzwingt das maschinell:

1. Jede committete Datei unter src/ tests/ spec/ results/ archive/ docs/
   muss eine SOURCES.md-Zeile haben, deren 16-Zeichen-SHA-256-Praefix die
   aktuellen Datei-Bytes matcht UND deren Ziel-Spalte auf das Verzeichnis
   der Datei zeigt (Haertung 2026-07-10: vorher war das Matching pfad-
   ungebunden - ein `git mv` oder eine byte-identische Kopie unter neuem
   Namen passierte das Gate ohne neue Zeile).
2. CRLF/LF(+BOM)-Varianten gelten NUR fuer Windows-Staging-Zeilen (Quelle
   mit Laufwerkspfad); repo-native Zeilen muessen byte-genau matchen
   (Haertung 2026-07-10: vorher wurden Zeilenenden fuer ALLE Zeilen
   gewaschen).
3. Kein getracktes .py darf ausserhalb der Scope-Verzeichnisse liegen
   (Haertung 2026-07-10: ein Root-Level-Helper haette sonst SOURCES-, Lint-
   und Compile-Gate gleichzeitig umgangen).
4. Der PHY042-Gate-Report ist zusaetzlich mit voller SHA-256 gepinnt
   (Issue #21).

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
SCOPE = ("src", "tests", "spec", "results", "archive", "docs")

PHY042_REPORT = (
    ROOT / "results" / "260707 PHY042 honeycomb wl-fss L24-32-48 gate report.json"
)
PHY042_SHA256_FULL = (
    "19a9ce3c799401dbb55e519c09b390bee9dfbe23792f99b2e7c650c3eeefa3cf"
)

_ROW_RE = re.compile(
    r"^\| ([0-9A-F]{16}) \| [^|]* \| ([^|]+) \| (.*) \|\s*$", re.MULTILINE)


def _rows() -> list[tuple[str, str, bool]]:
    """(prefix, ziel_normalisiert, ist_windows_staging) je SOURCES-Zeile."""
    text = SOURCES.read_text(encoding="utf-8")
    rows = []
    for m in _ROW_RE.finditer(text):
        prefix = m.group(1)
        ziel = m.group(2).strip().replace("\\", "/")
        quelle = m.group(3)
        is_windows = bool(re.search(r"\b[A-Z]:\\", quelle))
        rows.append((prefix, ziel, is_windows))
    assert rows, "keine Hash-Zeilen in SOURCES.md gefunden"
    return rows


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


def _is_covered(path: Path, rows: list[tuple[str, str, bool]]) -> bool:
    raw = path.read_bytes()
    parent = str(path.relative_to(ROOT).parent).replace("\\", "/")
    sha_raw = hashlib.sha256(raw).hexdigest().upper()[:16]
    sha_variants = None  # lazy: nur fuer Windows-Zeilen berechnen
    for prefix, ziel, is_windows in rows:
        if ziel != parent:
            continue
        if prefix == sha_raw:
            return True
        if is_windows:
            if sha_variants is None:
                sha_variants = {
                    hashlib.sha256(v).hexdigest().upper()[:16]
                    for v in _byte_variants(raw)
                }
            if prefix in sha_variants:
                return True
    return False


def test_every_scoped_file_has_matching_sources_row() -> None:
    rows = _rows()
    uncovered = [
        str(p.relative_to(ROOT))
        for p in _scoped_files()
        if not _is_covered(p, rows)
    ]
    assert not uncovered, (
        "Dateien ohne byte-genau passende, pfadgebundene SOURCES.md-Zeile "
        "(neue Zeile anhaengen, append-only; Ziel-Spalte = Verzeichnis): "
        + ", ".join(sorted(uncovered))
    )


def test_no_tracked_python_outside_scope() -> None:
    """Scope-Escape-Gate: getrackte .py nur innerhalb der Scope-Dirs."""
    out = subprocess.run(
        ["git", "ls-files", "-z", "--", "*.py"],
        cwd=ROOT,
        capture_output=True,
        check=True,
    ).stdout
    offenders = [
        p.decode("utf-8")
        for p in out.split(b"\0")
        if p and not p.decode("utf-8").startswith(tuple(d + "/" for d in SCOPE))
    ]
    assert not offenders, (
        "getrackte .py ausserhalb der Scope-Verzeichnisse (entziehen sich "
        "SOURCES-/Lint-/Compile-Gates): " + ", ".join(sorted(offenders))
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
    assert any(actual.upper()[:16] == pfx for pfx, _, _ in _rows()), (
        "16-Zeichen-Praefix des PHY042-Reports fehlt als SOURCES.md-Zeile"
    )
