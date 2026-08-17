"""SHA-First-Integritaets-Gate fuer SOURCES.md (Audit 2026-07-07, Issue #21;
gehaertet 2026-07-10, Code-Audit; plattform-portabel 2026-08-17).

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
5. Der Arbeitsbaum darf inhaltlich nicht vom Index abweichen (siehe unten).

Wird eine gedeckte Datei geaendert, schlaegt (1) bzw. (5) fehl, bis eine neue
Zeile angehaengt ist (append-only; Alt-Zeilen bleiben als Lineage stehen).

GEPRUEFTES OBJEKT (Haertung 2026-08-17)
---------------------------------------
(1), (4) und (5) hashen die Bytes aus dem **Git-Index**, nicht die des
Arbeitsbaums. Grund: SOURCES.md pinnt die *committeten* Bytes; der
Arbeitsbaum ist nur deren plattformabhaengige Darstellung. Unter Windows mit
`core.autocrlf=true` (und ohne `.gitattributes`) checkt git Textdateien mit
CRLF aus, waehrend der Index LF haelt — dieselbe Datei hat dann zwei
verschiedene SHA-256. Gemessen am 2026-08-17 auf b64a9c5: 88 von 88 gescopten
Dateien wichen Arbeitsbaum-vs-Index ab, **alle 88 ausschliesslich in den
Zeilenenden**, echte Byte-Abweichung 0. Das Gate meldete damit lokal auf
Windows 68 falsche Positive, waehrend Linux-CI gruen war — es mass das
falsche Objekt, nicht einen Integritaetsbruch.

Der Index-Bezug schwaecht die Regel-2-Haertung NICHT auf: Zeilenenden werden
weiterhin nirgends "gewaschen", der Vergleich bleibt byte-genau — er laeuft
nur gegen das Objekt, das SOURCES.md tatsaechlich bindet. Damit ein
uncommitteter Arbeitsbaum-Edit trotzdem sofort auffaellt (dieses Signal
haette der Index-Bezug allein gekostet), prueft (5) zusaetzlich Arbeitsbaum
gegen Index modulo Zeilenenden.
"""
from __future__ import annotations

import hashlib
import re
import subprocess
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ROOT / "SOURCES.md"
SCOPE = ("src", "tests", "spec", "results", "archive", "docs")

PHY042_REPORT = "results/260707 PHY042 honeycomb wl-fss L24-32-48 gate report.json"
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


def _lf(raw: bytes) -> bytes:
    """Zeilenenden-normalisierte Sicht — NUR fuer den Arbeitsbaum-Drift-Check."""
    return raw.replace(b"\r\n", b"\n")


def _index_entries() -> dict[str, str]:
    """rel-Pfad (posix, wie git ihn fuehrt) -> Blob-OID aus dem Git-Index."""
    out = subprocess.run(
        ["git", "ls-files", "-s", "-z", "--", *SCOPE],
        cwd=ROOT,
        capture_output=True,
        check=True,
    ).stdout
    entries: dict[str, str] = {}
    for rec in out.split(b"\0"):
        if not rec:
            continue
        meta, _, rel = rec.decode("utf-8").partition("\t")
        mode, oid, stage = meta.split()
        assert stage == "0", (
            f"Index-Eintrag mit Merge-Stage {stage} (ungeloester Konflikt): {rel}"
        )
        assert mode != "160000", (
            f"Submodul im Scope — hat keinen Blob, Gate nicht anwendbar: {rel}"
        )
        entries[rel] = oid
    assert entries, "git ls-files lieferte keine Dateien im Scope"
    return entries


def _blob_bytes(oids: list[str]) -> dict[str, bytes]:
    """OID -> Blob-Inhalt, in EINEM `git cat-file --batch`-Aufruf."""
    uniq = sorted(set(oids))
    if not uniq:
        return {}
    proc = subprocess.run(
        ["git", "cat-file", "--batch"],
        cwd=ROOT,
        input=("\n".join(uniq) + "\n").encode("ascii"),
        capture_output=True,
        check=True,
    )
    out, pos, blobs = proc.stdout, 0, {}
    for _ in uniq:
        nl = out.index(b"\n", pos)
        header = out[pos:nl].decode("ascii").split()
        assert len(header) == 3 and header[1] == "blob", (
            f"unerwarteter cat-file-Header: {header!r}"
        )
        start = nl + 1
        size = int(header[2])
        blobs[header[0]] = out[start:start + size]
        pos = start + size + 1  # cat-file haengt ein \n an den Inhalt
    assert len(blobs) == len(uniq), (
        f"cat-file lieferte {len(blobs)} statt {len(uniq)} Blobs"
    )
    return blobs


def _is_covered(rel: str, raw: bytes, rows: list[tuple[str, str, bool]]) -> bool:
    parent = str(PurePosixPath(rel).parent)
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
    entries = _index_entries()
    blobs = _blob_bytes(list(entries.values()))
    uncovered = [
        rel
        for rel, oid in entries.items()
        if not _is_covered(rel, blobs[oid], rows)
    ]
    assert not uncovered, (
        "Dateien ohne byte-genau passende, pfadgebundene SOURCES.md-Zeile "
        "(neue Zeile anhaengen, append-only; Ziel-Spalte = Verzeichnis): "
        + ", ".join(sorted(uncovered))
    )


def test_worktree_matches_index_modulo_line_endings() -> None:
    """Uncommitteter Arbeitsbaum-Edit faellt sofort auf (Regel 5).

    Ohne diesen Test wuerde der Index-Bezug von Regel (1)/(4) das lokale
    Sofort-Signal kosten: eine geaenderte, aber ungestagte Datei bliebe
    unsichtbar bis zum `git add`. Zeilenenden sind bewusst ausgenommen —
    genau sie sind die plattformabhaengige Darstellung, nicht der Inhalt.
    """
    entries = _index_entries()
    blobs = _blob_bytes(list(entries.values()))
    drift = []
    for rel, oid in entries.items():
        path = ROOT / rel
        if not path.is_file():
            drift.append(f"{rel} (im Index, fehlt im Arbeitsbaum)")
        elif _lf(path.read_bytes()) != _lf(blobs[oid]):
            drift.append(rel)
    assert not drift, (
        "Arbeitsbaum weicht inhaltlich vom Index ab (stagen und eine neue "
        "SOURCES.md-Zeile anhaengen): " + ", ".join(sorted(drift))
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
    entries = _index_entries()
    assert PHY042_REPORT in entries, (
        f"PHY042-Gate-Report nicht im Index-Scope: {PHY042_REPORT}"
    )
    raw = _blob_bytes([entries[PHY042_REPORT]])[entries[PHY042_REPORT]]
    actual = hashlib.sha256(raw).hexdigest()
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
