# ADR-001: Warum dieses Repo existiert

**Status:** accepted | **Datum:** 2026-06-04 | **Entscheider:** Marco (User1) + Vero

## Kontext
Phi-Hex-Material lag verstreut (Google Drive 0-Folder, Desktop Universe 2.0, coworker-dde).
Marco-Auftrag 2026-06-04: pro Forschungsthema ein sauberes privates GitHub-Repo; Daten zuvor
akribisch suchen, ordnen, SHA-dedupen (64 exakte Dubletten quarantaeniert, 24 dank SHA-Check
behalten); Quell-Ordner (Projekte/Projekte 2/Universe 2.0) bleiben unveraendert.

## Entscheidung
Eigenes Repo `phi-hex` nach `arbeitsschablone_forschungs-repo-anlage` (Phase 0-6): SHA-Provenance
in SOURCES.md, ehrliche Status-Zeilen, CI = Syntax-Gate, kanonische v3-Engine-Quelle bleibt
coworker-dde (byte-identisch verifiziert 2026-06-04, Equalita-Audit 4/4 PASS).

## Konsequenzen
+ Ein Ort pro Thema, nachvollziehbare Lineage, CI-Schutz.
- Engine-Doppelhaltung mit coworker-dde -> Sync-Disziplin noetig (AGENTS.md §Project).

## Nachtrag 2026-07-10 (Code-Audit; Historie oben unveraendert)

Zwei Aussagen dieses ADR sind ueberholt und werden hier datiert korrigiert
(nicht umgeschrieben - Lineage-Konvention):

1. **dde-Doppelhaltung:** AGENTS.md §Project sagt inzwischen ausdruecklich
   "KEINE dde-Doppelhaltung" fuer phi-hex (Quelle ist Drive `0Phi Hex` +
   `Universe 2.0\Hqstphi Hex`). Die im Konsequenzen-Punkt genannte
   Sync-Disziplin mit coworker-dde entfaellt damit fuer dieses Repo.
2. **"CI = Syntax-Gate":** Die CI ist seit den PRs #6/#15/#22 ein
   mehrstufiges Gate (ruff-Lint blockierend, compileall, pytest-Matrix
   3.10-3.13 mit Korrektheits-/Integritaets-Gates inkl. SHA-First-Gate
   `tests/test_sources_integrity.py`).
