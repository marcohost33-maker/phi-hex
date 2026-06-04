# Phi-Hex

> **Status:** Forschungs-Repo (privat) | Aktive Forschungslinie (letzte Haertung 2026-06-02). Kern-Engine kompiliert; Selftest-Evidenz vom 2026-06-02 liegt in results/.
> **Lineage/Provenance:** siehe `SOURCES.md` (SHA-256 je Quelldatei) | **Lizenz:** Apache-2.0

XY-Modell / BKT-Physik auf Dreiecks- und Honeycomb-Gittern: Helicity-Modulus, Nelson-Kosterlitz-Sprung, Wolff-Cluster, Finite-Size-Scaling.

## Kern

- **Engine:** `src/260602 PHI HEX core v2 2 hardened.py` (Kern, 2026-06-02) + PHY024-029-Experiment-Serie
- **Evidenz:** Selftest-Report + PHY025-029 Reports in `results/`; methodischer Audit (4 Responses) in `spec/`. Referenzlinie T_BKT(triangular) ~= 1.4.

## Struktur

```
src/        Engines + Experiment-Code
spec/       Specs, Theorie, Audits (gehaertete Versionen)
results/    Selftest-/Gate-Logs, Reports, Negativ-Results
archive/    Vorgaenger-Versionen (Lineage)
docs/ADR/   Architektur-Entscheide
SOURCES.md  Provenance: Quelle + SHA-256 + mtime je Datei
```

## Reproduzieren

Die Selftests laufen direkt in der Engine (`python <engine>.py --selftest` bzw. eingebaute Runner;
Details im Engine-Header). CI prueft Syntax (`compileall`) + Lint (non-blocking) — volle
physikalische Gate-Runs brauchen numpy/scipy und laufen lokal.

## Verwandte Repos

`coworker-dde` (kanonische v3-Engine-Quelle, 3DDE-Kanon) | `HQST` (p1_tasks-Floquet-Linie) |
`phi-hex` `hex-hqst` `u2-dualspec` `u6-bifurcation` `hexa-ntk` (diese Serie, 2026-06-04)

---
*Coworker Research / Coworkerz | Repo-Anlage 2026-06-04 nach arbeitsschablone_forschungs-repo-anlage*
