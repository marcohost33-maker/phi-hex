# Phi-Hex

> **Status:** Forschungs-Repo (privat) | Aktive Forschungslinie (letzte Haertung 2026-06-04). Kern-Engine kompiliert; Selftest- + T_BKT-Mess-Evidenz in results/. Testsuite 14/14 PASS (lokal, 2026-06-04; +7 Weber-Minnhagen-Log-Fit).
> **Lineage/Provenance:** siehe `SOURCES.md` (SHA-256 je Quelldatei) | **Lizenz:** Apache-2.0

XY-Modell / BKT-Physik auf Dreiecks- und Honeycomb-Gittern: Helicity-Modulus, Nelson-Kosterlitz-Sprung, Wolff-Cluster, Finite-Size-Scaling.

## Konvention: Helicity-Modulus per-Site (Stand 2026-06-04)

> **Konventionswechsel 2026-06-04 (Marco-Entscheid, Wissenschafts-Audit S2/S3).**
> Der Helicity-Modulus wird jetzt **per-Site** normiert (Division durch die
> Site-Zahl N), Literatur-Standard (arXiv:2406.12076). Spinwellen-Grenzwert
> Dreiecksgitter: **Υ(0) = 1.5 J** (exakt, getestet, tol 0.03).
>
> **Vorher (superseded):** Normierung durch die FLAECHE A = N·√3/2 -> Υ(0) = J·√3 ≈ 1.732.
> Das ueberschaetzte T_BKT um **~15%** und war durch die zu laxe Toleranz tol=0.15
> nie hart geprueft. BKT-Sprung-Kriterium: Nelson-Kosterlitz Υ(T_BKT) = 2·T_BKT/π.

### T_BKT-Mess-Stand (per-Site, Neumessung 2026-06-04)

Wolff-Sampling auf dem Dreiecks-Torus, Nelson-Kosterlitz-Crossing
(`src/260604 PHY030 ...`, Evidenz: `results/260604 PHY030 ... report.txt`):

| Gitter | T_BKT (per-Site, neu) | Referenz (arXiv:2501.07388) | Methode |
|---|---|---|---|
| triangular | **1.4007 ± 0.0081** (Weber-Minnhagen-Log-Fit, v02) | 1.418 | Wolff + WM-Log-Korrektur, chi^2-Fit |
| triangular | ~1.42 (v01-Schranken: L=19-Crossing 1.456; 1/lnL-Extrap. 1.382) | 1.418 | Wolff + NK-Crossing, FSS (konservativ) |
| square | 0.893 (PHY028, V&V) | 0.893 | Sandvik-Paar (C-frei) |
| honeycomb | nicht neu gemessen | 0.573 | offen |

**v02 (SOTA-Praezisierung, 2026-06-04):** Der Weber-Minnhagen-Log-Korrektur-Fit
`Υ(T_BKT,L) = (2 T_BKT/π)·(1 + 1/(2 ln L + C))` (Weber & Minnhagen PRB 37,
5986(R) (1988); Methodik arXiv:1302.2900) liefert einen Punktschaetzer
**1.4007 ± 0.0081 (−1.22%)** — naeher an der Referenz 1.418 als beide
konservativen v01-Schranken. Querscheck per C-freier Paar-Methode: 1.3801.
Beide v02-Schaetzer liegen knapp UNTER 1.418; bei nur drei kleinen Gittern
(L=9..19) sind sub-leading Log-Korrekturen eine reale Praezisions-Grenze
(ehrlich ausgewiesen). **v01 bleibt gueltig** als konservative Schranken-Methode.
Code: `src/260604 PHY030 ... v02 wm-logfit.py`, Evidenz:
`results/260604 PHY030 v02 wm-logfit report.txt`.

**Superseded:** alle vor 2026-06-04 publizierten T_BKT-Zahlen, die auf der
Flaechen-Normierung beruhen (u.a. PHY026-Report "Υ(T→0)=J·√3"), sind durch
diesen Mess-Stand ersetzt. Details in `CHANGELOG.md`.

## Kern

- **Engine:** `src/260602 PHI HEX core v2 2 hardened.py` (Kern) + PHY024-030-Experiment-Serie
- **Evidenz:** Selftest-Report + PHY025-030 Reports in `results/`; Testsuite in `tests/`; methodischer Audit (4 Responses) in `spec/`. Mess-Stand T_BKT(triangular) = 1.42 ± 0.04 (per-Site, 2026-06-04; Referenz 1.418).

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

Die Selftests laufen direkt in der Engine (`python <engine>.py`; Details im
Engine-Header). Die Testsuite (numpy/scipy/pytest noetig):

```
pytest                 # volle Suite (inkl. slow Mess-Lauf), 14/14 PASS
pytest -m "not slow"   # nur schnelle Korrektheits-Gates (inkl. WM-Synthetik-Orakel)
python "src/260604 PHY030 triangular tbkt per-site v01.py"            # T_BKT-Schranken (konservativ)
python "src/260604 PHY030 triangular tbkt per-site v02 wm-logfit.py"  # T_BKT Weber-Minnhagen-Log-Fit
```

CI prueft Syntax (`compileall`) + Lint (non-blocking) — die vollen
physikalischen Gate-Runs brauchen numpy/scipy und laufen lokal.

## Verwandte Repos

`coworker-dde` (kanonische v3-Engine-Quelle, 3DDE-Kanon) | `HQST` (p1_tasks-Floquet-Linie) |
`phi-hex` `hex-hqst` `u2-dualspec` `u6-bifurcation` `hexa-ntk` (diese Serie, 2026-06-04)

---
*Coworker Research / Coworkerz | Repo-Anlage 2026-06-04 nach arbeitsschablone_forschungs-repo-anlage*
