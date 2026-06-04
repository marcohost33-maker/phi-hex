# Changelog — Phi-Hex

Alle nennenswerten Aenderungen an Konventionen, Engine und Mess-Stand.
Format lose an Keep-a-Changelog angelehnt.

## [2026-06-04] Helicity per-Site + T_BKT-Neumessung

Wissenschafts-Audit 2026-06-04 (Befunde S2/S3), Marco-Entscheid 2026-06-04.

### Changed (KONVENTIONSWECHSEL — BREAKING fuer Mess-Zahlen)
- **Helicity-Modulus jetzt per-Site normiert** (Division durch Site-Zahl N)
  statt durch die FLAECHE A = N·√3/2 (Audit S2).
  - Spinwellen-Grenzwert Dreiecksgitter: **Υ(0) = 1.5 J** (vorher J·√3 ≈ 1.732).
  - Literatur-Standard (arXiv:2406.12076); konsistent mit der per-Site-Normierung
    in PHY028 (Quadratgitter, Υ(0)=1.0 J).
  - BKT-Sprung-Kriterium unveraendert: Nelson-Kosterlitz Υ(T_BKT) = 2·T_BKT/π.
- `helicity_from_ensemble` (core): Parameter `area_per_site` entfernt; teilt nun
  durch `n_nodes`.
- `validation_gates` (core): T_BKT-Toleranz `tol_T_bkt_rel` 0.15 → **0.03** geschaerft.
- PHY026 Teil-3-Gate: Referenz `J·√3` → `1.5·J` (per-Site).
- `T_BKT_REFERENCE` (core): Provenance arXiv:2501.07388 ergaenzt; honeycomb 0.575 → 0.573.

### Added
- `src/260604 PHY030 triangular tbkt per-site v01.py` — T_BKT(triangular)-
  Neumessung (Wolff + Nelson-Kosterlitz-Crossing, per-Site).
- `results/260604 PHY030 triangular tbkt per-site report.txt` — Gate-Evidenz.
- `tests/` — Korrektheits-Gates (pytest):
  - `test_helicity_per_site.py` — Υ(T→0)=1.5J, tol 0.03, groessen-unabhaengig.
  - `test_cross_validation_core_phy028.py` — core↔PHY028 gemeinsame per-Site-Defn (S3).
  - `test_tbkt_triangular.py` — NK-Crossing im physikalischen Band.
  - `conftest.py` — Modul-Loader fuer die Engine-Dateien mit Leerzeichen.
- `pytest.ini` — slow-Marker.

### Mess-Stand T_BKT (per-Site, NEU; ersetzt alle frueheren Zahlen)
Wolff-Sampling Dreiecks-Torus, L=9/13/19, n_measure=400, n_seeds=4, seed=42:

| Gitter | T_BKT (per-Site) | Referenz (arXiv:2501.07388) | Abweichung |
|---|---|---|---|
| triangular | ~1.42 ± 0.04 (L=19-Crossing 1.456; 1/lnL-Extrap. 1.382) | 1.418 | Ref liegt zwischen den Schaetzern |
| square | 0.893 (PHY028, V&V Sandvik-Paar) | 0.893 | <1% |
| honeycomb | nicht neu gemessen | 0.573 | offen |

Die Per-L-Crossings driften mit wachsendem L monoton von oben (1.483 → 1.472
→ 1.456) gegen die Referenz; die 1/ln(L)-Extrapolation liefert die untere
Schranke (1.382). Der Referenzwert 1.418 liegt zwischen beiden — physikalisch
konsistente BKT-Konvergenz. Unsicherheit ehrlich: kleine L + endliche
Wolff-Statistik; ein engeres 3%-Resultat braucht groessere L.

### Superseded
- Alle vor 2026-06-04 publizierten T_BKT-Zahlen, die auf der Flaechen-
  Normierung beruhen (Faktor √3/(2·… ) zu hoch, ~15% T_BKT-Ueberschaetzung).
- PHY026-Report-Aussage "Υ(T→0) = J·√3 = 1.732" — ersetzt durch 1.5·J per-Site.

## [2026-06-04] P0 edges-Fix (PR #1)
- `edges = sorted(edges_set)` im open-boundary-Pfad restauriert (war beim
  v2.0→v2.2-Hardening geloescht; Default-Selbsttest lief as-shipped nie). Audit S1.

## [2026-06-04] Initial-Anlage
- Forschungs-Repo-Schablone, SHA-Provenance in SOURCES.md, Haertungs-Welle (CI/SHA-Pins).
