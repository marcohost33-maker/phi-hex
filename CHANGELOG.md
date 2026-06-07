# Changelog — Phi-Hex

Alle nennenswerten Aenderungen an Konventionen, Engine und Mess-Stand.
Format lose an Keep-a-Changelog angelehnt.

## [2026-06-07] PHY031 — T_BKT(honeycomb) gemessen + Lint vollstaendig

Schliesst die als "offen" markierte Luecke T_BKT(honeycomb) und zieht die
Lint-Baseline auf den vollen pycodestyle-Default (E701/E702) durch.

### Added
- `src/260607 PHY031 honeycomb tbkt per-site v01.py`: T_BKT(honeycomb) per-Site.
  Honeycomb-Torus (2-atomige Basis A/B, z=3, N=2L^2), gitter-agnostische
  Wolff-Pipeline (wiederverwendet PHY026 `wolff_sweep` + core-Helicity), C-freie
  Sandvik-Paar-Schaetzung (L,2L) am Weber-Minnhagen-Punkt + per-L NK-Crossings.
- `results/260607 PHY031 honeycomb tbkt per-site report.txt`: Gate-Evidenz
  (deterministisch seed=42).
- `tests/test_tbkt_honeycomb.py` (+8): Geometrie-Orakel `Upsilon(0)=3/4 J`
  (analytisch exakt, z=3, groessen-unabhaengig), Koordinations-/Bondlaengen-
  Check, Sandvik-Synthetik-Orakel (rekonstruiert bekanntes T_true, fails-
  before-faehig) + slow Crossing-Smoke. conftest/SOURCES.md ergaenzt.

### Result (Evidenz: `results/260607 PHY031 honeycomb tbkt per-site report.txt`)
| Methode | T_BKT(honeycomb) | Abw. vs 0.573 |
|---|---|---|
| Sandvik-Paar (L=12,24), C-frei | **0.595** | +3.90% |
| Sandvik-Paar (L=6,12), C-frei | 0.602 | +5.06% |
| per-L NK-Crossing L=24 | 0.628 | +9.53% |
| per-L NK-Crossing L=6 | 0.678 | +18.31% |

Gate A (perfekt ausgerichtet, T->0): `Upsilon(0) = 0.750000` EXAKT (analytischer
per-Site-Spinwellen-Grenzwert 3/4 J; vgl. Dreieck 3/2, Quadrat 1). Die per-L-
Crossings driften monoton von oben gegen die Referenz 0.573 (textbook BKT-
Finite-Size); das groesste Sandvik-Paar (12,24) verankert T_BKT C-frei auf
+3.9%. Bei kleinen L ehrlich ausgewiesene Finite-Size-Grenze.

### Changed
- `ruff.toml`: E701/E702-Ignore entfernt — die kompakten `clean()`-JSON-
  Serialisierer (core/PHY026-029) + inline-`;`-Statements (PHY028/029) sind
  auf die mehrzeilige Form (wie PHY030/031 `_clean`) ausformatiert; der Stack
  ist jetzt gegen den vollen pycodestyle-Default sauber (verhaltensneutral).

## [2026-06-07] Code-Review-Haertung — CI-Test-Gating, Lint-Baseline, Dead-Code

Best-Practice-Welle (GitHub-Actions/pytest + ruff), keine Physik-Aenderung.
Mess-Zahlen unveraendert.

### Added
- `.github/workflows/ci.yml`: neuer `test`-Job laeuft die schnellen
  Korrektheits-Gates (`pytest -m "not slow"`, 13 Tests) ueber die Matrix
  Python 3.10/3.11/3.12 (`fail-fast: false`, pip-Cache). Bisher lief in CI
  KEIN Test — die Physik-Korrektheits-Gates konnten Regressionen nicht
  fangen (Luecke gegen das Reality-Anchor-Prinzip).
- `requirements-dev.txt`: reproduzierbare Dev-/CI-Abhaengigkeiten (numpy/
  scipy/pytest/ruff), Cache-Key-Quelle fuer `setup-python`.
- `ruff.toml`: explizite Lint-Baseline (die der CI-Kommentar als naechsten
  Schritt ankuendigte). E701/E702 bewusst ignoriert (intentionaler kompakter
  Stil der JSON-`clean()`-Serialisierer und der Acklam-`_norm_ppf`).

### Changed
- `ci.yml` ruff-Schritt von non-blocking (`|| true`) auf **blockierend**;
  Lint-Job getrennt vom Test-Job.
- README/Status: Testanzahl ehrlich von "14/14" auf **15/15** korrigiert
  (lokal verifiziert 2026-06-07; pytest sammelt 15 Items — die
  parametrisierte WM-Orakel-Funktion zaehlt 3-fach).

### Fixed (Dead-Code, ruff F401/F841/F541/E401 — verhaltensneutral)
- `core` `bca_bootstrap_ci`: ungenutzte `combined`/`full` entfernt.
- Ungenutzte Importe/Variablen in PHY026/PHY027/PHY028/PHY029, PHY024_R0,
  phy017 (tote `typing`-Importe, leerer f-string, Multi-Import-Zeile).
- Provenance: Aenderungen sind repo-native Reviews (Praezedenz PR #1);
  `SOURCES.md`-Eintraege bleiben unangetastet (append-only Lineage).

## [2026-06-04] PHY030 v02 — Weber-Minnhagen-Log-Korrektur-Fit (SOTA-Praezisierung)

Praeziserer T_BKT-Punktschaetzer ergaenzend zu den konservativen v01-Schranken.
v01 bleibt unveraendert gueltig.

### Added
- `src/260604 PHY030 triangular tbkt per-site v02 wm-logfit.py`: T_BKT via
  Weber-Minnhagen-Finite-Size-Form `Υ(T_BKT,L) = (2 T_BKT/π)·(1 + 1/(2 ln L + C))`
  (Weber & Minnhagen, Phys. Rev. B 37, 5986(R) (1988); Methodik arXiv:1302.2900,
  arXiv:2406.12076). Zwei Auswertungen derselben Wolff-Daten:
  - **(A) Fixed-T-Least-Squares:** je T 1-Parameter-Fit (C) gewichtet 1/sem^2;
    `T_BKT = argmin_T chi^2(T)`, parabolische Verfeinerung, 1-sigma aus chi^2_min+1.
  - **(B) Paar-C-Elimination** (C-frei, arXiv:1302.2900): T_BKT(L1,L2) via
    `1/R(T,L1) − 1/R(T,L2) = 2 ln(L2/L1)`.
  - Gleiche Pipeline/Params wie v01 (radii 4/6/9, n_measure=400, n_burn=300,
    n_seeds=4, master_seed=42), feineres T-Gitter (Δ=0.01, 1.36–1.46).
- `tests/test_tbkt_wm_logfit.py` (+7 Tests): Synthetik-Orakel (rekonstruiert
  bekanntes T_true aus WM-konformen Daten, fails-before-faehig) + Mini-Smoke.

### Result (Evidenz: `results/260604 PHY030 v02 wm-logfit report.txt`)
| Methode | T_BKT | Abw. vs 1.418 |
|---|---|---|
| v02 WM-Fit (A) | **1.4007 ± 0.0081** | −1.22% |
| v02 Paar-Mittel (B) | 1.3801 | −2.67% |
| v01 obere Schranke (L=19-Crossing) | 1.456 | +2.68% |
| v01 untere Schranke (1/lnL-Extrap.) | 1.382 | −2.54% |

Der WM-Fit liegt naeher an der Referenz als beide v01-Schranken (Gate <3% PASS),
Residuen je L < 0.002. Beide v02-Schaetzer liegen knapp UNTER 1.418 — bei nur
drei kleinen Gittern (L=9..19) sind sub-leading Log-Korrekturen eine reale
Praezisions-Grenze (arXiv:1302.2900), ehrlich ausgewiesen. Testsuite 14/14 PASS.

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
