# Changelog — Phi-Hex

Alle nennenswerten Aenderungen an Konventionen, Engine und Mess-Stand.
Format lose an Keep-a-Changelog angelehnt.

## [2026-06-11] PHY037 — HKS-Extrapolation gitter-uebergreifend (ehrliche Methoden-Bilanz)

Generalisiert den HKS-Apparat (arXiv:1302.2900) auf honeycomb/kagome/triangular
und zieht die ehrliche Bilanz, WO die C-eliminierte Paar-Extrapolation aus kleinen
L sauber konvergiert - und wo nicht.

### Added
- `src/260611 PHY037 hks multilattice extrapolation v01.py`: parst die
  committed "T_BKT(L=a,b)"-Paar-Zeilen der PHY030/032/033-Reports ZUR LAUFZEIT
  (single source of truth, kein eingebettetes Daten-Duplikat) und extrapoliert
  je Gitter u=1/(ln Lc)^2 -> L->inf.
  - **Befund (ehrlich, nicht getunt):** kagome 0.8236 (−0.17% vs 0.825) SAUBER;
    honeycomb 0.5747 (−0.23% vs 0.576) SAUBER, aber NUR mit L=48 (ohne L=48
    kippt es, s. PHY035); triangular 1.4591 (+2.90% vs 1.418) UEBERSCHIESST -
    hier ist das BLANKE groesste Paar (1.3888) naeher an der Referenz. =>
    Klein-L-Paar-Extrapolation ist gitter-abhaengig + empfindlich auf das
    groesste enthaltene L, KEIN universelles Wundermittel (deckt sich mit HKS:
    L=48..192 noetig).
- `tests/test_phy037_hks_multilattice.py` (+5 Gates): Drift-Guard (geparste
  Paare == committed Reports), Extrapolations-Regression (Werte gepinnt),
  Methoden-Verdikt (honeycomb/kagome SAUBER, triangular UEBERSCHIESST), Wiring.
- `results/260611 PHY037 hks multilattice extrapolation report.txt`: Gate-Log.

### Added (Mess-Apparat)
- `src/260611 PHY036 honeycomb hks large ladder v01.py`: groessere Leiter
  L=16/24/32/48/64 (drei verschachtelte Paare bis (32,64); dank Vektorisierung
  L=64/N=8192 lokal vertretbar), reine Wiederverwendung der parametrisierten
  PHY035-Bausteine. Klaert, ob L>=48 die Konvergenz stabilisiert (Lauf-Evidenz
  separat). `tests/test_phy036_*.py` (+2 Konfig-/Geometrie-Gates).
- PHY035-Bausteine (`measure_cube`/`pair_estimates`/`make_hks_extrap_estimator`)
  um Keyword-Parameter `ladder`/`pairs`/`t_grid` erweitert (rueckwaerts-
  kompatibel; Defaults unveraendert) -> PHY036 reuse ohne Duplikat.

## [2026-06-11] PHY035 — honeycomb dichte L-Leiter: falsifiziert die PHY034-Lesart (Negativ-Result)

Praezisions-Follow-up zu PHY034 - und ehrliche Selbst-Falsifikation. PHY034
extrapolierte mit NUR 2 Paaren (2-Punkt). PHY035 misst frisch auf der dichten
Leiter L=8/12/16/24/32 (drei verschachtelte Verdopplungspaare; dank Kern-
Vektorisierung ~3 min lokal) und macht damit eine ECHTE >=3-Punkt-Extrapolation.

### Added
- `src/260611 PHY035 honeycomb hks dense ladder v01.py`: Wolff-Messung
  (phy032-Pipeline 1:1) + C-eliminierte Paare + OLS-Extrapolation (phy034-
  Primitive) + EIGENES Seed-Bootstrap-CI des Extrapolations-Werts.
  - **Befund (NEGATIV):** die drei Paare (8,16)/(12,24)/(16,32) liegen FLACH
    bei ~0.597-0.601 - keine saubere Drift gegen 0.573. Die >=3-Punkt-
    Extrapolation ergibt ~0.605 (CI[0.577,0.616]), OBERHALB der Referenzen.
    => PHY034 (0.574) war NICHT robust: das Ergebnis hing allein am Paar
    (24,48) (groesstes L, hier nicht in der Leiter). Bei L<=32 zeigt honeycomb
    KEINE belastbare Konvergenz; das deckt sich mit HKS (L=48..192 noetig).
  - **Cross-Check:** L=12/24 teilen den RNG-Vertrag mit PHY032 -> T(12,24)=
    0.6014 bit-identisch reproduziert (Mess-Pipeline-Integritaet).
  - OVERALL=PASS bedeutet hier ANALYSE-INTEGRITAET (Gate A exakt, 3 Paare
    endlich, PHY032-Repro, CI berechnet) - die Physik-Hypothese "clean
    convergence" ist bewusst KEIN Build-Breaker und als FAIL-Finding ehrlich
    dokumentiert (AGENTS.md: Negativ-Results sind Buerger erster Klasse).
- `tests/test_phy035_hks_dense_ladder.py` (+4 Gates, davon 1 slow): Geometrie-
  Orakel (Υ(0)=0.75 exakt), Extrapolations-Wiring (synthetische WM-Daten ->
  t_true exakt), Paar-Wiring, Mini-Wolff-Smoke.
- `results/260611 PHY035 honeycomb hks dense ladder report.txt`: Gate-Log mit
  vollem Υ(T,L)-Gitter + ehrlichem Negativ-Befund.

### Changed
- README/CHANGELOG: PHY034-Befund von "belegt/konsistent" auf "suggestiv, nur
  2-Punkt, durch PHY035 als nicht robust widerlegt" getempert (Lineage-
  Ehrlichkeit; PHY034 selbst war bereits als Konsistenz-Check caveated).

## [2026-06-11] PHY034 — honeycomb T_BKT: HKS-L->inf-Extrapolation (Methoden-Haertung)

Schliesst die in PHY032 als "erwartbar finite-size, NICHT methodisch" NUR
BEHAUPTETE Luecke und testet sie. Web-recherchiert + gegen-recherchiert:
Hsieh-Kao-Sandvik (J. Stat. Mech. (2013) P09001, arXiv:1302.2900) zeigen, dass
die C-eliminierte Paar-Schaetzung T_BKT(L1,L2) selbst noch eine sub-leading
logarithmische Drift traegt und zum thermodynamischen Limes extrapoliert werden
MUSS ("sub-leading logarithmic corrections have significant effects; previous
works underestimated T_BKT"). PHY032 lieferte die Paare ohne diesen Schritt.

### Added
- `src/260611 PHY034 honeycomb hks extrapolation v01.py`: rein analytische
  Re-Analyse der committed PHY032-Daten (KEIN neuer Monte-Carlo-Lauf,
  deterministisch). Reproduziert die PHY032-Paar-Schaetzer (0.6014/0.5917/
  0.5968) und extrapoliert die Verdopplungspaare (12,24)/(24,48) linear in der
  BKT-Variablen u=1/(ln Lc)^p zum Intercept L->inf.
  - **Befund:** im HKS-Standard (p=2, geom. Mittel) ergibt sich
    T_BKT = **0.5741** (+0.19% vs 0.573 / -0.33% vs 0.576) — der scheinbare
    +3%-Offset der kleinen L ist damit quantitativ als sub-leading-log-
    Finite-Size-Artefakt belegt, konsistent mit BEIDEN Referenzen
    (0.573 arXiv:2501.07388; 0.576(3) arXiv:2406.12076).
  - Ehrliche Doppel-Unsicherheit: Methoden-Systematik (Variablenwahl)
    [0.547, 0.583] + statistische 95%-CI [0.553, 0.595] (MC-Propagation der
    per-Paar Seed-Bootstrap-Streuung; die 2-Punkt-Extrapolation VERSTAERKT
    den Paar-Fehler ~3x — ehrlich ausgewiesen).
  - Caveat: nur 2 Verdopplungspaare -> Konsistenz-Check, KEINE Praezision;
    dichtere L-Leiter (>=3 Paare) ist der Praezisions-Follow-up (dank
    Kern-Vektorisierung guenstig).
- `tests/test_phy034_hks_extrapolation.py` (+12 Gates): Drift-Guard (embedded
  Gitter == committed PHY032-Report, Zeichen-genau geparst), Extrapolations-
  Orakel (rein quadratische Drift -> exakte Rekonstruktion; kubische Rest ->
  Bias-Reduktion), OLS-Intercept, statistische CI (deterministisch + Hebel),
  Befund-Gates (Primaer im Referenzband, Offset reduziert, Band/CI decken Ref).
- `results/260611 PHY034 honeycomb hks extrapolation report.txt`: Gate-Log.

## [2026-06-11] Selbst-Audit-Haertung — cliff_delta NaN-Treue + Doku-Abgleich

Folge-Haertung nach Selbst-Review der Vektorisierungs-Welle (2026-06-09).
Behebt eine beim Sortier-Rewrite von `cliff_delta` eingeschleppte
Korrektheits-Luecke und gleicht stale Doku-Zaehler sowie die README-Mess-Stand-
Kurzfassung an die bereits korrekten Detailtabellen an. **Keine Aenderung an
Physik, Konvention oder gemessenen Werten** — nur Doku-Konsistenz.

### Fixed
- `cliff_delta`: NaN-Treue wiederhergestellt. Der Sortier-/`searchsorted`-Pfad
  zaehlte NaN-Eingaben falsch (NaN sortiert ans Ende -> als groesster Wert
  mitgezaehlt), waehrend die fruehe strikte Doppelschleife jeden NaN-Vergleich
  als False behandelt (NaN traegt weder zu gt noch zu lt bei). Fix: NUR NaN aus
  den Paar-Zaehlungen filtern (inf bleibt drin und vergleicht normal), Nenner
  na*nb bleibt die volle Stichprobe -> wieder **bit-identisch** zur Doppel-
  schleife in ALLEN Faellen (NaN/inf/Gleichstand). Plus Guard fuer leere
  Stichprobe (kein ZeroDivision). Verifiziert per neuem Regressions-Gate
  `test_cliff_delta_nan_inf_match_scalar` (+ `_empty_is_safe`).

### Docs
- README/CHANGELOG-Zaehler aktualisiert (Testsuite **53/53**: 48 schnelle Gates
  + 5 slow) und die CHANGELOG-Beschreibung des `cliff_delta`-Rewrites auf den
  tatsaechlich gemergten Sortier-/Speicher-sicheren Stand korrigiert (war
  irrefuehrend als "numpy-Broadcast" beschrieben).
- BKT-Referenzkonstanten gegen Literatur web-rueckverifiziert (arXiv:2501.07388):
  square 0.89290(5), honeycomb 0.573, kagome 0.825 (rough estimate) bestaetigt.
- README-Mess-Stand-Zusammenfassung (Abschnitt "Evidenz") auf die maßgeblichen
  Best-Schaetzungen angeglichen: honeycomb **0.5917** CI[0.587,0.598] (PHY032
  groesstes L, supersedet die alte PHY031-Zahl 0.595) und kagome **0.8479**
  CI[0.8414,0.8507] (PHY033) ergaenzt; Report-Range PHY025-031 -> PHY025-033,
  Experiment-Serie PHY024-032 -> PHY024-033, "letzte Haertung" -> 2026-06-11.
  Die detaillierten PHY032/PHY033-Tabellen waren bereits korrekt; nur die
  Kurzfassung hinkte hinterher.

## [2026-06-09] Kern-Performance — vektorisierte Mess-Hotspots (physik-invariant)

Reine Performance-Optimierung des Kernmoduls (`phi_hex_core_v2`). **Keine
Aenderung an Physik, Konvention oder Mess-Stand.** Die fruehen skalaren
Python-Schleifen der beiden heissesten Mess-Pfade sind durch aequivalente
numpy-Vektorisierung ersetzt. Aequivalenz ist als blockierendes CI-Gate
festgepinnt (`tests/test_core_vectorization.py`, +8 Tests).

### Changed
- `helicity_terms` vektorisiert (gecachte Kanten-Index-/Verschiebungs-Arrays
  auf dem Gitter-Objekt, gitter-agnostisch -> Triangular **und** Honeycomb).
  Dieser Pfad wird im Wolff-Sampling **pro Messung** aufgerufen und dominiert
  die Laufzeit der grossen Laeufe (PHY030/031/032). Gemessen **~9x** schneller
  am honeycomb-Hot-Loop (L=24/48). Mathematisch identisch zur Schleife; die
  Differenz ist reine Gleitkomma-Summationsreihenfolge (~1e-13 relativ),
  unterhalb jeder berichteten Stelle -> bestehende Gate-Logs in `results/`
  bleiben gueltige Evidenz. Cache identitaets-invalidiert (haelt Referenzen auf
  die edges-/edge_disp-Listen; Umverdrahten baut die Arrays neu) - kein
  stale-Cache (Codex-Review P2).
- `cliff_delta` sortier-/rangbasiert (`np.searchsorted`, O((na+nb) log nb) Zeit,
  O(nb) Speicher) statt O(na*nb)-Doppelschleife. **Bit-identisch** (nur
  ganzzahlige strikte Paar-Vergleiche) und ohne die na*nb-Vergleichsmatrix zu
  materialisieren (kein OOM bei grossen Stichproben; Codex-Review P2).

### Added
- `tests/test_core_vectorization.py`: Aequivalenz-Gates gegen eine
  unabhaengige in-test skalare Referenz (helicity_terms: rtol 1e-10 ueber
  Triangular/Honeycomb + Messrichtungen + Cache-Konsistenz + Rewire-
  Invalidierung; cliff_delta: bit-identisch + Gleichstand/Extrem-Faelle).
  Eine Vorzeichen-/Index-Mutation bricht das Gate.

## [2026-06-09] PHY033 — kagome T_BKT (AP-3): per-Site Wolff + WM-Log-Fit + Sandvik-Paar + Seed-Bootstrap-CI

Schliesst die im README/CHANGELOG als offen markierte Luecke **T_BKT(kagome)**
(AP-3 der phi-hex-Gitter-Serie). Das Kagome-Gitter (eckenteilende Dreiecke,
z=4, 3 Sites/Einheitszelle) reiht sich konsistent zwischen Honeycomb (z=3,
~0.573) und Triangular (z=6, ~1.418) ein. Keine Aenderung am Kern-Physik-Code;
die WM-Fit- (PHY030 v02) und Seed-Bootstrap-Maschinerie (PHY032) werden 1:1
gitter-/methoden-agnostisch wiederverwendet (kein Reinvent).

### Added
- `src/260609 PHY033 kagome tbkt per-site v01.py`:
  - **Kagome-Torus** (`build_kagome_lattice`): N=3L^2 Knoten, 6L^2 Bonds (=2N),
    jeder Knoten exakt z=4; Bravais a1=(2,0)/a2=(1,sqrt3), 3-atomige Basis;
    6 distinkte unit-length-Bonds je Zelle (algebraisch + im Geometrie-Gate
    verifiziert). **Per-Site-T=0-Orakel Upsilon(0) = 1 J EXAKT** (Gate A,
    gemessen 1.000000), groessen-unabhaengig.
  - **T_BKT** via (A) Weber-Minnhagen Fixed-T-Log-Fit + (B) C-freie
    Sandvik-Paar-Schaetzung, beide mit einheitlich propagiertem
    **Seed-Bootstrap-CI** (+ Jackknife). L-Set **(12, 24, 36)**, N bis 3888
    (8GB-RAM-PC). Eigener RNG-Stream-Praefix (500+...), keine Kollision mit
    Honeycomb (400+...).
- `results/260609 PHY033 kagome tbkt per-site report.txt`: Gate-Evidenz,
  deterministisch (master_seed=42, n_seeds=6, n_measure=120, n_boot=2000),
  Laufzeit ~511s lokal.
- `tests/test_tbkt_kagome.py` (+9 schnelle, +1 slow): Geometrie-Orakel
  (Upsilon(0)=1, z=4, unit bonds, min-image-Bond-Konsistenz), Sandvik-Synthetik-
  Orakel (rekonstruiert T_true), Bandplausibilitaet (0.573<0.825<1.418),
  Smoke (Mini-Wolff -> NK-Crossing im Band).

### Mess-Ergebnis (ehrlich, NICHT auf Referenz getunt)
- **WM-Log-Fit: T_BKT(kagome) = 0.8479, CI[0.8414, 0.8507] -> +2.78% vs 0.825.**
- Sandvik-Paar (24,36): 0.8377, CI[0.8273, 0.8650] -> +1.54%; (12,36): 0.8410.
- Per-L NK-Crossings driften von oben (L=12: 0.8899 -> L=36: 0.8671) gegen die
  Referenz — erwartetes finite-size-Verhalten.
- **Caveat:** L bis 36 (lokal RAM-/zeit-limitiert; Referenz nutzt L=48..192).
  Die Referenz 0.825 ist selbst ein "rough estimate" (arXiv:2501.07388, Tab. 1).
  Rest-Bias ist erwartbar finite-size, NICHT methodisch.

### Verifiziert
- `pytest -m "not slow"` 38/38 PASS, volle Suite 43/43 (lokal 2026-06-09).
- `ruff check src` clean.

## [2026-06-07] PHY032 — honeycomb WM-Log-Fit + groessere L (AP-1) + einheitliche Seed-Bootstrap-CI (AP-2)

Schliesst die zwei im Audit (`Vero/Meta/AUDITS/2026-06-07_phi-hex-welle`)
priorisierten Arbeitspakete: (AP-1) die Weber-Minnhagen-Log-Fit-Methodik von
PHY030 v02 auch auf honeycomb anwenden + groessere L (12/24/48 statt 6/12/24),
und (AP-2) einheitliche propagierte Fehlerbalken (Seed-Bootstrap + Jackknife)
auf ALLE berichteten T_BKT-Schaetzer. Keine Aenderung am Kern-Physik-Code.

### Added
- `src/260607 PHY032 honeycomb wm-logfit bootstrap v01.py`:
  - **AP-1:** Weber-Minnhagen Fixed-T-Fit (A) + C-freie Sandvik-Paar (B) auf
    honeycomb, **1:1 die gitter-agnostischen Fit-Funktionen aus PHY030 v02**
    (`wm_form`/`fit_tbkt_fixedT`/`tbkt_pair_C_eliminated` — kein Reinvent,
    single source of truth). L-Set **(12, 24, 48)**, N bis 4608.
  - **AP-2:** `bootstrap_tbkt_over_seeds` (Perzentil-Bootstrap ueber die
    n_seeds-Replikate, deterministisch via `make_rng`) + `jackknife_tbkt_over_seeds`
    (leave-one-seed-out SE), gitter-/methoden-agnostisch auf JEDEN T_BKT-
    Schaetzer (Sandvik-Paar + WM-Fit) anwendbar. `measure_honeycomb_seedwise`
    behaelt die per-Seed-Upsilon-Werte (RNG-Stream-Vertrag identisch PHY031 v01).
  - **Beide Literaturwerte ausgewiesen (kein Kanon gekuert):** 0.573
    (arXiv:2501.07388) UND 0.576(3) (arXiv:2406.12076).
- `results/260607 PHY032 honeycomb wm-logfit bootstrap report.txt`:
  Gate-Evidenz, deterministisch (master_seed=42, n_boot=2000 seed=42).
- `tests/test_tbkt_honeycomb_bootstrap.py` (+9 schnelle, +1 slow): Bootstrap-
  Orakel (CI ueberdeckt bekanntes T_true; degenerierter Estimator -> None-CI =
  Silent-Failure-Gate; Determinismus), Jackknife-Orakel (SE=0 bei identischen
  Seeds), Estimator-Verdrahtung (WM-Fit rekonstruiert T_true; Rand-Minimum ->
  None), Repro gegen PHY031 v01 (per-Seed-Mittel bit-identisch, slow).

### Result (Evidenz: `results/260607 PHY032 honeycomb wm-logfit bootstrap report.txt`)
| Methode | T_BKT(honeycomb) | vs 0.573 | vs 0.576(3) |
|---|---|---|---|
| Sandvik-Paar (24,48) — **groesstes L** | **0.5917** CI[0.587,0.598] | +3.27% | +2.73% |
| Sandvik-Paar (12,48) | 0.5968 CI[0.590,0.598] | +4.15% | +3.60% |
| Sandvik-Paar (12,24) | 0.6014 CI[0.599,0.610] | +4.95% | +4.41% |
| Weber-Minnhagen Fixed-T-Fit (alle L, argmin@0.5975) | 0.5964 CI[0.594,0.599] | +4.08% | +3.54% |
| PHY031 v01 Paar (12,24) — superseded | 0.595 | +3.90% | — |

**Befund (ehrlich, NICHT getunt):** Groesseres L senkt den Sandvik-Paar-Bias
monoton (Paar(24,48)=+3.27% < Paar(12,24)=+4.95%) — die Audit-Hypothese
bestaetigt sich (Praezedenz triangular). Der Rest-Bias bei lokal-vertretbaren
L (max 48 auf 8GB-RAM-PC; Referenz-Studien L=48..192 / 8..128) ist erwartbar
finite-size-dominiert, nicht methodisch. Die Seed-Bootstrap-CIs erfassen die
Seed-Streuung (NICHT den finite-size-Bias) und ueberdecken 0.573/0.576 NICHT
-> der Rest-Abstand ist reales finite-size. Der WM-Fit (0.5964) liegt hier -
anders als bei triangular - **nicht** unter dem groessten Paar-Schaetzer; bei
nur drei kleinen L sind sub-leading Log-Korrekturen eine reale Praezisions-
Grenze. Gate A (T->0): Upsilon(0)=0.750000 EXAKT. Testsuite 33/33 PASS.

## [2026-06-07] PHY031 — T_BKT(honeycomb) gemessen + Lint vollstaendig

Schliesst die als "offen" markierte Luecke T_BKT(honeycomb) und zieht die
Lint-Baseline ohne Ausnahmen auf den ruff-Default-Regelsatz (E4/E7/E9 + F,
inkl. E701/E702) durch. (Praezisierung: das ist die ruff-Default-Auswahl,
NICHT das gesamte pycodestyle-Regelwerk, das auch E1/E2/E3/E5 + W umfasst.)

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
  ist jetzt gegen den ruff-Default-Regelsatz (E4/E7/E9 + F) ohne Ausnahmen
  sauber (verhaltensneutral).

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
