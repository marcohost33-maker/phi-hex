# Phi-Hex

> **Stand 2026-09-28 (W4-v03 / PHY049/PHY050):** Jiang/PTEP primaertextlich nachgeprueft: der publizierte Honeycomb-Helicity-Modulus enthaelt den Faktor 4/(3 sqrt(3)); die per-Site-Erklaerung fuer 0.576(4) ist falsifiziert. Der normierungsfreie Correlation-Ratio-Quercheck ist vorregistriert; G4 ist deterministisch implementiert und testgegated, aber ein committierter `results/`-Gate-Log fehlt noch, daher wird hier noch kein repository-level VALIDATED-Status behauptet. Produktion ist nur nach gruenem Runtime-Preflight/VAL-BIT/G1-G4 zulaessig; Physikinterpretation bleibt bis zu realem G0/G5/G6-Evidenzlauf gesperrt.  
> **Status:** Forschungs-Repo (oeffentlich seit 2026-08-17) | XY/BKT-Physik auf Dreiecks-, Honeycomb- und Kagome-Gittern.  
> **Lizenz:** Apache-2.0 | **Lineage/Provenance:** siehe `SOURCES.md`.  
> **Stand 2026-09-26b (Selbstpruefung PR #50):** **Audit O1 entschieden - das NK-Kriterium verlangt Upsilon PRO FLAECHE** (PHY045: exakt + MC, per Site mit 22-54 sigma widerlegt). Alle per-Site-T_BKT-Werte auf triangular/honeycomb/kagome in diesem README sind damit konventions-verzerrt (Lineage bleibt; Neuauswertung separat). W4 laeuft nach Vorregistrierung v02 (Upsilon pro Flaeche, Wolff, HKS).  
> **Stand 2026-09-26 (Issue #45, PR #50):** Integritaets-Reparatur (PHY042-Domaenen-Semantik, Referenzband-Provenienz, PHY043-Text-Drift), W4-Vorregistrierung (vor der Kalibrierung committet) und blinde WL-Kalibrierung PHY044 mit bit-identischem Numba-Kernel: **W4-GO** mit Rezept "4x Produktion" (ohne Marge), kein T_BKT-Wert (Blind-Vertrag). Details: Abschnitt PHY044 unten.  
> **Aktueller Review-Stand:** Review-Nachtrag 2026-08-08 (`spec/260710 ... v01.md` §4): Messpipeline erneut gegen unabhaengige Orakel defektfrei; PHY040-M3-Guard, Lint-Baseline-Pin, O8/O9 inventarisiert. PHY043 (Audit O1): konventionsfreier Quercheck triangular — qualitativ konsistent mit der Referenz-Lage, keine 1%-Diskriminierung (siehe unten). Davor: Code-Audit 2026-07-10: P0-Geometrie-Fix im Quadratgitter (PHY028/039/040 neu gerechnet, V&V-Anker ehrlich auf ~1% herabgestuft), Pol-Guards in allen Paar-Schaetzern; PHY042/PR #20 als Pipeline-Finding, **kein neuer T_BKT-Bestwert** (Grenzen NR-PHY042-02/03).

Phi-Hex untersucht das 2D-XY-Modell und BKT-Physik auf periodischen Gittern: Helicity-Modulus, Nelson-Kosterlitz-Sprung, Wolff-Cluster, Wang-Landau-DOS und Finite-Size-Scaling.

## Was Phi-Hex ausdruecklich NICHT ist

- **Keine publizierte Arbeit und kein Referenzwert.** Die Tabelle unter *Mess-Stand*
  zeigt **interne** Staende, eingeordnet relativ zu Literaturankern — kein neuer
  T_BKT-Bestwert. Ein Bestwert-Claim braucht Cross-Family-Review (siehe letzter Abschnitt).
- **Keine Bibliothek.** Kein PyPI-Paket, keine stabile API, keine SemVer-Zusage; die
  Module sind Experiment-Skripte mit Datums-Praefix im Dateinamen.
- **Kein allgemeiner MC-Baukasten.** Der Code ist auf XY/BKT auf den hier genannten
  Gittern zugeschnitten.
- **Gruene Gates sind kein Physik-Beweis.** Die Pipeline-Gates pruefen Integritaet
  (Determinismus, Orakel, Rand-Leak, Domaenen-Guards), nicht die Richtigkeit einer
  physikalischen Aussage. Wo ein Gate ehrlich FAIL meldet, steht das so im README.
- **Kein Ersatz fuer die Specs.** Jede Zahl hier hat ihren Vertrag in `spec/` und ihren
  Lauf in `results/`; ohne diese beiden ist eine Zahl aus diesem README nicht zitierfaehig.

## Konvention: Helicity-Modulus — NK-Claims pro Flaeche

Audit O1/PHY045 hat die fruehere README-Konvention korrigiert: fuer den
Nelson-Kosterlitz-Sprung ist der Helicity-Modulus in der physikalischen
Flaechennormierung zu verwenden. Per-Site-Werte bleiben als historische
Lineage und als interne Geometrie-/Implementierungsorakel erhalten, duerfen
aber nicht unveraendert in einen NK-T_BKT-Claim eingehen.

Wichtige interne Orakel (per Site, **keine** NK-Claim-Normierung):

- Triangular: `Upsilon(0) = 1.5 J` per Site.
- Honeycomb: `Upsilon_2(0) = 3/4 J` per Site.
- Honeycomb ground state: `E0/N = -3J/2`.

## Honeycomb-Referenzen: Referenzband, kein Einzelanker

Ab PR #18/PHY041 gilt `spec/260703 PHI HEX honeycomb reference conventions audit v01.md` als Vertragsquelle. Honeycomb-Werte werden nicht mehr als ein einzelner harter Referenzwert gefuehrt, sondern als Band mit Quelle, Observable und beta/T-Konvention.

| Quelle | berichtete Groesse | T-Form fuer Vergleich | Rolle |
|---|---:|---:|---|
| arXiv:2501.07388 | T-Wert | 0.573 | Multi-Lattice-Anker |
| arXiv:2406.14812 (historische Pre-Journal-Version) | T_BKT,H = 0.571(8) | 0.571 +/- 0.008 | superseded_historical |
| arXiv:2406.14812 (historische Pre-Journal-Version) | T_BKT,H = 0.560(9) | 0.560 +/- 0.009 | superseded_historical |
| arXiv:2406.12076 (historische beta-Auswertung) | beta_BKT = 1.687(3) | 0.5928 +/- 0.0011 | superseded_historical |
| arXiv:2406.12076 (historische beta-Auswertung) | beta_BKT = 1.635(11) | 0.6116 +/- 0.0041 | superseded_historical |
| arXiv:2406.12076 (historische beta-Auswertung) | beta_BKT = 1.724(2) | 0.5800 +/- 0.0007 | superseded_historical |
| arXiv:2406.12076 (v4 / Phys. Scr. 100 065953) | T_BKT = 0.576(3) | 0.576 +/- 0.003 | Upsilon, Wang-Landau (direkt T) |
| arXiv:2406.12076 (v4 / Phys. Scr. 100 065953) | T_BKT = 0.575(8) | 0.575 +/- 0.008 | Upsilon, sim. annealing (direkt T) |
| PTEP 2024 103A02 (Journal-Fassung 2406.14812) | T_BKT = 0.576(4) | 0.576 +/- 0.004 | Helicity |
| PTEP 2024 103A02 (Journal-Fassung 2406.14812) | T_BKT = 0.572(3) | 0.572 +/- 0.003 | NN |

Konversion: `T = 1 / beta`, `sigma_T = sigma_beta / beta^2`.

**Provenienz-Nachtrag 2026-09-26 (Issue #45 §3):**
`spec/260926 PHI HEX honeycomb reference provenance addendum v01.md`.
0.576(3) ist der direkt berichtete WL-T-Wert von arXiv:2406.12076 (vorher als
"unattribuierter Legacy-Anker" gefuehrt) und stand bis dahin gar nicht im
Band; die Journal-Fassung von arXiv:2406.14812 berichtet andere Werte als
aeltere Versionen - beide bleiben getrennte Lineage. Die historischen
beta-Kanaele (0.580..0.612 in T-Form) sind fuer aktuelle Interpretation
**superseded_historical** und werden nicht mehr als aktive Referenzanker
verwendet. Jiang/PTEP sowie die direkten de-Andrade-v4-T-Werte wurden
primaertextlich verifiziert; aktuelle Auswertung filtert explizit nach dem
Provenienzstatus.

## Mess-Stand

**Konventions-Hinweis 2026-09-26 (Audit O1 entschieden, PHY045).** Die erste
Tabelle unten ist per Site normiert und auf triangular, honeycomb und kagome
konventions-verzerrt. Sie bleibt als Lineage stehen. Die per-Flaeche-
Neuauswertung derselben committeten Daten (PHY047, kleine L, ohne
Extrapolation; FINDING, kein Bestwert):

| Gitter | a_s | per Site (bisher) | per Flaeche (PHY047) | normierungsfreie Referenz |
|---|---:|---:|---:|---:|
| square | 1 | Paar (16,32) 0.8841 | identisch | 0.8929-0.8935 |
| triangular | 0.866 | WM 1.4008 (-4.4 %) | WM 1.4447 (-1.4 %) | 1.465 (HT-Reihe, arXiv:0806.1496) |
| honeycomb (WL, PHY042) | 1.299 | Paare 0.588-0.592 | Paare 0.558-0.560 | 0.573 (Corr.-Ratio, arXiv:2501.07388) |
| kagome | 1.155 | WM 0.850 (+3.0 %) | WM 0.819 (-0.7 %) | 0.825 (grob, arXiv:2501.07388) |

Pro Flaeche liegen alle vier Gitter bei kleinem L leicht UNTERHALB der
normierungsfreien Werte, im selben Vorzeichen wie die square-Kontrolle. Per
Site streute das gitterabhaengig um +/-3-4 %. Details:
`results/260926 PHY047 per-area reanalysis committed report.txt`.

| Gitter | aktueller interner Stand (per Site, Lineage) | externe Einordnung | Methode |
|---|---:|---|---|
| square | 0.8841 (Paar (16,32), -1.05% vs 0.8935) | V&V-Anker ~1% (nach edge_disp-Fix 2026-07-10; das fruehere +0.16% war Bug-Artefakt) | Sandvik-Paar / V&V, PHY028 |
| triangular | 1.4007 +/- 0.0081 | nahe 1.418 | Wolff + Weber-Minnhagen |
| honeycomb | 0.5917 CI[0.587,0.598] fuer Paar (24,48) | oberhalb der unteren Honeycomb-Anker, finite-size-sensitiv | Wolff + Sandvik-Paar / PHY032 |
| kagome | 0.8479 CI[0.8414,0.8507] | nahe rough estimate 0.825 | Wolff + WM-Fit / PHY033 |

**Audit-Hinweis 2026-07-10 (square-Familie):** Der Wrap-Bond-Vorzeichen-Bug
in `build_square_lattice` (PHY028, per Import auch PHY039/PHY040) hat
Upsilon auf dem Quadratgitter systematisch unterschaetzt; alle drei Module
wurden deterministisch (seed=42) neu gerechnet. Der PHY039-Y4-Dip liegt auf
den korrigierten Kurven AUSSERHALB des T-Fensters [0.85,0.95] (Dip-Gates
ehrlich FAIL; Fenster-Erweiterung = naechste Stufe). triangular/honeycomb/
kagome-Gitter waren nicht betroffen (Builder verifiziert). Details:
`spec/260710 PHI HEX code audit v01.md`.

PHY041 liefert fuer honeycomb mit Wang-Landau/1-t bei L<=24 die glatten Paarwerte:

| Paar | T_BKT |
|---|---:|
| (12,16) | 0.5951 |
| (12,24) | 0.6029 |
| (16,24) | 0.6087 |

Interpretation: Das ist ein **FINDING** zur Pipeline und zum kleinen-L-Verhalten, kein finaler T_BKT. Der Wert liegt im oberen Bereich des Referenzbands bzw. oberhalb der unteren Honeycomb-Anker. Die naechste Stufe (L=32/48 plus getrennte Upsilon_2-/Upsilon_4-FSS) ist mit PHY042 gelaufen (siehe unten).

## PHY042 — Voll-Lauf L=24/32/48 (Finding, kein Bestwert)

PHY042 fuehrt den PHY041-Kernel auf L=24/32/48 mit getrennten Upsilon_2-/Upsilon_4-FSS-Kanaelen und Multi-Walker-Systematik aus (L=24: 1 Walker; L=32/48: je 3). Gate-Log: `results/260707 PHY042 honeycomb wl-fss L24-32-48 gate report.json` (`overall_pass=True`, 9/9 Pipeline-Gates, Lauf 2026-07-06, ~3498 s wall, `master_seed=42`, `lnf_final=1e-5`). Die Pipeline-Gates pruefen **Integritaet, nicht Physik-Wahrheit**; der Physik-Befund ist ein FINDING.

| Gate / Diagnostik | Wert |
|---|---|
| A: aligned `Upsilon_2(0) = 3/4 J` exakt | PASS |
| max. kanonisches Rand-Leak | 1.6e-10 |
| max. unbesetzte kanonische Masse in-Domaene | 1.9e-10 |
| Validitaets-Domaene `T_max` (Walker-Spread < 0.04) | L24: 0.67, L32: 0.60, L48: 0.585 |

Upsilon_2-Paar-Schaetzer (C-eliminiert, Mittel-Kurven), Einordnung **nur relativ zum Referenzband**:

| Paar | T_BKT | belastbar? |
|---|---:|---|
| (24,32) | 0.5875 | ja — Crossing in beiden Validitaets-Domaenen |
| (24,48) | 0.5898 | **nein** — NR-PHY042-02 |
| (32,48) | 0.5916 | **nein** — NR-PHY042-02 |

Ehrliche Grenzen:

- **Kein neuer Bestwert / kein finaler T_BKT-Claim.** L=48 bleibt endlich; die Walker-Systematik misst nur den Sampler (g(E)-Bias), nicht den finite-size-Bias. Ein neuer Bestwert braucht Cross-Family-Review.
- **NR-PHY042-02:** Die Paare (24,48) und (32,48) crossen ausserhalb der gemeinsamen Validitaets-Domaene und sind bei diesem Statistik-Budget nicht belastbar (Walker-Spread 0.033 bzw. 0.039 in T_BKT vs. 0.013 fuer das quotable Paar (24,32)).
- **Erratum 2026-09-26 (Issue #45 §2):** L=24 lief als Einzel-Walker; der Report fuehrt dafuer `domain_tmax_spread004 = 0.67` und `walker_spread = 0` - beides Konstruktions-Artefakte (Gitterende/Null), **keine Messung**. Korrigierte Semantik (`null` + Grund) in `results/260926 PHY042 domain semantics erratum.json`, MC-frei aus den gespeicherten Kurven abgeleitet; alle Quotierbarkeits-Urteile bleiben gleich. Ehrliche Folge: das einzige quotierbare Paar (24,32) hat nur **einseitig** (L=32) gemessene Sampler-Evidenz (Basis `partial`). Der gepinnte Report bleibt byte-unveraendert.
- **NR-PHY042-03:** Der Upsilon_4-Dip ist bei L in {32,48} nicht walker-robust (L=24 robust bei T=0.65; L=32/48 walker-abhaengig). Upsilon_4-FSS oberhalb L=24 braucht mehr Produktion.

## PHY043 — Konventionsfreier Quercheck triangular (Audit O1; Finding, kein Bestwert)

Audit-Punkt O1 fragt, ob die per-Site-Normierung auf triangular
(Flaeche/Site = sqrt(3)/2 < 1) das Helicity-Crossing im Limes nach UNTEN
verschiebt — der interne per-Site-Wert 1.4007 liegt -1.22% unter der
Referenz 1.418, konsistent mit genau so einem Restbias. PHY043 misst
deshalb NORMIERUNGSFREIE Observablen (Binder U4, Korrelations-Ratio
xi_2/L aus dem Strukturfaktor) auf L=9/13/19/25 mit derselben
Wolff-Pipeline (8 Seeds, 800 Messungen, T=1.36..1.70, seed=42; Gate-Log:
`results/260808 PHY043 triangular convention-free crossing report.txt`,
OVERALL PASS 6/6, ~688 s).

Ergebnis (vorab festgelegter Interpretations-Vertrag, Spec §6):

- Alle xi_2/L-Splay-Temperaturen (persistenter 2-sigma-Splay) liegen bei
  T = 1.50..1.60, die U4-Splays bei 1.60..1.64 — durchweg OBERHALB von
  1.41; die Merge-Region (kritische Phase) ist bis dorthin intakt.
- KEIN Splay-Beginn unterhalb von 1.40: die konventionsfreien Signaturen
  stuetzen die Hypothese NICHT, das wahre T_BKT laege beim per-Site-Wert
  1.4007 oder darunter — qualitativ konsistent mit der Referenz-Lage.
- Adjacent-Crossings im Merge-Bereich sind ueberwiegend < 2 sigma
  (Rauschen, ehrlich mit Flanken-Signifikanz berichtet).
- **NR-PHY043-01 (Grenze):** Splay-Temperaturen sind obere Schranken der
  Merge-Region und driften logarithmisch; eine 1%-Diskriminierung
  zwischen 1.4007 und 1.418 leistet dieses Budget NICHT. O1 bleibt formal
  offen; naechster Pfad ist der Konventions-Nachweis je Referenz in
  SOURCES.md. Kein per-Site-Code-Fix ohne diesen Nachweis.

## W4 v02 — honeycomb T_BKT mit Upsilon pro Flaeche (PHY046, vorregistriert; FINDING)

Protokoll: `spec/260926 PHI HEX w4 honeycomb preregistration v02.md`
(committet vor der Datennahme). Kanonischer Wolff (Numba), Leiter 32..256,
16 T-Punkte, 8 Seeds, 2.9 CPU-h. Gate-Log: `results/260926 PHY046 honeycomb
w4 wolff area-helicity report.*`; Einordnung: `results/260926 PHY046 W4 v02
interpretation note.md`.

| Paar (L, 2L) | (32,64) | (48,96) | (64,128) | (96,192) | (128,256) |
|---|---:|---:|---:|---:|---:|
| T*, Upsilon pro Flaeche | 0.5658 | 0.5655 | 0.5657 | 0.5649 | 0.5663 |
| T*, per Site (Quercheck) | 0.5946 | 0.5916 | 0.5903 | 0.5873 | 0.5865 |

- **T_W4 = 0.5649, sigma_tot = 0.0040 -> vorregistriertes Verdikt CONSISTENT**
  mit B = [0.560, 0.580].
  - Flag: der HKS-Fit liegt am b-Rand, weil die Paare flach sind.
  - Die Varianten 0.5646..0.5673 stimmen ueberein.
- **Post-hoc-V&V (PHY048, square, gleiche Leiter):**
  - Die Paar-Crossings liegen auf <= 0.5 % bei T_BKT = 0.893.
  - Der 3-Parameter-HKS-Fit ist rauschverstaerkend (+0.45 %, sigma 0.009).
  - Lehre fuer kuenftige Protokolle: groesstes Paar bzw. Paar-Mittel als
    Primaer-Schaetzer.
- **Offene Spannung:** ~1.3 % (~2 sigma) unter den normierungsfreien
  Literaturwerten (0.572-0.573). Primaertext-Abgleich 2026-09-27: Jiang/PTEP
  verwendet fuer honeycomb explizit den Dichtefaktor 4/(3 sqrt(3)) im
  Helicity-Modulus und vergleicht diesen korrigierten Gamma mit 2T/pi.
  Die fruehere Erklaerung "Literaturwert ist per Site" ist fuer diesen
  0.576(4)-Kanal damit falsifiziert. Die Diskrepanz bleibt offen.
- Kein Bestwert. Vor einer externen Aussage braucht es einen eigenen
  normierungsfreien Schaetzer, den Primaertext-Abgleich und ein
  Cross-Family-Review.

## PHY044 — W4-Kalibrierung honeycomb (blind; Budget, kein T_BKT)

Vertrag: `spec/260926 PHI HEX w4 honeycomb preregistration v01.md` (separat
und **vor** den Laeufen committet). Gate-Log: `results/260926 PHY044 honeycomb
wl calibration report.{json,txt}` (OVERALL PASS 7/7, seed=42, BLAS auf 1 Thread
gepinnt, ~47 min auf 4 Kernen). **Blind:** der Report enthaelt keinen
T_BKT-Lagewert (maschinell geprueft), nur Kosten, Domaenen und Streuungen.

**Kernel.** `wl_entropic_fast` ist ein Numba-Zwilling des PHY041-Kernels mit
identischem RNG-Verbrauch und identischer Arithmetik je Update:

- **bit-identisch** zum Original auf derselben Maschine (kleine Gitter in CI;
  identischer L=24-Produktionsjob, Speedup 5.6x);
- reproduziert alle 7 committeten PHY042-Walker mit **exakt gleicher
  Trajektorie** (wl_sweeps, Bin-Belegung) - Kurven bis rtol 1e-12 (ULP-Ebene
  der numpy-Reduktionen, Plattform); PHY042 neu in 231 s (Originallauf
  ~3500 s auf ANDERER Hardware - kein Speedup-Mass; gleiche Maschine: 5.6x).

| K1 Kosten je Walker (PHY042-Rezept) | L=24 | L=32 | L=48 | L=64 |
|---|---:|---:|---:|---:|
| Spin-Updates | 1.5e8 | 3.0e8 | 1.3e9 | 3.6e9 |
| Wall Numba [s] | 25 | 48 | 215 | 600 |
| T_max der Walker-Spread-Domaene (3 Walker) | 0.575 | 0.630 | 0.610 | 0.615 |
| 1/t-Phase erreicht | 3/3 | 3/3 | **0/3** | **0/3** |

Kosten ~L^3.3 (lokal 3.6). **K3 Hebel bei L=48** (gepaart, gleiche g(E)):

| Variante | max. Walker-Spread | T_max | Wall/Walker |
|---|---:|---:|---:|
| Basis (prod 1x, lnf 1e-5) | 0.268 | 0.610 | 215 s |
| **prod 4x** | **0.074** | **0.620 = T_req** | 478 s |
| lnf 1e-6 (1/t greift) | 0.206 | 0.600 | 770 s |
| Zerlegung: 1 g(E), 3 Produktionen | 0.203 | 0.605 | - |

Befunde (ehrlich):

- **Der Walker-Spread ist produktions-dominiert** (Zerlegung: 94 % des Spreads
  bis T_req bei FESTER g(E)); laengere WL/1-t-Politur hilft nicht, mehr
  Produktion schon. Budget-Regel der Spec -> **W4-GO mit prod 4x** - bei L=48
  **ohne Marge** (T_max = T_req); ab L=64 muss der W4-Lauf T_req je L selbst
  messen (Spec §4.4), sonst Stop.
- **NR-PHY044-01:** mit lnf_final = 1e-5 erreicht WL bei L>=48 die 1/t-Phase
  nicht (endet im Halbierungs-Regime); betrifft auch 2 der 3 committeten
  PHY042-L=48-Walker. Nach K3 ist das nicht der limitierende Fehler.
- **NR-PHY044-02:** erstmals mit 3 Walkern gemessen liegt T_max(L=24) beim
  PHY042-Rezept bei 0.575 - **unter** dem committeten PHY042-Crossing
  (24,32) = 0.5875. Das einzige als QUOTIERBAR gefuehrte PHY042-Paar haelt der
  W4-Regel "beidseitig gemessene Domaene" nicht stand (Evidenz bleibt als
  Finding mit Basis `partial` stehen).
- Walker-Streuung der Paar-Crossings beim Basis-Rezept 0.013..0.025
  (-> sigma_sampler ~0.006..0.012, nahe an der Stop-Schwelle sigma_tot 0.010).
- Projektion W4-Rezept (Numba, je Walker): L=64 0.4 h, L=96 1.8 h, L=128 5.2 h
  (Extrapolation; Python-Kernel ~5.6x laenger).
- Reproduzierbarkeit: np.dot haengt ab n=12288 (L=64) bitweise von der
  BLAS-Thread-Zahl ab (verifiziert) -> PHY044 pinnt BLAS auf 1 Thread.

## Methodische Kernformeln

```text
H = -J * sum_<ij> cos(theta_i - theta_j)
T = 1 / beta
sigma_T = sigma_beta / beta^2
Upsilon(T_BKT) = 2 T_BKT / pi
Upsilon(T_BKT,L) = (2 T_BKT / pi) * (1 + 1 / (2 ln L + C))
A(E -> E') = min(1, g(E) / g(E'))
accept if log(u) <= log_g(E) - log_g(E')
w_T(E) = g(E) * exp(-E/T)
leak(T) = sum_edge w_T(E) / sum_all w_T(E)
```

## Reproduzieren

```bash
pytest -m "not slow"
pytest
python "src/260607 PHY032 honeycomb wm-logfit bootstrap v01.py"
python "src/260616 PHY040 wang-landau entropic helicity v01.py"
python "src/260702 PHY041 honeycomb wang-landau entropic helicity v01.py"
python "src/260706 PHY042 honeycomb wl fss v01.py"
python "src/260808 PHY043 triangular convention-free crossing v01.py"
python "src/260706 PHY042 honeycomb wl fss v01.py" --reanalyse   # Domaenen-Erratum
python "src/260926 PHY045 helicity normalization O1 test v01.py" all > ab.json   # O1-Nachweis
python "src/260926 PHY046 honeycomb w4 wolff area-helicity v01.py"               # W4 v02 (~45 min)
python "src/260926 PHY047 per-area reanalysis committed v01.py"                  # Neuauswertung
python "src/260926 PHY048 square pipeline validation v01.py"                     # post-hoc V&V
# PHY044 (numba empfohlen; Stufen je ~4..21 min auf 4 Kernen)
M="src/260926 PHY044 honeycomb wl calibration v01.py"
for st in ladder levers valbit speed; do python "$M" $st --out stages/; done
python "$M" report --out stages/
```

CI prueft Lint, Compile und schnelle Korrektheits-Gates. Slow-Messlaeufe bleiben lokal.

## Struktur

```text
src/        Engines + Experiment-Code
spec/       Specs, Theorie, Audits
results/    Gate-Logs, Reports, Negativ-Results
tests/      schnelle Gates + slow Mess-Smokes
archive/    Vorgaenger-Versionen
SOURCES.md  Provenance / SHA-256
```

## Naechste Stufe nach Selbstpruefung / W4 v02 (2026-09-26b)

1. **W4-v03a / PHY049+PHY050:** normierungsfreier Correlation-Ratio-Quercheck
   ist vorregistriert; der deterministische FSS-Estimator und G4-Synthetic-
   Recovery inklusive adversarial Nullkontrollen sind implementiert. G4
   autorisiert nur die Produktionsmessung: NO_PHYSICS_INTERPRETATION bleibt
   bis G0/G5/G6 auf echten PHY049-Daten und committetem Gate-Log bestehen.
2. **Primaertext-Abgleich:** Jiang/PTEP und de Andrade v4 sind fuer
   Honeycomb-Geometriefaktor und aktuelle direkte T-Kanaele geschlossen.
   Historische de-Andrade-beta-Werte bleiben versionsspezifische Lineage
   (v2 primaer verifiziert; 1.687/1.635 ohne aktuellen Referenz-Claim).
3. **triangular und kagome pro Flaeche neu messen** (Wolff-Numba, Leiter bis 256, T-Gitter um 1.465 bzw. 0.82); triangular gegen die normierungsfreie Reihe 1.465.
4. Protokoll-Lehre aus PHY048: in kuenftigen Vorregistrierungen das groesste Paar bzw. das Paar-Mittel als Primaer-Schaetzer, HKS-3-Parameter nur als Variante.

## Naechste Stufe nach Issue #45 / PHY044 (2026-09-26, teilweise ueberholt)

1. **W4-Produktionslauf** nach Vorregistrierung: Leiter {32, 48, 64, 96} (L=128 zulaessig, Projektion 5.2 h/Walker <= 24 h), je 3 Walker, prod 4x, Numba-Kernel mit BLAS-Pin; Stop-Regeln S1/S2 und T_req-Messung je L (Spec §4.4). Ein Budget-Nachtrag (z. B. prod 8x fuer L>=64) nur VOR Sicht auf W4-Daten und als v02 der Spec.
2. **Primaertext-Abgleich** der `search_corroborated`-Referenzkanaele (arXiv:2406.12076 v1..v4, PTEP-Abstract von 2406.14812) - in einer Umgebung mit arXiv-Zugang.
3. Walker-robuster Upsilon_4-Dip fuer L>=32 (NR-PHY042-03) im W4-Lauf als Sekundaer-Kanal mitfuehren.
4. Kein T_BKT-Bestwert-Claim ohne Cross-Family-Review.
5. Square-V&V-Re-Baseline: mehr Statistik/groessere L, damit das 1%-Gate wieder belastbar bindet; PHY039-Dip-Fenster erweitern (Audit O-Punkte).
6. O1 abschliessen: der MC-Quercheck-Teil ist mit PHY043 gelaufen (qualitativ konsistent mit Referenz-Lage, NR-PHY043-01: keine 1%-Diskriminierung); offen bleibt der Konventions-Nachweis je zitierter Referenz (per-Spin vs per-Flaeche) in SOURCES.md. Bis dahin kein per-Site-Code-Fix.
7. Re-Run PHY030v02/PHY032/PHY033 mit Pol-Guard bei der naechsten Produktions-Runde (Audit-Punkt O4).
8. Bei Neu-Produktion PHY040: L in den WL-Stream aufnehmen (O8) und PHY039 auf ungebuchte Stream-Basis 800+ heben (O3).

---
*Coworker Research / Coworkerz | Repo-Anlage 2026-06-04 nach arbeitsschablone_forschungs-repo-anlage*
