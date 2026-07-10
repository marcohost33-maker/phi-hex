# 260710 PHI HEX Code-Audit v01 — Befunde, Fixes, offene Punkte

Selbst-Audit / Coworkerz, 2026-07-10. Drei unabhaengige Review-Passes
(Kern+PHY026–033, PHY034–042/Wang-Landau, Tests/Infra/Doku); jeder Befund
vor dem Fix am Code bzw. numerisch verifiziert. Dieses Dokument ist die
Vertragsquelle fuer die Audit-Fixes; Gate-Evidenz in `results/` und
`tests/` (siehe CHANGELOG-Eintrag 2026-07-10).

## 1. Bestaetigte und GEFIXTE Defekte

### H1 — PHY028: Wrap-Bond-Vorzeichen im Quadratgitter (P0)

`build_square_lattice` sortierte Kanten mit `(min,max)`, liess `edge_disp`
aber bei festem `(+1,0)/(0,+1)`. Fuer die 2L Wrap-Bonds (x=L-1->0,
y=L-1->0) war der gespeicherte Bond-Vektor damit das NEGATIV des
min-image-Vektors des gespeicherten Paars; `sin(theta_i - theta_j)*proj`
ging mit falschem Vorzeichen ein. `T1` (cos, gerade) blieb korrekt — das
T=0-Gate KONNTE den Fehler nicht sehen.

Verifikation: uniformer Twist `theta=(2 pi/L) x` bei L=8: `sin_accum`
-33.94 statt exakt -45.25 (16/128 Bonds invertiert). Wirkung
thermodynamisch: `<sin_accum^2>` verfaelscht, Upsilon ~25% zu TIEF
(L=16, T=0.90: 0.486 buggy vs 0.646 fixed, seed=42) — der alte committete
PHY028-Report entspricht byte-fuer-Wert dem Bug-Code.

**Konsequenz (ehrlich):** Der alte V&V-Anker "T_BKT(16,32)=0.8950, +0.16%"
war ein ZUFALLSTREFFER auf systematisch falschen Upsilon-Werten. Nach dem
Fix (deterministischer Re-Run seed=42): (8,16)=0.8867 (-0.76%),
(16,32)=0.8841 (-1.05%) — das 1%-Gate ist knapp FAIL. Das Gate wurde NICHT
aufgeweicht; der V&V-Anker ist auf "~1% bei L<=32, 12 Seeds" herabgestuft
(Re-Baseline mit mehr Statistik/groesseren L = naechste Stufe). Betroffen
und neu gerechnet: PHY028, PHY039, PHY040 (beide importieren
`build_square_lattice` aus PHY028). PHY039-Re-Run: Y4-Dip ausserhalb des
T-Fensters (Dip-Gates FAIL, Dip-Lage 0.9316 zurueckgezogen). PHY040-Re-Run:
WL-vs-Wolff-Gates PASS (intern konsistent), Paare (12,16)=0.9005 (+0.85%),
(12,24)=0.9117 (+2.10%), (16,24) ohne Crossing — Bias-Gate ehrlich FAIL;
Instrumentierung bestaetigt "1/t-Phase=NIE" live. NICHT betroffen:
Kern-Torus (bond-weise verifiziert 0/243 inkonsistent), honeycomb/
kagome-Builder (natuerliche Orientierung), PHY041/042.

Fix: Kanten in natuerlicher Orientierung speichern (kein min/max).
Gate: `tests/test_phy028_square_geometry.py` (min-image-Konsistenz aller
Bonds + exaktes Uniform-Twist-Orakel).

### M1 — Pol-Artefakt in allen C-freien Paar-Root-Findern (P1)

`1/R(T,L1) - 1/R(T,L2) - target` wurde auf Vorzeichenwechsel gescannt,
aber R crosst am finite-L-NK-Crossing selbst null: der Pol von 1/R wurde
als "Nulldurchgang" interpoliert (verifiziert: synthetisches Szenario ohne
WM-Wurzel lieferte T_BKT=1.449). In PHY032/033 liefen solche Artefakte als
"valide" Bootstrap-Samples ins CI. Auf dem physikalischen WM-Ast gilt
`1/R = 2 ln L + C > 0` — Punkte mit `R <= 0` werden jetzt uebersprungen
(PHY028, PHY029, PHY030 v01+v02, PHY031, PHY040; via Import wirksam in
PHY032/033/034/035/041/042). Gates: Pol-Rejektions-Orakel in
`tests/test_phy028_square_geometry.py`.

**Nachschaerfung 2026-07-10b (Codex-Review PR#23, verifiziert):** das
blosse Ueberspringen der R<=0-Punkte interpolierte noch QUER ueber die
Pol-Luecke (Vorzeichenwechsel zwischen nicht-benachbarten Gitterpunkten =
disconnected branch; Repro: Luecken-Szenario lieferte 1.4311). Jetzt wird
ein Crossing nur zwischen BENACHBARTEN Gitterpunkten akzeptiert (alle 6
Root-Finder). Committete v02-Evidenz unveraendert (PHY028-Paare und
PHY040-Paare aus den Report-JSONs bit-genau reproduziert); Gap-Brueckungs-
Orakel als CI-Gates ergaenzt.

Hinweis Evidenz: die committeten PHY030v02/PHY032/PHY033-Reports sind vor
dem Guard entstanden; ihre Scan-Fenster lagen ueberwiegend im R>0-Bereich,
Bootstrap-Tail-Resamples koennen aber betroffen sein. Re-Run bei der
naechsten Produktions-Runde (offener Punkt O4).

### WL-H1 — PHY040: 1/t-Phase war toter Code (P1, Method-Drift)

Die Umschaltbedingung `lnf <= 1/sweeps` (nur am Halbierungs-Event geprueft)
war mit den committeten Parametern (lnf_final=2e-4, flat=0.7) nie wahr:
die Laeufe endeten in reinem Standard-WL mit bekannter Fehler-Saettigung
(~0.02 in <E>/N; PHY041 hat das unabhaengig gefunden und dort die echte
B&P-Phase implementiert). Der Method-String claimte trotzdem "WL + 1/t".
Fix: Instrumentierung `WLResult.one_over_t_engaged` + korrigierter
Method-String + Header-Nachtrag; Sampler-Verhalten UNVERAENDERT (bewusst:
PHY041 ist der korrigierte Kernel, PHY040 bleibt die dokumentierte Stufe).

### WL-H2 — PHY041: Coverage-Massen-Luecke (P1)

`_canonical_weights` maskiert unbesetzte Bins und renormiert still; das
Rand-Leak-Gate deckte nur die aeussersten 3 Bins je Fensterrand. Kanonische
Masse auf unbesetzten Bins WEITER INNEN (committeter Lauf: 491/512 Bins
besetzt bei L=24) war ungegatet. Fix: `uncovered_canonical_mass` +
`PASS_NO_UNCOVERED_CANONICAL_MASS` (Backport des PHY042-Gates; PHY042 hat
fuer denselben L=24-Kernel 1.9e-10 gemessen — benigne, aber jetzt gegatet).

### Weitere gefixte Punkte (P2)

- **M3 (PHY034–038):** `None`-Paar-Schaetzer crashten den Lauf VOR dem
  Gate-Report (TypeError statt FAIL-Evidenz) — None-sichere Formatierung/
  Extrapolation; FAIL-Pfad schreibt jetzt Report + Exit 1.
- **WL-M1 (PHY041/042):** PHY032-Drift-Guard passte vacuous, wenn der
  Report-Parser nichts fand — Parse-Vertrag (8 T-Punkte je L) ergaenzt.
- **WL-M2 (PHY042):** Wolff-Y2-Gate passte vacuous, wenn KEIN
  Validierungs-T in der Validitaets-Domaene lag — `y2_any_in_domain`.
- **WL-L2 (PHY042):** `max()` ueber leeren Generator bei n_walkers=1
  crashte vor dem Report — `default` ergaenzt.
- **L3 (PHY026):** `PASS_CLUSTER_NONTRIVIAL` war per Konstruktion immer
  wahr — jetzt 0.05 < frac < 0.95 (committed: 0.38..0.50).
- **L4 (PHY030v02):** Residuen mischten zwei Temperaturen (Daten@T_min vs
  Modell@T_bkt) — jetzt beide bei T_min; plus `C_at_bound`-Diagnose (L5).
- **L6 (PHY028/031):** fails-open-Gates (fehlendes Paar liess OVERALL
  passen; "LARGEST_PAIR" war das letzte konvergierte) — fails-closed und
  explizit gebunden.
- **L4-WL (PHY040/041):** Anneal-/WL-Schleifen ohne Abbruchkriterium
  konnten endlos haengen — Caps mit RuntimeError.
- **Vacuous-/schwache Tests:** Referenzkonventions-Gate testete nur
  test-lokale Literale (jetzt Produktions-Konstanten); Edge-Rejektions-Gate
  akzeptierte fast alles (jetzt garantierter Rand-argmin -> None);
  REF_BAND-Substring-Check gehaertet.
- **Coverage:** phy017, PHY027, PHY029, phy024 R0/R2 hatten NULL Tests
  (nur compile+lint) — Import- + Orakel-Smokes ergaenzt; 6 Sub-4s-Tests
  faelschlich als `slow` markiert (CI lief ohne die staerksten Gates,
  inkl. RNG-Stream-Bit-Vertrag) — de-slowed.
- **SOURCES-Gate:** Matching war pfad-ungebunden (git mv/Kopie entkam),
  CRLF-Waesche galt fuer alle Zeilen, Root-Level-.py entkam allen Gates —
  pfadgebunden + Windows-only-Varianten + Scope-Escape-Gate; `docs/` in
  den Scope aufgenommen.
- **Infra:** Dependabot ohne pip-Ecosystem trotz requirements-dev.txt;
  ruff lintete tests/ nicht; CI-Matrix ohne 3.13; ADR-001 widersprach
  AGENTS.md (dde) — alles nachgezogen.

## 2. OFFENE Punkte (dokumentiert, bewusst NICHT still gefixt)

### O1 — Konventionsfrage: per-Site-Upsilon vs NK-Gerade auf nicht-quadratischen Gittern

Die universelle Sprung-Bedingung `Upsilon(T_BKT) = 2 T_BKT/pi` ist fuer die
Steifigkeit der Kontinuums-Gradienten-Wirkung (per FLAECHE) hergeleitet.
Mit Bond-Laenge a=1 gilt per-Site = per-Flaeche * (Flaeche/Site):
sqrt(3)/2 ~ 0.866 (triangular), 3*sqrt(3)/4 ~ 1.299 (honeycomb),
2*sqrt(3)/3 ~ 1.155 (kagome), exakt 1 (square — deshalb kann PHY028 die
Konventionen NICHT diskriminieren). Fuer Faktoren > 1 konvergiert das
Crossing trotzdem gegen T_BKT (Upsilon springt oberhalb auf 0); fuer
triangular (0.866 < 1) laege das per-Site-Crossing im Limes strikt
UNTERHALB des universellen Sprungs. Gegenevidenz: mehrere zitierte
Referenzen definieren Upsilon per Spin, und die per-Site-Resultate treffen
unabhaengige Referenzen auf drei Gittern (triangular 1.4007 vs 1.418 —
konsistent mit einem kleinen Restbias in genau dieser Richtung).
**Naechster Schritt:** konventionsfreier Quercheck (Binder/Korrelations-
Ratio) auf triangular ODER Konventions-Nachweis je Referenz in SOURCES.md.
KEIN Code-Fix ohne diesen Nachweis (Reality-Anchor).

### O2 — RNG: T-Scans teilen den Stream je (L, Seed)

Alle Wolff-Module starten je (L, seed) fuer JEDES T dieselbe RNG-Sequenz
und Startkonfiguration. Das ist als Common-Random-Numbers-Technik fuer
Crossings varianzREDUZIEREND und bit-reproduzierbar committed — aber die
T-Punkte eines Scans sind korreliert (SEM je T bleibt korrekt; die
Interpolation ueber T-Punkte behandelt sie implizit als unabhaengig).
Dokumentiert als bewusste Design-Entscheidung; bei einer kuenftigen
Neu-Produktion T-Index in den Stream aufnehmen (invalidiert alle
Bit-Repro-Vertraege — nur mit vollem Re-Run aller Reports).

### O3 — Stream-Basis 300+ doppelt gebucht (PHY028 + PHY039)

Gleiche Basis, gleiche Gitterfamilie, ueberlappende L; heute trennen NUR
die disjunkten T-Gitter die Ketten. In PHY039/PHY041/PHY042 als Kommentar
inventarisiert; bei T-Gitter-Aenderung ZUERST PHY039 auf ungebuchte Basis
(z.B. 800+) heben.

### O4 — Re-Runs nach M1-Guard

PHY030v02/PHY032/PHY033(+PHY034/035-Ketten) bei naechster Produktions-
Runde mit Pol-Guard neu rechnen und Reports neu pinnen. Teil-Entwarnung
(verifiziert 2026-07-10): die deterministische Re-Analyse der committeten
PHY035-Vollparameter (8/140/60, seed=42) MIT Guard reproduziert die
committeten Punktschaetzer EXAKT ((8,16)=0.5969, (12,24)=0.6014,
(16,32)=0.5993, Extrapolation 0.6052) — die Punkt-Evidenz war nicht
pol-beeinflusst. Offen bleiben die Bootstrap-Tail-Resamples (CIs).
Gegenteiliges Mini-Statistik-Beispiel als Beleg der Guard-Wirkung: der
PHY035-Smoke mit 4/50/25 bezog das (16,32)-"Crossing" aus der Pol-Region
(Smoke-Statistik deshalb auf 6/100/40 angehoben, Test-Kommentar).

### O5 — PHY025-Generator nie migriert

`results/260602 PHY025 fss report v02.txt` + Spec existieren; der
erzeugende Code liegt nur in der Drive-Lineage (SOURCES-Zeile), nicht im
Repo. Bewusst offen (Report ist Evidenz, Methode in PHY027+ ersetzt).

### O6 — PHY034-CI approximativ

Die MC-Fehlerfortpflanzung behandelt (12,24)/(24,48) als unabhaengig,
obwohl beide die L=24-Daten teilen (jetzt im Code/Report ausgewiesen);
die korrekte gemeinsame Seed-Bootstrap-Propagation macht PHY035.

### O7 — Kern-Diagnostik-Pfad (LOW, nicht messungsrelevant)

Im Langevin-/Defekt-Diagnostik-Pfad des Kerns (NICHT in der Wolff-/WL-
Messpipeline): (a) `phase_correlation` und `DefectTrackerV2` nutzen auf dem
Torus rohe (ungewrappte) Positionen — g1-Bins und Trajektorien-Matching
ueber die Naht sind verzerrt; (b) die open-boundary-`cells` sind
Wheel-Loops (Zentrum+6 Nachbarn), ein Vortex kann als bis zu 3 gleich
geladene Defekte gezaehlt werden (Torus-Pfad nutzt korrekte Elementar-
Dreiecke); (c) `integrated_autocorr_time` normiert die FFT-Autokovarianz
nicht per Lag ((n-w)/n-Bias, bei w << n vernachlaessigbar). Helicity ist
in allen Faellen unberuehrt (nutzt edge_disp). Fix nur zusammen mit einem
Re-Run des Kern-Selftest-Reports.

## 3. Referenz der Review-Verdikte (sauber befundene Kerne)

Wolff-Update (Bond-Wahrscheinlichkeit, Spiegelung, Frontier, detailed
balance), Kern-Helicity-Maschinerie (Vektorisierung, Cache-Invalidierung,
Torus-edge_disp), honeycomb/kagome-Geometrie inkl. T=0-Orakel, WL-Akzeptanz
und Bin-Randbehandlung in allen drei Kerneln, log-sum-exp-Rueckgewichtung,
PHY039-F2/F4-Kumulanten (gegen unabhaengige Herleitung + numerisches
Orakel), PHY037-Parser/Drift-Guard (fails-closed), PHY042-Walker-Streams
(kollisionsfrei) — alle drei Passes ohne Defektbefund.
