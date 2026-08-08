# 260808 PHI HEX PHY043 — konventionsfreier Quercheck triangular (Methode) v01

Coworker Research / Coworkerz, 2026-08-08.
Vertragsquelle fuer `src/260808 PHY043 triangular convention-free crossing v01.py`.

## 1. Zweck und Einordnung (Audit-Punkt O1)

Das Code-Audit 2026-07-10 (`spec/260710 PHI HEX code audit v01.md`, §2 O1)
laesst offen, ob die per-Site-Normierung des Helicity-Modulus auf
nicht-quadratischen Gittern die universelle Sprung-Bedingung
`Upsilon(T_BKT) = 2 T_BKT / pi` verzerrt: per-Site = per-Flaeche *
(Flaeche/Site), fuer triangular ist der Faktor sqrt(3)/2 ~ 0.866 < 1,
das per-Site-Crossing laege im Limes strikt UNTERHALB des universellen
Sprungs. Der interne per-Site-Bestwert (PHY030 v02, WM-Fit) ist
1.4007 +/- 0.0081 gegen die Referenz 1.418 (arXiv:2501.07388) — eine
Abweichung von -1.22 %, konsistent mit genau so einem Restbias.

O1 fordert als naechsten Schritt einen KONVENTIONSFREIEN Quercheck
(Binder/Korrelations-Ratio) auf triangular, BEVOR ein Code-Fix oder ein
neuer triangular-Claim erfolgt. PHY043 ist dieser Quercheck.

**Kein T_BKT-Bestwert-Claim.** PHY043 liefert ein FINDING zur Lage der
konventionsfreien Signaturen relativ zu (a) dem internen per-Site-Wert
1.4007 und (b) der Referenz 1.418. Die Pipeline-Gates pruefen
Integritaet, nicht Physik-Wahrheit (gleiches Prinzip wie PHY042).

## 2. Observablen (dimensionslos, normierungsfrei)

Beide Observablen sind Verhaeltnisse und damit unabhaengig von jeder
per-Site-/per-Flaeche-Konvention:

1. **Binder-Ratio** `U4(T, L) = <|m|^4> / <|m|^2>^2` mit
   `m = (1/N) sum_j exp(i theta_j)`.
   Grenzwerte: geordnet -> 1; ungeordnet (iid-Limes) -> exakt `2 - 1/N`
   (komplexer Gauss-Limes; exaktes Orakel fuer Tests, Herleitung:
   `<|m|^2> = 1/N`, `<|m|^4> = (2N-1)/N^3` fuer iid-uniforme Winkel).

2. **Korrelations-Ratio** `xi_2/L` aus dem Strukturfaktor
   `S(k) = (1/N) |sum_j exp(i theta_j) exp(-i k . r_j)|^2`:

   `xi_2 / L = sqrt(S(0)/S(k1) - 1) * sqrt(3) / (4 pi)`.

   Begruendung Vorfaktor: zweite-Moment-Definition
   `xi_2 = (1/|k1|) sqrt(S(0)/S(k1) - 1)` mit `|k1| = 4 pi / (sqrt(3) L)`
   auf dem triangular-Torus (Bond-Laenge a=1). Der Vorfaktor ist
   L-unabhaengig und kuerzt sich aus jeder Crossing-/Splay-Analyse —
   nur die Kurvengeometrie zaehlt, nicht die absolute Skala.

**Erlaubte Wellenvektoren / Modenbuchhaltung.** Knoten liegen auf axialen
Koordinaten `(q, r)` in `0..L-1` (Rhombus-Torus, `_build_triangular_torus`).
Die exakt torus-periodischen ebenen Wellen sind
`exp(2 pi i (n1 q + n2 r) / L)`, n1, n2 ganz. Mit Primitivvektoren im
60-Grad-Winkel gilt `|k(n1,n2)|^2 propto n1^2 + n2^2 - n1 n2`; die drei
symmetrieaequivalenten Minimalmoden sind `(1,0)`, `(0,1)`, `(1,1)`
(alle Norm 1; `(1,-1)` hat Norm 3). `S(k1)` ist das Mittel ueber die drei
Minimalmoden.

**Bewusst NICHT verwendet:** der universelle Wert `(xi_2/L)* ~ 0.7507`
(Hasenbusch, J. Phys. A 38 (2005) 5869) gilt fuer den QUADRATISCHEN Torus
(Modularparameter tau = i). Unser LxL-Rhombus-Torus hat tau = exp(i pi/3);
der universelle Wert ist geometrieabhaengig und darf hier nicht als Anker
benutzt werden. PHY043 nutzt deshalb ausschliesslich Crossing-/Splay-
Signaturen zwischen L-Paaren derselben Geometrie.

## 3. BKT-Signaturen: Splay statt einfacher Crossing-Annahme

Unterhalb T_BKT ist die ganze Phase kritisch: `xi_2/L`- und `U4`-Kurven
verschiedener L MERGEN (bis auf Korrekturen), oberhalb splayen sie
(groesseres L: kleineres `xi_2/L`, groesseres `U4`). Ein sauberer
Einzel-Crossing-Punkt wie bei gewoehnlichen kritischen Punkten existiert
nicht zwingend; endliche-L-Crossings existieren, driften aber
logarithmisch. Deshalb zwei Auswertungen je Paar (L1 < L2) und
Observable:

1. **Adjacent-Crossings** `T_x(L1,L2)`: Nulldurchgaenge von
   `D(T) = O(T,L2) - O(T,L1)` NUR zwischen benachbarten Gitterpunkten
   (gleicher Vertrag wie der M1/P2-gehaertete Root-Finder-Stack; keine
   Interpolation ueber Luecken). Je Crossing wird die Signifikanz
   `max(|D|/sigma_D)` der flankierenden Punkte mitberichtet — Crossings
   im Merge-Bereich koennen Rauschen sein.
2. **Splay-Temperatur** `T_splay(L1,L2)`: kleinstes Gitter-T, ab dem
   `D(T)` fuer ALLE folgenden Gitterpunkte das erwartete Splay-Vorzeichen
   traegt UND `|D| > z * sigma_D` (z = 2) erfuellt
   (`sigma_D = sqrt(sem1^2 + sem2^2)`, Seed-SEMs). Persistenz-Vertrag:
   ein einzelner 2-sigma-Ausreisser genuegt nicht.

Unsicherheiten: parametrischer Bootstrap — jeder Kurvenpunkt wird als
`Normal(mean, sem)` resampelt (n_boot = 300, deterministisch
`make_rng(master_seed, stream=99043)`), Splay/Crossings je Resample neu
bestimmt, Perzentil-CI [2.5 %, 97.5 %] ueber die Nicht-None-Resamples
plus None-Quote (fails-closed ausgewiesen).

## 4. Sampling-Design

- Gitter: triangular-Torus, radii (4, 6, 9, 12) -> L in {9, 13, 19, 25} —
  die PHY030-Leiter (9/13/19, damit der Quercheck dieselben Groessen
  betrifft wie der per-Site-Befund) plus L=25 fuer ein staerkstes Paar.
- T-Gitter: 1.36 .. 1.70, Delta T = 0.02 (18 Punkte). Fenster bewusst
  DEUTLICH breiter als PHY030: die Splay-Signaturen liegen oberhalb von
  T_BKT, und nahe des Uebergangs ist xi(T) exponentiell gross — das
  Fenster muss bis in den Bereich xi < L_min reichen, sonst mergen alle
  Kurven im ganzen Fenster.
- Wolff: n_measure = 800, n_burn = 300, n_seeds = 8, master_seed = 42.
- **Budget-Herleitung (Pilot, dokumentiert):** ein Pilot mit exakt dem
  PHY030-v02-Budget (radii 4/6/9, T <= 1.58, n_measure=400, n_seeds=4,
  Laufzeit 76 s) lieferte Seed-SEMs in xi_2/L (~0.01-0.02) in derselben
  Groessenordnung wie die L-Differenzen im ganzen Fenster — kein
  persistenter 2-sigma-Splay, Bootstrap-none_frac 0.93-1.00. Das ist
  keine Ergebnis-Selektion, sondern eine Power-Korrektur VOR dem
  Evidenz-Pin: Statistik x2 Seeds, x2 Messungen, Fenster +0.12,
  ein zusaetzliches L. Der Pilot ist damit superseded; committete
  Evidenz ist ausschliesslich der finale Lauf.
- **RNG-Vertrag (Audit O2/O3 umgesetzt):** PHY043 ist Neu-Produktion ohne
  Bit-Repro-Altlast, daher geht der T-Index in den Stream ein:
  `stream = 900 + s + 1000*L + 100000*t_idx`. T-Punkte eines Scans sind
  damit UNABHAENGIG (kein Common-Random-Numbers-Sharing; die Kurvenpunkte
  duerfen in der Interpolation als unabhaengig behandelt werden). Basis
  900 ist im Stream-Inventar frei (100 PHY026, 300 PHY028/039, 400
  PHY031/032, 500 PHY033/040, 600/660 PHY041, 700 PHY040, 99001/99043
  Bootstrap); 800+ bleibt fuer die PHY039-Umbuchung (O3) reserviert.

## 5. Gates (Pipeline-Integritaet, fails-closed)

| Gate | Bedingung |
|---|---|
| PASS_INPUT_COMPLETE | alle (T, L)-Punkte vorhanden, finite, S(k1) > 0 |
| PASS_U4_RANGE | alle U4 in [0.9, 2.1] |
| PASS_HIGH_T_ORDERING_XI | bei T_max: xi_2/L strikt fallend in L |
| PASS_HIGH_T_ORDERING_U4 | bei T_max: U4 strikt steigend in L |
| PASS_LOW_T_MERGE_XI | bei T_min: max. Paar-Abweichung xi_2/L < max(8 % rel., 3 sigma) |
| PASS_SPLAY_XI_LARGEST_PAIR | T_splay(19,25) fuer xi_2/L existiert im Fenster |

Der Vergleich der Splay-/Crossing-Lagen mit 1.4007 (intern) und 1.418
(Referenz) ist DIAGNOSTIK/FINDING, kein Gate (Physik-Wahrheit wird nicht
gegatet; Reality-Anchor: der Report in `results/` ist die Evidenz).

## 6. Interpretations-Vertrag (vorab festgelegt)

- Liegen die xi_2/L-Splay-Temperaturen der groesseren Paare oberhalb von
  ~1.41 und ist der Merge-Bereich bis dorthin intakt, ist das QUALITATIV
  konsistent mit der Referenz-Lage T_BKT ~ 1.418 und liefert KEINE
  Stuetzung fuer die Hypothese, das wahre T_BKT laege beim per-Site-Wert
  1.4007 oder darunter. Umgekehrt waere ein belastbarer Splay-Beginn
  deutlich UNTER 1.40 ein Signal gegen die Referenz-Lage.
- Bei diesem Budget (4 Seeds, L <= 19, Log-Drift der Crossings) ist
  KEINE Diskriminierung auf 1-%-Niveau zu erwarten. Ist das Ergebnis
  nicht diskriminierend, wird das als Negativ-Result **NR-PHY043-01**
  ausgewiesen (Buerger erster Klasse) und O1 bleibt offen — dann ist der
  zweite O1-Pfad (Konventions-Nachweis je Referenz in SOURCES.md) der
  naechste Schritt.
- In JEDEM Fall gilt weiter: kein per-Site-Code-Fix ohne O1-Nachweis.

## 7. Evidenz

Deterministischer Gate-Report nach
`results/260808 PHY043 triangular convention-free crossing report.txt`
(+ JSON-Block auf stdout). Schnelle Orakel-Gates (iid-U4, Spin-Wellen-
Modenbuchhaltung, Splay-/Crossing-Vertraege, Mini-MC-Smoke) in
`tests/test_phy043_convention_crossing.py`.
