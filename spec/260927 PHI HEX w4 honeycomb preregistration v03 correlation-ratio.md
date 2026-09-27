# W4-Vorregistrierung v03: honeycomb T_BKT - normierungsfreier Correlation-Ratio-Quercheck

> **Status:** VORREGISTRIERT am 2026-09-27, vor jeder PHY049-Datennahme.
> Diese Spezifikation aendert oder ersetzt W4-v02/PHY046 nicht. Sie definiert
> einen unabhaengigen, helicity-normierungsfreien Diskriminationstest fuer die
> verbleibende Spannung zwischen PHY046 (~0.565) und der externen
> Correlation-Ratio-Lage (~0.573).
>
> **Primaerquellen-Anker:** Okabe/Otsuka, arXiv:2501.07388/J. Phys. A 58
> 065003 (2025), definieren
> R(T)=<g(L/2)>/<g(L/4)> und nutzen L={48,72,96,144,192};
> fuer honeycomb berichten sie einen groben T_BKT-Wert 0.573.
> Jiang, PTEP 2024 103A02, definiert den Honeycomb-Helicity-Modulus dagegen
> mit dem expliziten Geometriefaktor 4/(3 sqrt(3)); dessen 0.576(4)-Kanal ist
> daher NICHT als per-Site-Erklaerung fuer die W4-v02-Spannung zu behandeln.

Coworker Research / Coworkerz, 2026-09-27.

## 1. Forschungsfrage

Diskriminiert ein von der Helicity-Normierung unabhaengiger Correlation-Ratio-
Schaetzer zwischen den beiden vorab feststehenden Lagen

- H_A: T_BKT in der W4-v02-Lage [0.557, 0.573] (PHY046: 0.5649 +/- 0.0040),
- H_B: T_BKT in der Literatur-Lage [0.570, 0.580] (Okabe/Otsuka 0.573 grob;
  Jiang 0.572(3) NN / 0.576(4) Helicity)?

Die Intervalle ueberlappen. v03 ist daher kein Ja/Nein-Test einer einzelnen
Zahl, sondern ein unabhaengiger Schaetzer mit vorab festgelegtem
Unsicherheits- und Inkonklusiv-Regime.

## 2. Observable

Fuer den Honeycomb-Torus mit Zellindizes (i,j), Subgitter s in {A,B}:

g_d(r) = mean_{i,j,s} cos(theta[i,j,s] - theta[i+r_d,j+r'_d,s])

mit den beiden primitiven Zellrichtungen d=a1,a2. Pro Konfiguration wird

g(r) = (g_a1(r) + g_a2(r))/2

gebildet. Der **Correlation Ratio** ist ein Verhaeltnis von
Ensemble-Mitteln, nicht das Mittel einzelner Quotienten:

R_L(T) = <g(L/2)> / <g(L/4)>.

Nur L durch 4 teilbar ist zulaessig. Gleiche Subgitter werden verglichen;
dadurch entsteht kein A/B-Basisversatz. Der absolute metrische Faktor der
primitiven Honeycomb-Vektoren kuerzt sich im Distanzverhaeltnis heraus.

Diese Observable benutzt weder Upsilon_site noch Upsilon_area und ist damit
gegen den O1-Normierungskonflikt orthogonal.

## 3. Sampler und Vorverarbeitung

- kanonischer Wolff-Single-Cluster-Kernel aus PHY045, unveraendert;
- geordneter Start;
- ein Sweep: Cluster-Updates bis mindestens N geflippte Spins;
- BLAS-Threads = 1;
- keine Parallel-Tempering-Abhaengigkeit; jede (L,T,Seed)-Kette unabhaengig;
- Messung erst nach Thermalisierung;
- Correlations werden aus der aktuellen Konfiguration direkt berechnet.

## 4. Leiter, Temperaturgitter und Statistik

Primaere Leiter:
L = {48, 72, 96, 144, 192}.

Sie repliziert die Groessen der externen Correlation-Ratio-Studie und erlaubt
den direktesten methodischen Vergleich. Alle L sind durch 4 teilbar.

T-Gitter:
0.540 .. 0.610 in Schritten von 0.0025 (29 Punkte).

Statistik:
- n_seeds = 12;
- n_therm = 1000 Sweeps;
- n_meas = 4000 Sweeps;
- Correlation-Ratio wird aus Seed-weisen Ensemble-Mitteln aufgebaut;
- Seed-Vertrag:
  49_000_000 + 1000*L + 100*t_idx + s.

Budget-Stop:
Wenn die vorregistrierte Leiter auf 4 Prozessen > 24 h Wall benoetigt, werden
fehlende Groessen als ungemessen markiert. Parameter werden nicht nach Sicht
auf R(T) geaendert.

## 5. Primaerauswertung

### 5.1 Kurven-Splay

Im BKT-Bereich sind R_L(T)-Kurven fuer verschiedene L naeherungsweise
groessenunabhaengig; oberhalb T_BKT trennen sie sich.

Fuer jedes geordnete Paar L1<L2:
D(T)=R_L2(T)-R_L1(T).

T_splay ist der kleinste T-Gitterpunkt, ab dem D(T) fuer mindestens drei
aufeinanderfolgende Punkte und bis zum oberen Fensterrand das **vor Daten fest
fixierte negative Hochtemperatur-Vorzeichen** besitzt: fuer L1<L2 gilt oberhalb
T_BKT bei endlicher Korrelationslaenge R_L2<R_L1, also D(T)<0, mit
|D| > 2 sigma_D. Das Vorzeichen darf nicht aus dem beobachteten Tail gelernt
werden. Fehlende, nicht-endliche oder ungueltige Punkte brechen die Persistenz
(fail closed).

### 5.2 FSS-Schaetzer — durch v03a mechanisch supersediert

Der folgende Absatz dokumentiert den urspruenglich vorregistrierten
Spline-Ansatz. **Operativ bindend ist seit dem pre-data Addendum
`260927 PHI HEX w4 v03a fss estimator hardening.md` der dort definierte
deterministische symmetric cross-size collapse score.** Kein Produktionslauf
darf auf die hier historische Spline-Mechanik zurueckfallen.

Primaer wird die Okabe/Otsuka-Form verwendet:
R(T,L) = F(X),  X = L / exp(c / sqrt(T-T_BKT)), fuer T>T_BKT.

Gemeinsamer nichtlinearer Fit aller L und T in [T_BKT, 0.610]:
- T_BKT frei in [0.540,0.610],
- c > 0,
- F als kubische Spline mit gemeinsamem Glattheitsparameter,
- Seed-Block-Bootstrap (1000 Replikate) fuer Unsicherheit.

Weil eine flexible Spline Fit-Bias erzeugen kann, ist dieser Fit nur
quotierbar, wenn die synthetischen Recovery-Gates in §7 bestanden sind.

### 5.3 Robuste Sekundaerschaetzer

- groesste zwei L: Beginn des persistenten 2-sigma-Splays;
- Leave-one-L-out FSS;
- T-Fenster [0.545,0.605] und [0.550,0.610];
- lineare statt kubische lokale Darstellung von F(X).

Kein einzelner Sekundaerschaetzer darf nachtraeglich zum Primaerschaetzer
erklaert werden.

## 6. Unsicherheit und Entscheidung

sigma_boot = Bootstrap-SE des primaeren FSS-Schaetzers.

sigma_model = halbe Spannweite der quotierbaren vorregistrierten Varianten,
mit Untergrenze 0.002.

sigma_tot = hypot(sigma_boot, sigma_model).

Ausgabe:
- SUPPORTED_LOW, wenn T_hat + 2 sigma_tot < 0.570;
- SUPPORTED_LITERATURE, wenn T_hat - 2 sigma_tot > 0.573;
- OVERLAP, wenn das 95%-Intervall den gesamten vorregistrierten
  Hypothesen-Ueberlapp [0.570, 0.573] schneidet;
- INCONCLUSIVE, wenn ein Gate aus §7 scheitert oder sigma_tot > 0.010.

**Korrektur 2026-09-27 vor Produktionsdaten:** Die fruehere Ein-Punkt-Regel bei
0.570 war logisch zu permissiv und konnte Werte innerhalb H_A∩H_B faelschlich
als SUPPORTED_LITERATURE etikettieren. Diese Fassung supersediert jene Regel.

Diese Labels sind rein interne Diskriminationssemantik, kein externer
Bestwert-Claim.

## 7. Pflicht-Gates vor Physikinterpretation

G0 INPUT:
- exakt die vorregistrierten Metadaten n_seeds=12, n_therm=1000,
  n_meas=4000;
- exakt die persistierten Preflight-Gates
  {VAL_BIT_numba,G1_geometry,G2_aligned_limit,G3_seed_unique,G4_fss_recovery}
  mit literal `true` fuer jeden Eintrag und ohne Zusatz-Gates;
- exakt 5*29*12 erwartete Rohzeilen, keine Zusatz-/Fremdzeilen;
- `complete is True` und `unmeasured == []`;
- alle erwarteten (L,T,Seed)-Jobs mit exakter (L,t_idx,T,s,seed)-Identitaet;
- endliche, numerisch parsebare g(L/4), g(L/2); malformed Werte fail closed
  statt den Adjudikator zu crashen;
- Nenner |<g(L/4)>| > 1e-6.

VAL-BIT BACKEND:
- wenn Numba verfuegbar ist, muss ein identischer Tiny-Wolff-Lauf bei gleichem
  Seed Python und Numba bit-identische Korrelationsreihen liefern; andernfalls
  ist Produktion gesperrt;
- diese Sperre bindet auch den direkten `_job()`-/Worker-Pfad: blosse
  Numba-Verfuegbarkeit darf den Produktionskernel nicht autorisieren;
- der oeffentliche `produce()`-Pfad muss vor Start irgendeines Messjobs den
  vollstaendigen Preflight VAL-BIT + G1-G4 erfolgreich absolvieren und die
  Gate-Ergebnisse im Produkt persistieren.

G1 GEOMETRY:
- Translation um L liefert identische Siteindizes;
- r=L/4 und L/2 sind gleiche-Sublattice-Verschiebungen;
- a1/a2-Symmetrie auf ausgerichteter Konfiguration exakt.

G2 LIMIT:
- bei perfekt ausgerichteten Spins R_L = 1 fuer jedes L.

G3 SEED:
- Seed-Vertrag kollisionsfrei.

G4 SYNTHETIC RECOVERY:
- generierte BKT-artige Testdaten mit bekanntem T0 muessen T0 innerhalb
  0.003 reproduzieren;
- Null-/Lueckenfaelle muessen fail closed bleiben.

G5 POWER:
- sigma_tot <= 0.010.

G6 ROBUSTNESS:
- Leave-one-L-out darf den Primaerschaetzer um hoechstens 0.008 verschieben;
  sonst INCONCLUSIVE.

## 8. Claim-Decke

PHY049 ist ein unabhaengiger normierungsfreier Quercheck. Selbst bei
bestandenen Gates:
- kein "neuer Weltbestwert";
- kein Wegmitteln mit PHY046 ohne Cross-Family-Modell;
- keine Aussage, dass eine Literaturquelle "falsch" sei;
- externe Aussage erst nach Cross-Family-Review.

## 9. Lineage

- W4-v01 / PHY044: WL-Kosten-/Validitaetskalibrierung.
- W4-v02 / PHY046: per-area Helicity, Wolff, vorregistrierter Finding.
- W4-v03 / PHY049: Correlation Ratio, normierungsfrei, unabhaengiger
  Diskriminationstest.


## 10. Preproduction-Integritaetsnachtrag 2026-09-27

Vor jeder PHY049-Produktion gelten zusaetzlich fail-closed:

- PHY050-`assess_production()` akzeptiert fuer eine Produktionsentscheidung
  **exakt 1000** Bootstrap-Replikate. Kleinere Testbudgets bleiben zulaessig
  fuer interne Unit-Tests von `bootstrap_tbkt()`, koennen aber weder G5 noch
  G6 noch eine Physikentscheidung freischalten.
- G0 bindet neben L/T/Seed auch n_therm und n_meas sowie die exakte
  Rohzeilenzahl. Ein formal komplettes Produkt mit abweichendem
  Simulationsbudget ist kein W4-v03-Produkt.
- Persistierte Rohmomente werden als untrusted evidence behandelt:
  fehlende/nichtnumerische/NaN/Inf-Werte fuehren zu G0=FAIL, nicht zu einer
  Exception ausserhalb des Adjudikators.


## 11. Operativer Produktionsvertrag / Checkpointing 2026-09-27

Dieser Nachtrag aendert **keine** Physik-, Leiter-, Seed-, Sweep-, Bootstrap-
oder Entscheidungsregel. Er macht den bereits vorregistrierten 24-h-Lauf
transaktional ausfuehrbar:

- der offizielle CLI-Produktionspfad verwendet ausschliesslich die exakten
  W4-v03-Defaults und schreibt in eine explizite JSON-Ergebnisdatei;
- persistiert wird atomar via Temp-Datei + fsync + Replace, sodass ein
  Prozess-/Hostabbruch kein halb geschriebenes Evidenzartefakt hinterlaesst;
- Checkpoints werden nur nach einem **vollstaendigen L-Block** committed.
  Ein partieller L-Block ist beim Resume unzulaessig und wird nie als Evidenz
  uebernommen;
- Resume akzeptiert nur einen byte-semantisch identischen
  `campaign_contract` (Leiter, T-Gitter, Seeds, Thermalisierung, Messungen,
  Workerzahl, Wall-Budget) und identische gruene Preflight-Gates;
- ein terminaler `WALL_BUDGET_STOP` darf nicht per Resume in eine laengere
  Kampagne verwandelt werden. Fehlende Jobs bleiben explizit ungemessen;
- ein bereits vollstaendiges Checkpoint wird ohne Neuberechnung wiedergegeben;
- der Produktionsrunner selbst erzeugt **keine Physikinterpretation**.
  Erst ein exaktes, vollstaendiges Produkt darf an PHY050
  `assess_production()` uebergeben werden.

Operative Aufrufe:

`python "src/260927 PHY049 honeycomb correlation ratio v01.py"`
fuehrt nur den Preflight aus.

`python "src/260927 PHY049 honeycomb correlation ratio v01.py" --production`
startet die exakte Kampagne und schreibt atomare Whole-L-Checkpoints.

`--resume` ist nur fuer einen nichtterminalen, exakt passenden
Whole-L-Checkpoint zulaessig.
