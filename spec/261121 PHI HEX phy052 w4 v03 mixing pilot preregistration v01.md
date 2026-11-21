# PHY052 - Vorregistrierung: Mixing-/Thermalisierungs-Pilot vor W4 v03

> **Status:** VORREGISTRIERT am 2026-11-21, **vor** jeder Pilot-Datennahme und
> vor jeder PHY049-Produktionsmessung. Aendert W4 v03 / v03a nicht: Leiter,
> T-Gitter, Seeds, Sweeps und Adjudikation von PHY049/PHY050 bleiben exakt.
>
> **Anlass:** Cross-Family-Gesamtaudit vom 2026-10-04 (007acc2,
> `ANTWORT_007acc2_PhiHex_Gesamtaudit_Hardening_2026-10-04`, Roadmap P1 6-8):
> W4 v03 startet jede Kette geordnet. Ein gemeinsamer Thermalisierungs-Bias
> aller Ketten waere fuer Diagnostik zwischen Ketten unsichtbar. Vor der
> 5x29x12-Kampagne ist deshalb ein Pilot mit kalten UND heissen Starts
> verlangt; nur bei PASS darf PHY049 messen.

Coworker Research / Coworkerz, 2026-11-21.

## 1. Frage

Erreichen Ketten des unveraenderten PHY049-Kernels (Wolff, PHY045) bei den
vorregistrierten W4-v03-Sweepzahlen 1000 + 4000 aus geordnetem und aus
gleichverteilt zufaelligem Start dieselbe stationaere Verteilung des
entscheidungsrelevanten Funktionals R = <g(L/2)> / <g(L/4)>, mit genug
effektiver Stichprobe?

## 2. Design (fix)

- Zellen: L in {48, 96, 192} x T in {0.560, 0.575, 0.595} (9 Zellen; kleinste,
  mittlere, groesste Leiterstufe; unterhalb, nahe und oberhalb der
  Hypothesenlagen). Alle T liegen auf dem W4-v03-Gitter.
- Je Zelle 4 unabhaengige Ketten: s = 0, 2 geordnet (theta = 0),
  s = 1, 3 gleichverteilt in [-pi, pi).
- Je Kette 1000 Thermalisierungs- + 4000 Mess-Sweeps, Messung g(L/4), g(L/2)
  wie PHY049 (identischer, VAL-BIT-gebundener Kernel).
- Seeds: Mess-Seed `79_000_000 + 1000*L + 100*t_idx + s`,
  Init-Seed `89_000_000 + 1000*L + 100*t_idx + s` (t_idx = Index im
  W4-v03-Gitter). Beide Namensraeume sind von der Produktion
  (`49_000_000 + ...`) disjunkt: Pilotdaten sind fuer PHY050 nicht verwertbar.
- BLAS 1 Thread; Prozess-Pool ueber `PHY045.map_jobs`.

## 3. Diagnostik (fix)

Gate-Observable je Zelle: die Einflussreihe des Ratio-of-Means
u_t = (h_t - R q_t) / <q>, q = g(L/4), h = g(L/2), R und <q> gepoolt ueber
alle vier Ketten. Var(R_hat) ist in erster Ordnung Var(mean u); u ist also
genau die Groesse, deren Durchmischung PHY050 voraussetzt.

- **rhat:** Maximum aus rank-normalisiertem split-Rhat und gefaltetem
  split-Rhat (Vehtari et al. 2021), Definition und Rechenweg exakt wie
  `stan-dev/posterior` `rhat()` (Normal-Score (r - 3/8)/(S + 1/4), Faltung um
  den Median aller Ziehungen, Split in Haelften).
- **ess_bulk:** Multi-Ketten-ESS der rank-normalisierten Split-Ketten mit
  Geyers initial-positive und initial-monotone Folge und "improved
  estimate", exakt wie `posterior` `ess_bulk()`; Autokovarianz biased (1/N).
- Ungueltige oder konstante Reihen -> NaN -> Gate rot (fail-closed).

Korrektur gegenueber dem Audit-Bundle v03 (vor Daten, nicht nach Sicht):
dessen Normal-Score-Nenner S - 1/4 und Autokovarianz 1/(N-k) weichen von der
Referenz ab; massgeblich ist die Referenz. Die Portierung ist gegen ArviZ
(arviz-stats 1.3.3) auf 24 synthetischen Faellen bis auf Rundung gleich
(Ausnahme: ungerade Kettenlaenge, dort folgt PHY052 `posterior`
`.split_chains`, das die mittlere Ziehung verwirft; im Pilot ist N = 4000).

## 4. Entscheidungsregel (fix)

Eine Zelle besteht, wenn **rhat(u) <= 1.01** und **ess_bulk(u) >= 400**
und **|z_kh| <= 3.5**.

- 1.01 und 400 (= 100 je Kette bei 4 Ketten) sind die Empfehlung von
  Vehtari et al. (2021).
- z_kh = (R_kalt - R_heiss) / hypot(SE_kalt, SE_heiss): Ratio-of-Means der
  zwei geordnet bzw. der zwei zufaellig gestarteten Ketten; SE je Gruppe =
  sd(u_Gruppe) / sqrt(ess_bulk(u_Gruppe)). Das ist der klassische
  Zweistart-Test der statistischen Physik. Er ist noetig, weil rhat einen
  Versatz einer Kette erst ab etwa 0.3 Standardabweichungen je Ziehung
  ueber 1.01 hebt (gemessen vor den Daten: Versatz 0.15 sigma in einer von
  vier Ketten a 4000 -> rhat 1.0027), waehrend fuer W4 v03 Bias relativ zum
  Standardfehler des Mittelwerts zaehlt. 3.5: zweiseitig, unter H0
  Fehlalarm je Zelle 4.7e-4, ueber 9 Zellen ~0.4 %. Nicht endliches z ->
  Zelle rot.
- Aufloesung (ehrlich): der Pilot schliesst einen Thermalisierungs-Bias
  oberhalb von etwa 3.5 * hypot(SE_kalt, SE_heiss) aus; kleinere Biases
  liegen unter seiner Aufloesung. Die Produktion (12 Seeds) hat einen
  kleineren Standardfehler; der Pilot ist daher eine notwendige, keine
  hinreichende Bedingung. Die tatsaechliche Schranke je Zelle steht im
  Report (`cold_minus_hot_se`).

- `MIXING_PILOT_PASS`: exakt das Design aus Abschnitt 2 (keine Zusatz- oder
  Fehlzellen, keine verkleinerten Laeufe) und alle 9 Zellen bestehen.
- sonst `MIXING_PILOT_FAIL`: PHY049 startet nicht. Eine Aenderung von
  n_therm/n_meas ist dann nur ueber eine neue, datierte Vorregistrierung
  (W4 v03b) moeglich; die Schwellen hier werden nach Sicht der Pilotdaten
  nicht gelockert.

Nur berichtet, nicht entscheidend: rhat/ess_bulk von g(L/4) und g(L/2)
einzeln.

Nachtrag am selben Tag, VOR jeder Pilot-Datennahme: der Zweistart-Test war
im ersten Entwurf nur berichtet; ein synthetischer Selbsttest (siehe oben)
hat gezeigt, dass rhat allein zu grob ist. Er ist deshalb Gate.

## 5. Freigabe der Produktion (fix)

`ELIGIBLE_TO_START_PHY049_MEASUREMENT` genau dann, wenn
(1) PHY052 `MIXING_PILOT_PASS` und (2) der native PHY049-Preflight exakt die
fuenf Gates VAL_BIT_numba, G1_geometry, G2_aligned_limit, G3_seed_unique,
G4_fss_recovery literal true traegt und `production_measurement_eligible`
true ist. Sonst `DO_NOT_START_PHY049_MEASUREMENT`. Die Freigabe ist kein
T_BKT-Urteil und schaltet PHY050 nicht frei.

## 6. Claim-Decke

Ein PASS belegt Stationaritaet und Durchmischung des Kernels bei diesen
Sweepzahlen in 9 Zellen, nicht die Lage von T_BKT und nicht die Gueltigkeit
der FSS-Form. Pilotdaten gehen in keine Physikzahl ein.

## 7. Maschinelle Bindung

Die Konstanten `PILOT_*`, `MAX_RHAT`, `MIN_ESS_BULK`, `MAX_COLD_HOT_Z` in
`src/261121 PHY052 honeycomb correlation-ratio mixing pilot v01.py` werden per
Test gegen diese Spec und gegen W4 v03 (n_therm, n_meas, T-Gitter, Leiter,
Seed-Disjunktheit) geprueft.

## 8. Quellen

- Vehtari, Gelman, Simpson, Carpenter, Buerkner, "Rank-normalization,
  folding, and localization: An improved R-hat for assessing convergence of
  MCMC", Bayesian Analysis 16, 667-718 (2021), doi:10.1214/20-BA1221.
- Geyer, "Practical Markov Chain Monte Carlo", Statistical Science 7,
  473-483 (1992).
- Referenzimplementierung `stan-dev/posterior`, R/convergence.R (abgerufen
  2026-11-21): z_scale, backtransform_ranks, fold_draws, .split_chains,
  .rhat, .ess, rhat, ess_bulk.
