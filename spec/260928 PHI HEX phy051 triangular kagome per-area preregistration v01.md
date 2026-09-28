# PHY051-Vorregistrierung v01: triangular + kagome T_BKT mit Upsilon PRO FLAECHE - O1-Diskriminator gegen normierungsfreie Referenzen

> **Status:** VORREGISTRIERT am 2026-09-28, **vor** jeder PHY051-Datennahme
> (Commit-Zeitstempel = Beleg). Setzt den in PR #50 offen gelassenen Punkt
> "triangular und kagome pro Flaeche neu messen" (README "Naechste Stufe
> 2026-09-26b", Punkt 3; O1-Entscheid §5) um. Vor diesem Commit lief nur ein
> reiner Timing-Pilot (Wall je Sweep, keine Upsilon-Werte; §4).
>
> Vertragsquellen: `spec/260926 PHI HEX O1 helicity normalization decision
> v01.md` (Normierung), `spec/260926 PHI HEX w4 honeycomb preregistration
> v02.md` (Protokoll-Vorlage), `results/260926 PHY046 W4 v02 interpretation
> note.md` §2 (Lehre PHY048: Paar-Mittel statt 3-Parameter-HKS als Primaer).

Coworker Research / Coworkerz, 2026-09-28.

## 1. Frage

Audit O1 ist auf der Ebene der Steifigkeits-Definition entschieden (PHY045:
per Flaeche). Offen ist die LAGE-Konsequenz auf den beiden Gittern, deren
Referenzwerte NICHT aus einer Helicity-Messung stammen:

- **triangular** (a_s = sqrt(3)/2 = 0.8660 < 1): Die per-Site-Pipeline gab
  1.4007 (PHY030 v02, WM-Fit) und traf damit die Helicity-Referenz 1.418(2)
  (Sorokin, zitiert ueber arXiv:2305.00651; deren Normierung ist nicht
  belegt). Die normierungsFREIEN Referenzen liegen deutlich hoeher:
  Hochtemperatur-Reihe Butera & Pernici (arXiv:0806.1496) J_c = 0.6824(8),
  also T_BKT = 1.4654(17); MC-Quercheck arXiv:1010.3075 J_c = 0.6833(6),
  also 1.4635(13). PHY047 (per Flaeche, nur L <= 19) gab 1.4447-1.4496.
- **kagome** (a_s = 2/sqrt(3) = 1.1547 > 1): per Site 0.8479 (PHY033);
  normierungsfreie Referenz 0.825 (arXiv:2501.07388 Tab. 1, "rough
  estimate", ohne Fehlerbalken); PHY047 per Flaeche (L <= 36) 0.819-0.820.

Beide Gitter bilden zusammen einen Vorzeichen-Test: per Site verzerrt
triangular nach UNTEN (a_s < 1) und kagome nach OBEN (a_s > 1). Liegt der
per-Flaeche-Kanal auf beiden Gittern im normierungsfreien Band und der
per-Site-Kanal ausserhalb, ist der O1-Entscheid auf der Lage-Ebene
korroboriert. Der umgekehrte Fall waere ein Signal GEGEN den O1-Entscheid.

**Claim-Decke:** FINDING relativ zu vorregistrierten Baendern; kein
T_BKT-Bestwert; keine Aussage ueber die Richtigkeit einzelner
Literaturwerte; externe Aussage erst nach Cross-Family-Review.

## 2. Observable (identisch zu W4 v02)

Upsilon_A = [(T1_x - beta S_x^2) + (T1_y - beta S_y^2)] / (2 A), J = 1,
NN-Abstand 1, per-Bond-Twist-Schaetzer (PHY045-Kernel). Torus-Flaeche
A = N * a_s. Per Site (Upsilon_A * a_s) wird aus denselben Samples als
Konventions-Quercheck mitgefuehrt und dient dem O1-Diskriminator (§6).

Exaktes T=0-Orakel (Gate A, MC-frei): Upsilon_site(0) = 1.5 (triangular),
1.0 (kagome); Upsilon_A(0) = sqrt(3) bzw. sqrt(3)/2.

## 3. Sampler, Leitern, T-Gitter, Statistik

- **Sampler:** kanonischer Wolff-Single-Cluster (Numba-Kernel PHY045,
  Fallback reines Python), 1 Sweep = Cluster bis >= N Spins geflippt;
  geordneter Start. BLAS auf 1 Thread (AGENTS.md).
- **triangular** (Core-Torus `build_triangular_lattice(r, periodic=True)`,
  L = 2r+1 ungerade, Rhombus): Leiter L in {33, 49, 65, 97, 129, 193, 257},
  Paare (L1, L2) = (33,65), (49,97), (65,129), (97,193), (129,257). Die
  Paare sind (L, ~2L) mit L2/L1 in [1.97, 1.99]; der C-eliminierte
  Paar-Schaetzer nutzt ln(L2/L1) exakt, exakte Verdopplung ist nicht
  noetig. T-Gitter 1.380 .. 1.510, Schritt 0.005 (27 Punkte): deckt den
  per-Flaeche-Erwartungsbereich (1.44 .. 1.49) UND den per-Site-Bereich
  (~1.38 .. 1.42) ab.
- **kagome** (`build_kagome_lattice(L)`, N = 3 L^2):
  Leiter L in {32, 48, 64, 96, 128, 192}, Paare (32,64), (48,96),
  (64,128), (96,192).
  T-Gitter 0.780 .. 0.860, Schritt 0.005 (17 Punkte).
- **Statistik je (L, T):** 8 unabhaengige Seeds, je 500 Thermalisierungs-
  und 1500 Mess-Sweeps (12000 Mess-Sweeps je Punkt, wie W4 v02).
- **Seed-Vertrag:** `51_000_000 + 1_000_000*lat_id + 1000*L + 10*t_idx + s`
  mit lat_id = 0 (triangular), 1 (kagome). Disjunkt zu PHY045 (45M),
  PHY045-C (46M) und PHY046 (47M); < 2^32.
- **Fehler:** Jackknife ueber Seeds (leave-one-seed-out, alle L gemeinsam)
  fuer jeden abgeleiteten Wert.
- **Reihenfolge / Budget (S0):** je Gitter L-Bloecke absteigend; ueberschreitet
  die Wall-Zeit 3.0 h auf 4 Prozessen, werden die restlichen L als
  `unmeasured` berichtet (NaN), es wird nicht umgeplant. Ladder unvollstaendig
  -> NEGATIVE_RESULT (Regel S0).

## 4. Timing-Pilot (blind, 2026-09-28, diese Maschine, 1 Prozess)

Wall je Sweep mit dem PHY045-Numba-Kernel (30+30 Sweeps, keine Upsilon-
Werte ausgewertet): triangular L=65 0.93 ms, L=129 3.6 ms, L=193 7.1 ms,
L=257 15.8 ms; kagome L=32 0.49 ms, L=64 1.8 ms, L=128 8.1 ms, L=192
17.4 ms. Projektion fuer §3: triangular ~3.6 CPU-h, kagome ~2.5 CPU-h,
zusammen ~1.6 h Wall auf 4 Prozessen - innerhalb S0.

## 5. Schaetzer (Lehre PHY048 umgesetzt)

- **Paar-Crossing T*(L1, L2):** C-eliminierte WM-Bedingung
  (`tbkt_pair_from_curves`, PHY040; R > 0-Guard, nur benachbarte
  Gitterpunkte) auf den Seed-Mittelkurven von Upsilon_A.
- **Primaer-Schaetzer T_P:** inverse-varianz-gewichtetes Mittel der
  Paar-Crossings mit L1 >= L_min (triangular L_min = 65: Paare (65,129),
  (97,193), (129,257); kagome L_min = 64: Paare (64,128), (96,192)),
  Gewichte = Paar-Jackknife-SE aus dem Vollsatz. Fehlt eine SE (oder ist
  sie 0), wird ungewichtet gemittelt. Weniger als 2 solche Paare ->
  fail-closed NEGATIVE_RESULT (Regel S1b).
- **FSS-Varianten fuer sigma_FSS** (vorab fixiert): (v1) gewichtetes
  Mittel ALLER Paare; (v2) groesstes Paar allein; (v3) HKS
  T*(L1) = T_c + a/ln^2(b L1), 3 Parameter, b in [0.1, 10]; (v4) HKS mit
  b = 1; (v5) Weber-Minnhagen-Fit mit freiem C ueber alle L bei festem T
  (chi^2-Minimum ueber das T-Gitter). Alle Varianten 1:1 aus PHY046.

## 6. sigma-Budget, Baender, Entscheidungsregeln

- sigma_sampler = Jackknife-SE von T_P.
- sigma_FSS = halbe Spannweite ueber {T_P, v1..v5}, Untergrenze **0.5 %
  von T_P** (relatives Aequivalent des v02-Floors 0.003/0.57).
- sigma_tot = hypot(sigma_sampler, sigma_FSS).
- **Baender (literaturgestuetzt, normierungsfrei):**
  - triangular B_tri = **[1.450, 1.480]** (~ +/- 1 % um 1.465; enthaelt
    1.4654(17) und 1.4635(13)).
  - kagome B_kag = **[0.808, 0.842]** (~ +/- 2 % um 0.825, weil die Referenz
    eine grobe Schaetzung ohne Fehlerbalken ist).
- Regeln je Kanal, in dieser Reihenfolge: **S0** Ladder unvollstaendig;
  **S1** weniger als 3 Paare mit Crossing; **S1b** weniger als 2 Primaer-
  Paare; **S2** sigma_tot > 1 % von T_P -> NEGATIVE_RESULT. Sonst **C**:
  [T_P - 2 sigma_tot, T_P + 2 sigma_tot] schneidet B -> CONSISTENT;
  **I**: disjunkt -> INCONSISTENT.
- **O1-Diskriminator (Wahrheitstafel, je Gitter):**

  | per Flaeche | per Site | Label |
  |---|---|---|
  | CONSISTENT | INCONSISTENT | `O1_PER_AREA_CORROBORATED` |
  | INCONSISTENT | CONSISTENT | `O1_CHALLENGED` |
  | CONSISTENT | CONSISTENT | `NON_DISCRIMINATING` |
  | INCONSISTENT | INCONSISTENT | `TENSION_BOTH_CHANNELS` (NR) |
  | NEGATIVE_RESULT in einem Kanal | - | `NEGATIVE_RESULT` |

  Zusaetzlich (Diagnostik, nicht entscheidend): z-Abstand von T_P (per
  Flaeche) zur HT-Reihe 1.4654(17) bzw. zu 0.825, und z-Abstand des
  per-Site-Kanals zur Helicity-Referenz 1.418(2) (triangular).

## 7. Pipeline-Gates (Integritaet, fails-closed, KEINE Physik-Wahrheit)

| Gate | Bedingung |
|---|---|
| PASS_INPUT_COMPLETE | alle (L, T, Seed)-Werte vorhanden und finite, Form = Plan |
| PASS_T0_GEOMETRY_ORACLE | exaktes Upsilon_site(0) je Leiter-L (1.5 / 1.0) auf 1e-12; a_s wie §2 |
| PASS_LOW_T_MERGE | bei T_min: max. relative Paar-Abweichung der Upsilon_A-Mittelkurven ueber alle L < 8 % |
| PASS_HIGH_T_SPLAY | bei T_max: Upsilon_A(L_max) < Upsilon_A(L_min) |
| PASS_CURVE_SE_BOUNDED | max. Seed-SE von Upsilon_A ueber das Gitter < 2 % von 2 T_max / pi |
| PASS_MIN_PAIRS | >= 3 Paare mit Crossing (per Flaeche) |

OVERALL = alle Gates. Ein FAIL schreibt trotzdem Report + Evidenz.

## 8. Evidenz und maschinelle Bindung

- Reports: `results/260928 PHY051 triangular area-helicity report.{json,txt}`
  und `results/260928 PHY051 kagome area-helicity report.{json,txt}`; JSON
  enthaelt die per-Seed-Rohdaten beider Kanaele.
- `tests/test_phy051_area_crosscheck.py`: Konstanten zeichengenau an diese
  Spec gebunden (Leitern, Paare, T-Gitter, Seeds, Sweeps, Seed-Vertrag,
  Baender, Floors); Seed-Disjunktheit; T=0-Orakel; Schaetzer-Orakel auf
  synthetischen WM-Kurven; Wahrheitstafel §6 fails-closed; Mini-Produktion
  end-to-end; committete Reports (falls vorhanden) gegen den Plan.
- Kein Parameter dieses Vertrags darf nach Sicht der PHY051-Daten geaendert
  werden (sonst v02 als datierter post-hoc-Nachtrag).
