# W4-Vorregistrierung: honeycomb T_BKT gegen ein fixes Literatur-Konsistenzband

> **Status:** VORREGISTRIERT am 2026-09-26, **vor** der Kalibrierung PHY044
> (K1-K3) und vor jedem W4-Produktionslauf. Diese Datei ist per SHA-256 in
> `SOURCES.md` gepinnt; jede Aenderung erzeugt eine neue, sichtbare
> SOURCES-Zeile. Aenderungen nur als versionierter Nachtrag (v02, ...) mit
> Begruendung und dem Etikett **post-hoc**, falls W4-Produktionsdaten schon
> gesehen wurden. Nie still.
>
> Anlass: Issue #45 §1 (Machbarkeitsstudie W4, 2026-09-11). Vertragsquellen
> fuer Referenzen: `spec/260703 ... reference conventions audit v01.md` +
> `spec/260926 ... reference provenance addendum v01.md`.

Coworker Research / Coworkerz, 2026-09-26.

## 1. Frage (umformuliert)

Das urspruengliche Ziel (i) - 0.573 (arXiv:2501.07388, ohne Fehlerbalken) von
0.576(3) (arXiv:2406.12076, WL) mit 3 sigma zu trennen - ist **gestrichen**: es
braucht sigma_T <= 0.001, waehrend schon die FSS-Modellwahl bei L <= 128
sigma_T >~ 0.003 erzeugt (HKS arXiv:1302.2900: Polynom- vs ln^2-Form 7.5e-4
relativ selbst bei L <= 512; arXiv:2406.12076: Upsilon 0.576(3) vs Upsilon_4
0.568(1) in derselben Arbeit).

**W4-Frage (ii):** Ist der vorab festgelegte Upsilon_2-Schaetzer (§3) bei
explizit ausgewiesenem Gesamt-sigma (§6) mit dem festen Konsistenzband B (§2)
vertraeglich?

## 2. Konsistenzband B (fix)

**B = [0.560, 0.580]** (T-Form, J = k_B = 1).

- Untere Grenze 0.560: arXiv:2406.14812 v1, NN-Kanal (kleinster direkt
  berichteter T-Wert im Band).
- Obere Grenze 0.580: Binder-beta-Kanal von arXiv:2406.12076 in T-Form
  (0.5800(7)); alle direkt berichteten T-Werte (0.560 .. 0.576, Journal-
  Fassungen 0.572 .. 0.576) liegen darin.
- **Nicht** in B: die beta-Kanaele 0.5928(11) (Upsilon) und 0.6116(41)
  (Upsilon_4). Sie sind eine quellen-interne Spannung (Addendum §2.1). Ein
  W4-Ergebnis, das mit 0.5928 vertraeglich, mit B aber nicht vertraeglich ist,
  wird **genau so** berichtet - B wird deswegen nicht nachtraeglich verschoben.

## 3. Primaer-Observable und Primaer-Schaetzer (fix)

- **Observable:** Upsilon_2 per Site (Orakel honeycomb Upsilon_2(0) = 3/4 J,
  Gate A aus PHY042), Wang-Landau/1-t-Kernel PHY041 (oder sein bit-identischer
  Numba-Zwilling aus PHY044, VAL-BIT-Gates vorausgesetzt).
- **Kurve je L:** Walker-Mittel der Upsilon_2(T)-Kurven (PHY042).
- **Primaer-Schaetzer T_W4:** C-eliminiertes WM-Paar (`tbkt_pair_from_curves`,
  PHY040) des **groessten quotierbaren Paares** (§4) - "groesst" = groesstes
  L_a, bei Gleichstand groesstes L_b.
- **Sekundaer (Querchecks, keine Band-Entscheidung):** Upsilon_4-Dip-Lage je L
  (nur walker-robuste Dips, Kriterium PHY042: Walker-Dip-Spread <= 0.005);
  konventionsfreie xi_2/L-Splays (PHY043-Methodik) in einem separaten Modul.

**L-Leiter:** {32, 48, 64, 96}. L = 128 nur, wenn die K1-Projektion (§5) die
Kosten je Walker auf <= 24 h (Numba, Referenzmaschine des Kalibrierlaufs)
beziffert. L <= 24 ist **nicht** Teil der Primaer-Leiter (PHY041-Befund:
kleine-L-Bias dominiert).

## 4. Walker-, Domaenen- und Quotierbarkeits-Regeln (fix)

1. **>= 3 unabhaengige g(E)-Walker an JEDEM L** der Leiter. Maschinell:
   `_walker_plan(..., min_walkers=3)` bzw. ein W4-Walker-Plan, der dieselbe
   Pruefung fail-closed vor jeder Rechnung macht. Kein konstruierter
   spread = 0, kein "Domaene formal voll" (Issue #45 §2).
2. **Validitaets-Domaene je L:** zusammenhaengend ab dem unteren T-Rand,
   solange der Walker-Spread von Upsilon_2 < 0.04 bleibt (PHY042-Semantik,
   `_walker_domain`).
3. **Quotierbar** ist ein Paar nur mit Basis `both_measured` und Crossing
   <= min(T_max(L_a), T_max(L_b)) (`_pair_quotable`). Basis `partial`,
   `unmeasured` oder `empty` ist in W4 **nie** quotierbar.
4. **Domaenen-Anforderung an das Budget:** T_max(L) >= **T_req = 0.62** fuer
   jedes L der Leiter. Begruendung: 0.62 liegt ueber allen Referenzkanaelen
   inkl. 0.6116 + 0.0041 - die Anforderung haengt damit nicht davon ab, wo das
   Crossing liegt.

## 5. Kalibrier-Protokoll PHY044 (blind, vor W4)

Ziel: Budgets festlegen, nicht Physik messen. Modul
`src/260926 PHY044 honeycomb wl calibration v01.py`.

- **VAL-BIT:** der Numba-Kernel reproduziert (S) kleine Gitter bit-exakt und
  (P) alle 7 committeten PHY042-Walker. Ohne PASS kein Einsatz des
  Numba-Kernels in W4.
  *Praezisierung 2026-09-26, nach dem ersten VAL-BIT-P-Lauf und VOR K1-K3
  (keine Physik gesehen):* gegen die auf einer anderen Plattform committete
  PHY042-Evidenz sind die Sampler-Trajektorien exakt gleich (wl_sweeps und
  Bin-Belegung aller 7 Walker), die Y2/Y4-Kurven weichen aber auf ULP-Ebene
  ab (Y2 <= 1 ulp, Y4 rel. ~1e-13) - Reduktionsreihenfolge von np.dot/np.sum
  je Plattform/numpy-Version, nicht Sampler. (P) prueft deshalb Trajektorie
  exakt und Kurven bis rtol 1e-12; die strikte Bit-Gleichheit Original vs
  Numba wird auf DERSELBEN Maschine an einem identischen L=24-Produktionsjob
  gepinnt (Stufe `speed`).
- **K1 Kosten:** L in {24, 32, 48, 64}, je 3 Walker, PHY042-Budget-Rezept
  (lnf_final = 1e-5, prod = max(30000, 60 * nbins)); gemessen: wl_sweeps,
  Spin-Updates, Wall-Zeit (Numba) je Walker; Potenzgesetz-Fit der Spin-Updates
  in L; Python-Wall-Zeit ueber den gemessenen Speedup auf einem identischen
  L=24-Job.
- **K2 Praezision:** T_max(L) je L; Walker-Streuung (max - min) der Paar-
  Crossings ueber alle Walker-Kombinationen - **nur Streuung, nie Lage**.
- **K3 Hebel (L = 48):** (a) 4x Produktion, (b) lnf_final = 1e-6, (c) Zerlegung:
  eine g(E), 3 unabhaengige Produktionen (Produktions-Anteil des Spreads).
- **Blind-Vertrag:** kein T_BKT-Lagewert in Konsole, JSON oder Text-Report;
  maschinell `assert_blind` (fails-closed) vor dem Schreiben.
- **Budget-Entscheidungsregel:** Budget(L) := kleinstes kalibriertes Budget mit
  T_max(L) >= T_req. Erreicht bei L = 48 **kein** kalibriertes Budget T_req,
  ist W4-Produktion bei L >= 48 **NO-GO** (Stop-Regel S0), bis ein weiterer
  Hebel kalibriert ist. Budgets fuer L = 64/96 werden aus dem K3-Hebel mit dem
  K1-Exponenten projiziert und im ersten W4-Lauf durch dessen eigene
  Domaenen-Messung geprueft (Regel 4.4 bleibt bindend).

## 6. sigma-Budget (Gesamt-Unsicherheit, fix)

sigma_tot = sqrt(sigma_sampler^2 + sigma_FSS^2)

- **sigma_sampler:** halber Walker-Kombinations-Spread des Primaer-Paares
  (PHY042-Konvention).
- **sigma_FSS (Modellwahl):** halbe Spannweite ueber drei vorab fixierte
  Schaetzer-Varianten: (a) T_W4; (b) alle anderen quotierbaren Paare;
  (c) WM-Fit mit freiem C ueber alle L der Leiter mit gemessener Domaene.
  **Untergrenze sigma_FSS >= 0.003** (Issue #45 §1).
- Seed-/Bootstrap-CIs werden **nie** als Gesamt-CI bezeichnet.

## 7. Entscheidungsregeln (fix, maschinell `w4_verdict`)

In dieser Reihenfolge:

1. **NEGATIVE_RESULT (S1):** weniger als **3** quotierbare Paare (Basis
   `both_measured`) -> keine Band-Aussage.
2. **NEGATIVE_RESULT (S2):** sigma_tot > **0.010** -> keine Band-Aussage.
3. **CONSISTENT:** [T_W4 - 2 sigma_tot, T_W4 + 2 sigma_tot] schneidet B.
4. **INCONSISTENT:** das Intervall ist disjunkt zu B -> FINDING (Spannung),
   ohne Nachjustieren von B, Leiter oder Schaetzer.

Zusaetzlich berichtet (nicht entscheidend): Lage relativ zu den beta-Kanaelen
0.5928(11) / 0.6116(41); die Paar-Sequenz in L (Drift-Diagnose).

## 8. Claim-Decke und Verbote

- Kein "Bestwert", keine Diskriminierung 0.573 vs 0.576.
- Kein Mischen von beta- und T-Kanaelen (260703-Regel).
- Keine Parameter-Aenderung nach Sicht auf W4-Produktionsdaten ausser per
  versioniertem post-hoc-Nachtrag; post-hoc-Analysen separat berichtet.
- Kanaele mit Status `search_corroborated`/`version_unclear` (Addendum §4)
  tragen allein keinen Claim; vor einer externen Aussage Primaertext-Abgleich
  und Cross-Family-Review.

## 9. Maschinelle Bindung

- Konstanten `W4_*` und `w4_verdict` in `src/260926 PHY044 honeycomb wl
  calibration v01.py`; `tests/test_phy044_wl_calibration.py` prueft, dass die
  Zahlenwerte dieser Spec (B, 3 Walker, 3 Paare, 0.010, 0.003, 0.04, 0.62)
  zeichengenau mit den Konstanten uebereinstimmen, und testet `w4_verdict`
  gegen Orakel-Faelle jeder Regel.
- PHY042 `_walker_plan(min_walkers)`, `_walker_domain`, `_pair_quotable`
  (Issue #45 §2, fail-closed) sind die Domaenen-Primitiven.
