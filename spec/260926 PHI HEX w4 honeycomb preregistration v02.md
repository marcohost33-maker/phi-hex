# W4-Vorregistrierung v02: honeycomb T_BKT - Upsilon pro Flaeche, kanonischer Wolff, HKS-Paare

> **Status:** VORREGISTRIERT am 2026-09-26, **vor** jeder W4-Datennahme.
> Ersetzt v01 (`spec/260926 PHI HEX w4 honeycomb preregistration v01.md`) in
> den unten genannten Abschnitten. v01 bleibt als Lineage stehen; es wurden
> keine W4-Daten unter v01 genommen.
>
> **Anlass (Selbstpruefung, `spec/260926 PHI HEX PR50 self-review v01.md`):**
> (1) v01 verwendete Upsilon pro Site. Audit O1 ist entschieden
> (`spec/260926 PHI HEX O1 helicity normalization decision v01.md`): das
> NK-Kriterium verlangt Upsilon pro Flaeche. (2) Die blinde Effizienzmessung
> K4 (PHY045 Teil C) zeigt, dass kanonischer Wolff fuer die W4-Genauigkeit um
> Groessenordnungen billiger ist als der Wang-Landau-Pfad. (3) Der Primaer-
> Schaetzer von v01 (groesstes Paar, ohne Extrapolation) entspricht nicht der
> Best Practice (HKS 2013).

Coworker Research / Coworkerz, 2026-09-26.

## 1. Frage und Band (unveraendert aus v01, Begruendung korrigiert)

Frage (ii): Ist der vorab festgelegte Schaetzer (§4) bei ausgewiesenem
Gesamt-sigma (§5) mit **B = [0.560, 0.580]** vertraeglich?

*Korrigierte Begruendung von B:* B stammt aus dem Vorschlag in Issue #45.
Seine Grenzen fallen mit arXiv:2406.14812 v1 (NN 0.560, Journal-Fassung
0.572) und dem Binder-beta-Kanal von arXiv:2406.12076 (T-Form 0.5800)
zusammen. v01 begruendete die obere Grenze mit einem beta-Kanal und schloss
zugleich die beiden anderen beta-Kanaele aus - das war inkonsistent. B bleibt
unveraendert, weil es literatur- und nicht datengestuetzt ist. Zusaetzlich
(nicht entscheidend) wird die Lage relativ zu den Journal-Werten
0.572-0.576, zum normierungsfreien Correlation-Ratio-Wert 0.573
(arXiv:2501.07388) und zu den beta-Kanaelen berichtet.

## 2. Observable (neu)

**Upsilon pro Flaeche:**
Upsilon_A = [(T1_x - beta S_x^2) + (T1_y - beta S_y^2)] / (2 A),
mit A = L^2 * 3 sqrt(3)/2 (honeycomb-Torus, NN-Abstand 1) und J = 1.
Der per-Sample-Schaetzer ist erwartungstreu, weil beide Terme lineare
Erwartungswerte sind. Pro Site wird nur als Konventions-Quercheck berichtet,
um die Groesse des Normierungseffekts zu zeigen.

## 3. Sampler, Leiter, T-Gitter, Statistik (neu)

- **Sampler:** kanonischer Wolff-Single-Cluster (Numba-Kernel PHY045),
  1 Sweep = so viele Cluster, bis >= N Spins geflippt sind. Validierung:
  T->0-Orakel auf allen Gittern <= 0.05 %; gegen den Python-Wolff (PHY031)
  0.6 sigma (honeycomb L=12, T=0.60). tau_int(Upsilon) ~ 0.6 Sweeps (K4).
- **Leiter:** L in {32, 48, 64, 96, 128, 192, 256}. Das ergibt die
  HKS-Paare (L, 2L) = (32,64), (48,96), (64,128), (96,192), (128,256).
- **T-Gitter:** 0.545 .. 0.620, Schritt 0.005 (16 Punkte).
- **Statistik je (L, T):** 8 unabhaengige Seeds, je 500 Thermalisierungs-
  und 1500 Mess-Sweeps aus geordnetem Start (12000 Mess-Sweeps je Punkt;
  K4-Bedarf fuer SE(Upsilon_A) = 0.002 bei L <= 128: <= 8200).
  Seed-Vertrag: `47_000_000 + 1000*L + 10*t_idx + s`.
- **Fehler:** Jackknife ueber Seeds (leave-one-seed-out, fuer alle L
  gemeinsam) fuer jeden abgeleiteten Wert.
- **BLAS:** 1 Thread (Reproduzierbarkeit, AGENTS.md).

## 4. Schaetzer (neu)

- **Paar-Crossing T*(L, 2L):** C-eliminierte WM-Bedingung
  (`tbkt_pair_from_curves`, PHY040) auf den Upsilon_A-Seed-Mittelkurven.
  Lineare Interpolation nur zwischen benachbarten Gitterpunkten.
- **Primaer-Schaetzer T_W4:** HKS-Extrapolation
  T*(L) = T_c + a / ln^2(b L) ueber alle Paare mit Crossing (L = kleineres
  L des Paares); 3 Parameter (T_c, a, b), b in [0.1, 10].
- **FSS-Varianten fuer sigma_FSS** (vorab fixiert):
  - (v1) HKS mit b = 1 fest (2 Parameter);
  - (v2) HKS ohne das kleinste Paar;
  - (v3) groesstes Paar allein;
  - (v4) Weber-Minnhagen-Fit mit freiem C ueber alle L bei festem T
    (chi^2-Minimum ueber das T-Gitter).

## 5. sigma-Budget und Entscheidungsregeln

- **sigma_sampler** = Jackknife-SE von T_W4.
- **sigma_FSS** = halbe Spannweite ueber {T_W4, v1, v2, v3, v4}, mit
  Untergrenze **0.003**.
- **sigma_tot** = hypot(sigma_sampler, sigma_FSS).

Die Regeln aus v01 §7 gelten unveraendert, in dieser Reihenfolge:
- **S1:** weniger als 3 Paare mit Crossing -> NEGATIVE_RESULT.
- **S2:** sigma_tot > 0.010 -> NEGATIVE_RESULT.
- **C:** das Intervall [T_W4 - 2 sigma_tot, T_W4 + 2 sigma_tot] schneidet B.
- **I:** das Intervall ist disjunkt zu B.

Maschinell gilt `w4_verdict` (PHY044). **S0** (Budget) bezieht sich jetzt auf
die Rechenzeit: Ueberschreitet die Produktion 24 h Wall auf 4 Kernen, werden
die fehlenden L als nicht gemessen berichtet; es wird nicht umgeplant.

## 6. Ersetzt / bleibt

- **Ersetzt aus v01:** §3 (Observable pro Site, WL-Kernel), §4 (Walker-
  Domaenen, T_req als Budget-Kriterium), §5 (WL-Kalibrier-Budget als
  W4-Rezept), der Primaer-Schaetzer.
- **Bleibt:** Band B, die Regeln S1/S2/C/I, sigma_FSS-Floor 0.003, die
  Claim-Decke (kein Bestwert, Cross-Family-Review vor externer Aussage), der
  Blind-Grundsatz (Protokoll vor Daten).
- Die PHY044-WL-Kalibrierung bleibt gueltige Evidenz fuer den WL-Pfad; der
  WL-Pfad ist fuer W4 nur noch optionaler Quercheck.

## 7. Maschinelle Bindung

Die Konstanten `W4V2_*` im W4-Modul werden per Test gegen die Zahlen dieser
Spec geprueft (Leiter, T-Gitter, Seeds, Sweeps, Normierung, Seed-Vertrag).
