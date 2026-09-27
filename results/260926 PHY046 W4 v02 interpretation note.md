# W4 v02 - Ergebnis und Einordnung (honeycomb T_BKT, Upsilon pro Flaeche)

Coworker Research / Coworkerz, 2026-09-26. Gate-Logs:
`results/260926 PHY046 honeycomb w4 wolff area-helicity report.*`
(vorregistriert), `results/260926 PHY048 square pipeline validation
report.*` (post-hoc). Die Vorregistrierung v02 wurde in 82e2738 committet,
vor jeder W4-Datennahme.

## 1. Vorregistriertes Ergebnis (PHY046, unveraendert gemaess Spec v02)

| Paar (L, 2L) | T*(L, 2L), Upsilon pro Flaeche | Jackknife-SE |
|---|---:|---:|
| (32, 64) | 0.5658 | 0.0004 |
| (48, 96) | 0.5655 | 0.0009 |
| (64, 128) | 0.5657 | 0.0003 |
| (96, 192) | 0.5649 | 0.0010 |
| (128, 256) | 0.5663 | 0.0020 |

- **T_W4** (HKS, 3 Parameter) = **0.5649**.
  - Varianten: v1 0.5652, v2 0.5646, v3 0.5663, v4 0.5673.
- **Unsicherheit:**
  - sigma_sampler = 0.0026;
  - sigma_FSS = 0.0014, auf die Untergrenze 0.003 angehoben;
  - sigma_tot = **0.0040**.
- **Verdikt (Regel C): CONSISTENT** mit B = [0.560, 0.580]. Das Intervall
  +/- 2 sigma_tot ist [0.557, 0.573].
- **Diagnose-Flag:** der HKS-Fit liegt mit b = 10 am Rand des vorregistrierten
  Bereichs. Die Paare sind flach (keine L-Drift aufloesbar), die
  ln^2-Korrektur ist daher nicht bestimmbar. Die Varianten stimmen auf 0.0027
  ueberein.
- Per Site (nur Quercheck): HKS 0.5713; die Paare fallen von 0.5946 (32, 64)
  auf 0.5865 (128, 256).

## 2. Post-hoc-V&V derselben Pipeline auf dem Quadratgitter (PHY048, NICHT vorregistriert)

Referenz T_BKT = 0.8929-0.8935. Paar-Crossings bei L = 32..256: 0.8919,
0.8916, 0.8888, 0.8935, 0.8933. Das ist -0.46 % bis +0.06 %.

- **Die Paar-Crossings sind bei diesen L auf <= 0.5 % unverzerrt.** Ein
  systematischer Tiefbias von ~1 %, wie ihn PHY028 bei (16, 32) zeigte, ist
  bei L >= 32 nicht zu sehen.
- **Der 3-Parameter-HKS-Fit ist rauschverstaerkend:** 0.8969 +/- 0.0093, die
  Varianten reichen bis 0.905. Mit diesen Regeln waere square an S2
  gescheitert (sigma_tot ~ 0.011). Robust sind das groesste Paar (0.8933) und
  der WM-Fit mit freiem C (0.8922). Auf honeycomb war der Fit nur deshalb
  stabil, weil die Paare flach und praeziser sind.
- **Lehre fuer kuenftige Protokolle** (v03, vor neuen Daten): als
  Primaer-Schaetzer das gewichtete Mittel der Paare bzw. das groesste Paar
  nehmen, HKS nur als Variante.

## 3. Einordnung (ehrlich)

- Honeycomb pro Flaeche liegt bei **~0.565 +/- 0.004**.
  - Die direkt berichteten Literaturwerte liegen hoeher: Journal 0.572-0.576;
    normierungsfrei Correlation-Ratio 0.573 (grob, ohne Fehlerbalken) und
    NN 0.572(3).
  - Die Spannung zu den normierungsfreien Werten betraegt ~1.3 % (~2 sigma).
    Die square-V&V erklaert sie nicht.
- **Primaertext-Korrektur 2026-09-27:** die fruehere Hypothese, die
  Uebereinstimmung der internen per-Site-Linie mit Jiang koenne auf einer
  per-Site-Literaturnormierung beruhen, ist fuer die publizierte PTEP-Fassung
  falsifiziert. Jiang definiert Gamma mit dem Honeycomb-Dichtefaktor
  4/(3 sqrt(3)) und vergleicht diesen korrigierten Modulus mit 2T/pi.
  - Der publizierte Helicity-Wert 0.576(4) ist daher kein roher per-Site-Wert.
  - Der normierungsfreie NN-Wert 0.572(3) bleibt ebenfalls oberhalb PHY046.
  - Primaertext-Nachtrag 2026-09-27: auch de Andrade v4 multipliziert
    Upsilon/Upsilon_4 fuer Honeycomb explizit mit 4/(3 sqrt(3)) und verwendet
    danach das Nelson-Kosterlitz-Kriterium. Die aktuelle direkte T-Linie ist
    damit normierungsseitig geklaert; historische beta-Linien bleiben
    versionsspezifische Lineage.
- **Kandidaten fuer die Rest-Spannung:**
  - (a) honeycomb-spezifische Finite-Size-Korrekturen, groesser als auf square;
  - (b) Biases der normierungsfreien Literaturwerte (grobe Schaetzung; NN);
  - (c) Torus-Form-/Windungs-Effekte auf den Sprung. Fuer den Rhombus-Torus
    werden sie als klein erwartet, sind aber nicht nachgerechnet.
- **Naechster diskriminierender Schritt:** W4-v03/PHY049 ist am 2026-09-27
  vorregistriert: eigener normierungsfreier Correlation-Ratio-Schaetzer
  R=<g(L/2)>/<g(L/4)> auf L={48,72,96,144,192}. Die Mess-Engine ist
  fail-closed. W4-v03a/PHY050 implementiert inzwischen das deterministische
  FSS-Recovery-Gate mit positiven und adversarial synthetischen Kontrollen.
  Ein gruener G4 autorisiert nur die Produktion; Physikinterpretation bleibt
  bis G0/G5/G6 auf echten PHY049-Daten gesperrt.

## 4. Claim-Decke

FINDING relativ zum vorregistrierten Band. Kein Bestwert, keine Aussage
0.573 vs 0.576. Vor jeder externen Aussage braucht es den Primaertext-Abgleich
der Literatur-Normierungen und ein Cross-Family-Review.
