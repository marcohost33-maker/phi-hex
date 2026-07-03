# Methoden-Spec: PHY041 - honeycomb entropic T_BKT

> **Status:** SUPERSEDED-IN-PLACE fuer Referenzkonventionen.  
> Diese Datei bleibt als historische PHY041-Methoden-Spec erhalten, aber alle Honeycomb-Referenzvergleiche sind seit 2026-07-03 durch `spec/260703 PHI HEX honeycomb reference conventions audit v01.md` zu lesen.

## Aktueller Vertragsstand

PHY041 validiert den honeycomb-Wang-Landau/1-t-Apparat fuer kleine Gitter L=12/16/24. Das Ergebnis ist ein Pipeline-Finding, kein neuer unendlicher-L-Bestwert.

Die externe Honeycomb-Referenzlage wird als Band gefuehrt:

- arXiv:2501.07388: Multi-Lattice-Anker T=0.573.
- arXiv:2406.14812: direkter Honeycomb-Helicity-Anker T=0.571 +/- 0.008; NN-Zusatzanker T=0.560 +/- 0.009.
- arXiv:2406.12076: berichtet beta-Werte. Fuer den T-Vergleich gilt `T = 1 / beta`, `sigma_T = sigma_beta / beta^2`; daraus folgen ca. 0.5928, 0.6116 und 0.5800 fuer die dortigen Kanaele.

## Verbotene Lesart

Die alte Kurzform eines einzelnen dedizierten Honeycomb-Ankers ist nicht mehr mergefaehig. Insbesondere darf kein Wert aus beta-Angaben ohne beta/T-Konversion als direkter T-Wert behandelt werden.

## Methode, weiterhin gueltig

- Gitter-agnostischer Kernel auf Adjazenz, Kanten und Bond-Projektionen.
- Belardinelli-Pereyra-1-t-Phase gegen klassische Wang-Landau-Fehler-Saettigung.
- Auto-Energiefenster aus Wolff-Ankern.
- Edge-Leak-Gate fuer jedes Analyse-T.
- Wolff-Quervalidierung und Drift-Guard gegen committed PHY032-Daten.

## PHY041-Zahlen

| Paar | T_BKT |
|---|---:|
| (12,16) | 0.5951 |
| (12,24) | 0.6029 |
| (16,24) | 0.6087 |

Interpretation: kleine-L-Finding im oberen Referenzband; kein Bestwert.

## Naechste Stufe

1. PR #18 bleibt Draft, bis CI und Doku-Drift-Gate gruen sind.
2. L=32/48 mit demselben Kernel.
3. Upsilon_2- und Upsilon_4-FSS getrennt.
4. Mehrere g(E)-Walker nur als Systematik-Check.
