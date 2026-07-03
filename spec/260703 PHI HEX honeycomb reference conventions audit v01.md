# Honeycomb Reference-Conventions-Audit fuer PHY041

Status: Vertragsquelle fuer PR #18 ab 2026-07-03.

## Regel

Honeycomb-Referenzen werden als Band mit Quelle, Observable und berichteter Variable gefuehrt. beta_BKT und T_BKT duerfen nicht vermischt werden.

Konversion:

```text
T = 1 / beta
sigma_T = sigma_beta / beta^2
```

## Referenzband

| Quelle | berichtete Groesse | T-Form fuer Vergleich | Kanal |
|---|---:|---:|---|
| arXiv:2501.07388 | T-Wert | 0.573 | Multi-Lattice |
| arXiv:2406.14812 | T_BKT,H = 0.571(8) | 0.571 +/- 0.008 | Helicity |
| arXiv:2406.14812 | T_BKT,H = 0.560(9) | 0.560 +/- 0.009 | NN/MC |
| arXiv:2406.12076 | beta_BKT = 1.687(3) | 0.5928 +/- 0.0011 | Upsilon |
| arXiv:2406.12076 | beta_BKT = 1.635(11) | 0.6116 +/- 0.0041 | Upsilon_4 |
| arXiv:2406.12076 | beta_BKT = 1.724(2) | 0.5800 +/- 0.0007 | Binder |
| Legacy-Anker | fruehere interne Kurznotation | 0.576 +/- 0.003 | nur noch Band-/Legacy-Kontext |

## Konsequenz

PHY041-L<=24 ist ein Pipeline-Finding. Die Paare 0.5951, 0.6029 und 0.6087 werden gegen dieses Referenzband eingeordnet. Daraus folgt kein neuer T_BKT-Bestwert.

## Naechste Stufe

1. Doku-Drift-Gate gruen halten.
2. L=32/48 mit demselben Kernel laufen lassen.
3. Upsilon_2- und Upsilon_4-FSS getrennt auswerten.
4. Danach gemeinsame Interpretation.
