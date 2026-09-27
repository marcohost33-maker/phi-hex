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

## Nachtrag 2026-09-26 (Issue #45 §3)

Die Zeile "Legacy-Anker 0.576 +/- 0.003" ist ueberholt: 0.576(3) ist der direkt
berichtete Wang-Landau-T-Wert von arXiv:2406.12076 (Abstract v4 / Physica
Scripta 100 065953) und steht jetzt als eigener T-Kanal `upsilon_wl_T` im Band.
Die Journal-Fassung von arXiv:2406.14812 (PTEP 2024 103A02) berichtet 0.576(4)
(Helicity) und 0.572(3) (NN). Die Regel dieser Spec (beta und T nie vermischen)
gilt unveraendert. Details, Fassungen und Beleg-Status:
`spec/260926 PHI HEX honeycomb reference provenance addendum v01.md`.


## Nachtrag 2026-09-27 - Versionsstatus supersediert

Diese Datei bleibt als historischer Vertragsstand vom 2026-07-03 erhalten.
Fuer **aktuelle** Literaturvergleiche supersediert jedoch
`spec/260926 PHI HEX honeycomb reference provenance addendum v01.md`
die oben stehende Banddarstellung.

Insbesondere gilt ab diesem Nachtrag:

- publizierte Jiang/PTEP-Kanaele: NN `0.572(3)`, Helicity `0.576(4)`;
- de-Andrade v4: direkte T-Kanaele `0.575(8)` (SA) und `0.576(3)` (WL);
- der Honeycomb-Geometriefaktor `4/(3 sqrt(3))` ist fuer Jiang/PTEP und
  de-Andrade v4 primaertextlich verifiziert;
- fruehe Jiang-v1- und de-Andrade-beta-Kanaele bleiben ausschliesslich
  **superseded_historical** Lineage und duerfen nicht als aktuelle
  Vergleichswerte in einen Akzeptanzbereich eingehen;
- Okabe/Otsuka `0.573` bleibt ein normierungsfreier Methoden-/Lageanker ohne
  erfundene Unsicherheit.

Neue Auswertungspfade muessen die versionierten Statusfelder in PHY042
(`REF_PROVENANCE`, `REF_CURRENT_KEYS`) oder den Provenienz-Nachtrag
verwenden, nicht die historische Tabelle dieser Datei isoliert.
