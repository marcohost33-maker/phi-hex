# Honeycomb Reference-Conventions-Audit fuer PHY041

Status: Vertragsquelle fuer PR #18 ab 2026-07-03.  
Nachtrag 2026-09-27 / Issue #45: Quellenrevisionen und Publikationsfassungen
neu abgeglichen; aktuelle und supersedierte Werte werden strikt getrennt.

## Regel

Honeycomb-Referenzen werden mit **Quelle, Version/Publikationsstand, Observable,
berichteter Variable und Unsicherheit** gefuehrt. Ein historischer arXiv-Wert
darf nicht still als aktueller Literaturanker weiterleben. beta_BKT und T_BKT
duerfen nicht vermischt werden.

Konversion fuer historische beta-Angaben:

```text
T = 1 / beta
sigma_T = sigma_beta / beta^2
```

## Aktueller Evidenz-Ledger (Stand 2026-09-27)

Der Ledger ist **kein** durch Min/Max konstruiertes Akzeptanzintervall. W4
vergleicht seine Schaetzung gegen die einzelnen Quellen/Kanaele und behandelt
Methodenstreuung explizit.

| Key | Quelle / Stand | berichtete Groesse | T-Form fuer Vergleich | Rolle |
|---|---|---:|---:|---|
| okabe_otsuka_2501_v1_rough | arXiv:2501.07388v1 / J. Phys. A 2025 | rough T_BKT | 0.573 | Multi-Lattice, keine zitierte Unsicherheit |
| jiang_ptep_nn | PTEP 2024, 103A02 (publizierte Fassung) | T_BKT,H (NN) | 0.572 +/- 0.003 | direkter Honeycomb-Zusatzanker |
| jiang_ptep_helicity | PTEP 2024, 103A02 (publizierte Fassung) | T_BKT,H (Helicity) | 0.576 +/- 0.004 | direkter Honeycomb-Helicity-Anker |
| andrade_v4_helicity_sa | arXiv:2406.12076v4 | T_BKT (Upsilon, SA) | 0.575 +/- 0.008 | Y2 / aktuelle Fassung |
| andrade_v4_helicity_wl | arXiv:2406.12076v4 | T_BKT (Upsilon, WL) | 0.576 +/- 0.003 | Y2 / aktuelle Fassung |
| andrade_v4_upsilon4_sa | arXiv:2406.12076v4 | T_BKT (Upsilon_4, SA) | 0.551 +/- 0.011 | Y4-Diagnostik |
| andrade_v4_upsilon4_wl | arXiv:2406.12076v4 | T_BKT (Upsilon_4, WL) | 0.568 +/- 0.001 | Y4-Diagnostik |

**Primaerset fuer W4:** Okabe/Otsuka 0.573, Jiang 0.572(3)/0.576(4) und
de Andrade/Jorge/DaSilva v4 Y2 0.575(8)/0.576(3). Die Y4-Werte bleiben
Diagnostik, weil dieselbe Quelle die Minimum-/FSS-Systematik diskutiert.

## Superseded / historische Werte (nicht als aktueller W4-Anker)

Diese Werte bleiben fuer Reproduzierbarkeit alter Reports erhalten:

| Quelle / alte Fassung | alter Wert | historische T-Form | Status |
|---|---:|---:|---|
| arXiv:2406.14812v1 | T_BKT,H = 0.571(8) | 0.571 +/- 0.008 | **superseded** durch publizierte PTEP-Fassung |
| arXiv:2406.14812v1 | T_BKT,H = 0.560(9) | 0.560 +/- 0.009 | **superseded** durch publizierte PTEP-Fassung |
| fruehe 2406.12076-Fassung | beta_BKT = 1.687(3) | 0.5928 +/- 0.0011 | **superseded** fuer aktuellen Vergleich |
| fruehe 2406.12076-Fassung | beta_BKT = 1.635(11) | 0.6116 +/- 0.0041 | **superseded** fuer aktuellen Vergleich |
| fruehe 2406.12076-Fassung | beta_BKT = 1.724(2) | 0.5800 +/- 0.0007 | **superseded**; Binder-Ergebnis nicht in v4 verifiziert |
| Legacy-Anker | fruehere interne Kurznotation | 0.576 +/- 0.003 | nur Lineage |

## Konsequenz fuer bestehende Laeufe

PHY041/PHY042 sind datierte Pipeline-Findings. Ihre historischen Reports werden
nicht umgeschrieben. Neue Auswertungen duerfen die supersedierten Werte jedoch
nicht als aktuelle Literaturreferenzen verwenden.

Fuer den committed PHY042-Report gilt zusaetzlich die Erratum-Semantik aus
Issue #45: `domain_tmax_spread004["24"] = 0.67` entstand aus einem
Ein-Walker-Sonderpfad und ist **keine gemessene Walker-Spread-Domaene**.
Aktueller Code serialisiert einen solchen Fall als `null` und haelt einen
separaten Evidenzgrund fuer alternative Guards fest.

## Naechste Stufe

1. W4 vorregistriert und fail-closed ausfuehren; mindestens 3 Walker je
   verwendeter L-Groesse.
2. Vor Grossproduktion zuerst L=64-Kalibrierlauf fuer Laufzeit,
   WL-Konvergenz und Sampler-Domaene.
3. Y2 und Y4 getrennt auswerten; Y4 nicht als stillen Y2-Ersatz verwenden.
4. Kein neuer T_BKT-Bestwert ohne Modellwahl-/Finite-Size-Systematik und
   Cross-Family-Review.
