# Reference-Conventions-Audit: Honeycomb `T_BKT` fuer PHY041

> **Status:** Vertrags-/Korrekturquelle fuer PR #18 / PHY041.  
> **Datum:** 2026-07-03.  
> **Scope:** Honeycomb-Referenzen fuer das 2D-XY-Modell, insbesondere die Unterscheidung `beta_BKT` vs. `T_BKT = 1/beta_BKT`, Observable-Abhaengigkeit und Repo-Konvention `J=1`, Helicity-Modulus per Site.

## 1. Warum dieses Audit noetig ist

PR #18 vergleicht die PHY041-Wang-Landau/Honeycomb-Paare gegen zwei externe Anker:

- `0.573` aus der Multi-Lattice-Linie `arXiv:2501.07388`.
- `0.576(3)` als dedizierte Honeycomb-Linie.

Die Web-/Quellenpruefung zeigt: Der dedizierte Honeycomb-Stand ist **kein einzelner kanonischer Skalar**. Mindestens zwei relevante 2024-Arbeiten berichten unterschiedliche Groessen/Observables:

1. `arXiv:2406.12076` berichtet primaer **inverse** Uebergangstemperaturen `beta_BKT` aus `Upsilon`, `Upsilon_4` und Binder-Crossing. Diese Werte muessen erst zu `T=1/beta` umgerechnet werden.
2. `arXiv:2406.14812` berichtet direkt `T_BKT,H`, u.a. `0.571(8)` aus dem Helicity-Modulus und `0.560(9)` aus Neural-Network-Auswertung.

Damit ist die fruehere Kurzform `0.576(3) (arXiv:2406.12076)` **zu weich / potentiell falsch attribuiert**. Bis die Original-PDF-Tabellen und Konventionsdetails vollstaendig nachgezogen sind, gilt: `0.576(3)` darf nicht als gesicherter Einzelwert aus `2406.12076` behauptet werden.

## 2. Referenztabelle

Alle Werte unten sind fuer Vergleichszwecke in `T_BKT` notiert. Wenn eine Quelle `beta_BKT` berichtet, wurde linear propagiert:

`T = 1 / beta`, `sigma_T = sigma_beta / beta^2`.

| Label im Repo | Quelle | Berichtete Groesse | Konvertiertes `T_BKT` | Observable / Methode | Status fuer PHY041 |
|---|---|---:|---:|---|---|
| `REF_MULTI = 0.573` | Okabe/Otsuka, `arXiv:2501.07388` | `T`-Wert aus Multi-Lattice-MC/ML-Uebersicht | `0.573` | XY/Clock auf verschiedenen 2D-Gittern | gueltiger externer Vergleichsanker, aber keine dedizierte Honeycomb-WL-Methodenkopie |
| `JIANG_HELICITY` | Jiang, `arXiv:2406.14812` | `T_BKT,H = 0.571(8)` | `0.571 +/- 0.008` | Helicity-Modulus | beste direkte `T`-Formulierung fuer Honeycomb-Helicity-Vergleich |
| `JIANG_NN` | Jiang, `arXiv:2406.14812` | `T_BKT,H = 0.560(9)` | `0.560 +/- 0.009` | Neural Network + MC | Zusatzanker, nicht direkt gleicher Estimator wie PHY041 |
| `AJD_Y2_BETA` | de Andrade/Jorge/DaSilva, `arXiv:2406.12076` | `beta_BKT = 1.687(3)` | `0.5928 +/- 0.0011` | zweiter Helicity-Modulus `Upsilon` | relevante WL/Helicity-Referenz, aber als `beta` berichtet und oberhalb 0.576 |
| `AJD_Y4_BETA` | de Andrade/Jorge/DaSilva, `arXiv:2406.12076` | `beta_BKT = 1.635(11)` | `0.6116 +/- 0.0041` | vierte Ordnung `Upsilon_4` | relevante WL/Y4-Referenz, klar observablespezifisch |
| `AJD_BINDER_BETA` | de Andrade/Jorge/DaSilva, `arXiv:2406.12076` | `beta_BKT = 1.724(2)` | `0.5800 +/- 0.0007` | Binder-Crossing | Zusatzanker, nicht derselbe Observable-Kanal wie PHY041-Y2-Paare |
| `LEGACY_DEDIC = 0.576(3)` | fruehere Repo-/Drive-Notation | `T`-Wert, genaue Herkunft in PR #18 nicht hart genug belegt | `0.576 +/- 0.003` | vermutlich dedizierte Honeycomb-Linie / Konventionsmix | **gelb:** nur noch als Legacy-/Band-Anker, nicht als hart attribuierter `2406.12076`-Wert |

## 3. Konsequenz fuer PR #18

PR #18 darf weiterhin sagen:

- PHY041 validiert eine **Pipeline**: WL-DOS, B&P-1/t, Auto-Fenster, Leak-Gate, Wolff-Quervalidierung, Drift-Guard gegen PHY032.
- Die kleinen-L-Paare `(12,16)/(12,24)/(16,24)` liegen bei `0.5951 / 0.6029 / 0.6087`.
- Das ist ein **FINDING**, kein neuer Bestwert.

PR #18 darf bis zur naechsten Haertung **nicht** mehr sagen:

- `0.576(3)` sei eindeutig der dedizierte `arXiv:2406.12076`-Honeycomb-WL-Wert.
- Die Abweichung gegen genau `0.576(3)` sei die einzige oder beste dedizierte Referenzdiagnose.

Korrekte Formulierung:

> Gegen das externe Honeycomb-Referenzband liegt PHY041 bei L<=24 im oberen Bereich: nahe an `arXiv:2406.12076`-Y2/Y4-konvertierten Werten, oberhalb der Jiang-/Okabe-Honeycomb-Anker um ca. 4-7 %. Daraus folgt kein neuer `T_BKT`-Bestwert; die Aussage bleibt finite-size-/Observable-/Konventions-sensitiv.

## 4. Numerische Guard-Rules fuer Tests

Diese Umrechnungen sind klein genug, um als Fast-Gate festgehalten zu werden:

```text
1 / 1.687 = 0.5928  (sigma ~ 0.0011)
1 / 1.635 = 0.6116  (sigma ~ 0.0041)
1 / 1.724 = 0.5800  (sigma ~ 0.0007)
```

Ein Fast-Test soll verhindern, dass `0.576(3)` wieder stillschweigend als direkter `2406.12076`-Y2-Wert in den Methodenvertrag rutscht.

## 5. Naechste wissenschaftliche Stufe

1. `README`/Report-Sprache auf **Referenzband** statt Einzelanker umstellen.
2. L=32/48 mit PHY041-Kernel laufen lassen.
3. Y2- und Y4-FSS getrennt berichten, danach erst gemeinsame Interpretation.
4. Bei Diskrepanz gegen Jiang/Okabe nicht sofort "Samplerfehler" behaupten: zuerst Konventionen, Observable, L-Leiter und `beta`/`T`-Achse pruefen.

## 6. Quellenanker

- de Andrade, Jorge, DaSilva: `arXiv:2406.12076`, Honeycomb-XY, Standard-MC + Wang-Landau, `beta_BKT` aus `Upsilon`, `Upsilon_4`, Binder.
- Jiang: `arXiv:2406.14812`, Honeycomb-XY, NN + MC, direkter `T_BKT,H`-Bericht.
- Okabe/Otsuka: `arXiv:2501.07388`, Multi-Lattice-XY/Clock-Vergleich, externer `0.573`-Anker.
- Hsieh/Kao/Sandvik: `arXiv:1302.2900`, logarithmische FSS-/Paar-Methodik fuer BKT.
- Belardinelli/Pereyra: `cond-mat/0702414` / JCP 127, 184105 (2007), `1/t`-WL gegen Fehler-Saettigung.
