# Phi-Hex

> **Status:** Forschungs-Repo (privat) | XY/BKT-Physik auf Dreiecks-, Honeycomb- und Kagome-Gittern.  
> **Lizenz:** Apache-2.0 | **Lineage/Provenance:** siehe `SOURCES.md`.  
> **Aktueller Review-Stand:** PHY042/PR #20 (Voll-Lauf L=24/32/48) ist ein Pipeline-Finding, **kein neuer T_BKT-Bestwert** (interne Gates PASS; Grenzen NR-PHY042-02/03).

Phi-Hex untersucht das 2D-XY-Modell und BKT-Physik auf periodischen Gittern: Helicity-Modulus, Nelson-Kosterlitz-Sprung, Wolff-Cluster, Wang-Landau-DOS und Finite-Size-Scaling.

## Konvention: Helicity-Modulus per Site

Seit dem Audit vom 2026-06-04 wird der Helicity-Modulus per Site normiert. Das ersetzt die fruehere Flaechen-Normierung, die T_BKT um etwa 15 % ueberschaetzt hatte.

Wichtige Orakel:

- Triangular: `Upsilon(0) = 1.5 J` per Site.
- Honeycomb: `Upsilon_2(0) = 3/4 J` per Site.
- Honeycomb ground state: `E0/N = -3J/2`.

## Honeycomb-Referenzen: Referenzband, kein Einzelanker

Ab PR #18/PHY041 gilt `spec/260703 PHI HEX honeycomb reference conventions audit v01.md` als Vertragsquelle. Honeycomb-Werte werden nicht mehr als ein einzelner harter Referenzwert gefuehrt, sondern als Band mit Quelle, Observable und beta/T-Konvention.

| Quelle | berichtete Groesse | T-Form fuer Vergleich | Rolle |
|---|---:|---:|---|
| arXiv:2501.07388 | T-Wert | 0.573 | Multi-Lattice-Anker |
| arXiv:2406.14812 | T_BKT,H = 0.571(8) | 0.571 +/- 0.008 | direkter Honeycomb-Helicity-Anker |
| arXiv:2406.14812 | T_BKT,H = 0.560(9) | 0.560 +/- 0.009 | NN/MC-Zusatzanker |
| arXiv:2406.12076 | beta_BKT = 1.687(3) | 0.5928 +/- 0.0011 | Upsilon / WL-Honeycomb |
| arXiv:2406.12076 | beta_BKT = 1.635(11) | 0.6116 +/- 0.0041 | Upsilon_4 |
| arXiv:2406.12076 | beta_BKT = 1.724(2) | 0.5800 +/- 0.0007 | Binder |

Konversion: `T = 1 / beta`, `sigma_T = sigma_beta / beta^2`.

## Mess-Stand

| Gitter | aktueller interner Stand | externe Einordnung | Methode |
|---|---:|---|---|
| square | 0.893 | Goldstandard/Referenzanker | Sandvik-Paar / V&V |
| triangular | 1.4007 +/- 0.0081 | nahe 1.418 | Wolff + Weber-Minnhagen |
| honeycomb | 0.5917 CI[0.587,0.598] fuer Paar (24,48) | oberhalb der unteren Honeycomb-Anker, finite-size-sensitiv | Wolff + Sandvik-Paar / PHY032 |
| kagome | 0.8479 CI[0.8414,0.8507] | nahe rough estimate 0.825 | Wolff + WM-Fit / PHY033 |

PHY041 liefert fuer honeycomb mit Wang-Landau/1-t bei L<=24 die glatten Paarwerte:

| Paar | T_BKT |
|---|---:|
| (12,16) | 0.5951 |
| (12,24) | 0.6029 |
| (16,24) | 0.6087 |

Interpretation: Das ist ein **FINDING** zur Pipeline und zum kleinen-L-Verhalten, kein finaler T_BKT. Der Wert liegt im oberen Bereich des Referenzbands bzw. oberhalb der unteren Honeycomb-Anker. Die naechste Stufe (L=32/48 plus getrennte Upsilon_2-/Upsilon_4-FSS) ist mit PHY042 gelaufen (siehe unten).

## PHY042 — Voll-Lauf L=24/32/48 (Finding, kein Bestwert)

PHY042 fuehrt den PHY041-Kernel auf L=24/32/48 mit getrennten Upsilon_2-/Upsilon_4-FSS-Kanaelen und Multi-Walker-Systematik aus (L=24: 1 Walker; L=32/48: je 3). Gate-Log: `results/260707 PHY042 honeycomb wl-fss L24-32-48 gate report.json` (`overall_pass=True`, 9/9 Pipeline-Gates, Lauf 2026-07-06, ~3498 s wall, `master_seed=42`, `lnf_final=1e-5`). Die Pipeline-Gates pruefen **Integritaet, nicht Physik-Wahrheit**; der Physik-Befund ist ein FINDING.

| Gate / Diagnostik | Wert |
|---|---|
| A: aligned `Upsilon_2(0) = 3/4 J` exakt | PASS |
| max. kanonisches Rand-Leak | 1.6e-10 |
| max. unbesetzte kanonische Masse in-Domaene | 1.9e-10 |
| Validitaets-Domaene `T_max` (Walker-Spread < 0.04) | L24: 0.67, L32: 0.60, L48: 0.585 |

Upsilon_2-Paar-Schaetzer (C-eliminiert, Mittel-Kurven), Einordnung **nur relativ zum Referenzband**:

| Paar | T_BKT | belastbar? |
|---|---:|---|
| (24,32) | 0.5875 | ja — Crossing in beiden Validitaets-Domaenen |
| (24,48) | 0.5898 | **nein** — NR-PHY042-02 |
| (32,48) | 0.5916 | **nein** — NR-PHY042-02 |

Ehrliche Grenzen:

- **Kein neuer Bestwert / kein finaler T_BKT-Claim.** L=48 bleibt endlich; die Walker-Systematik misst nur den Sampler (g(E)-Bias), nicht den finite-size-Bias. Ein neuer Bestwert braucht Cross-Family-Review.
- **NR-PHY042-02:** Die Paare (24,48) und (32,48) crossen ausserhalb der gemeinsamen Validitaets-Domaene und sind bei diesem Statistik-Budget nicht belastbar (Walker-Spread 0.033 bzw. 0.039 in T_BKT vs. 0.013 fuer das quotable Paar (24,32)).
- **NR-PHY042-03:** Der Upsilon_4-Dip ist bei L in {32,48} nicht walker-robust (L=24 robust bei T=0.65; L=32/48 walker-abhaengig). Upsilon_4-FSS oberhalb L=24 braucht mehr Produktion.

## Methodische Kernformeln

```text
H = -J * sum_<ij> cos(theta_i - theta_j)
T = 1 / beta
sigma_T = sigma_beta / beta^2
Upsilon(T_BKT) = 2 T_BKT / pi
Upsilon(T_BKT,L) = (2 T_BKT / pi) * (1 + 1 / (2 ln L + C))
A(E -> E') = min(1, g(E) / g(E'))
accept if log(u) <= log_g(E) - log_g(E')
w_T(E) = g(E) * exp(-E/T)
leak(T) = sum_edge w_T(E) / sum_all w_T(E)
```

## Reproduzieren

```bash
pytest -m "not slow"
pytest
python "src/260607 PHY032 honeycomb wm-logfit bootstrap v01.py"
python "src/260616 PHY040 wang-landau entropic helicity v01.py"
python "src/260702 PHY041 honeycomb wang-landau entropic helicity v01.py"
```

CI prueft Lint, Compile und schnelle Korrektheits-Gates. Slow-Messlaeufe bleiben lokal.

## Struktur

```text
src/        Engines + Experiment-Code
spec/       Specs, Theorie, Audits
results/    Gate-Logs, Reports, Negativ-Results
tests/      schnelle Gates + slow Mess-Smokes
archive/    Vorgaenger-Versionen
SOURCES.md  Provenance / SHA-256
```

## Naechste Stufe nach PR #20

1. Mehr Produktions-Statistik, damit die Crossings (24,48)/(32,48) in-Domaene und belastbar werden (NR-PHY042-02 aufloesen).
2. Walker-robuster Upsilon_4-Dip fuer L>=32 (NR-PHY042-03 aufloesen), dann Dip-FSS Richtung T_BKT.
3. Kein T_BKT-Bestwert-Claim ohne Cross-Family-Review.

---
*Coworker Research / Coworkerz | Repo-Anlage 2026-06-04 nach arbeitsschablone_forschungs-repo-anlage*
