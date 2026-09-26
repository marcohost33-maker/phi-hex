# Audit O1 - Entscheid: Helicity-Normierung im Nelson-Kosterlitz-Kriterium

> **Status:** ENTSCHIEDEN am 2026-09-26 (PHY045). Ersetzt die Konventions-
> Entscheidung "per Site" des Wissenschafts-Audits 2026-06-04 (Befund S2)
> fuer jede NK-/Weber-Minnhagen-Auswertung auf nicht-quadratischen Gittern.
> Die Entscheidung stuetzt sich auf Theorie, einen exakten Gitter-Nachweis
> und einen MC-Nachweis ohne T_BKT-Lagewert; sie wurde VOR jeder
> Neuauswertung von T_BKT-Daten festgelegt und committet.

Coworker Research / Coworkerz, 2026-09-26.

## 1. Frage

Seit 2026-06-04 teilt `helicity_from_ensemble` durch die Site-Zahl N
(Upsilon_site), nicht durch die Flaeche A (Upsilon_area). Mit NN-Abstand 1
gilt Upsilon_site = a_s * Upsilon_area, a_s = A/N:

| Gitter | a_s | Upsilon_site(0) | Upsilon_area(0) |
|---|---:|---:|---:|
| square | 1 | 1 | 1 |
| triangular | sqrt(3)/2 = 0.8660 | 1.5 | sqrt(3) = 1.7321 |
| honeycomb | 3 sqrt(3)/4 = 1.2990 | 0.75 | 1/sqrt(3) = 0.5774 |
| kagome | 2/sqrt(3) = 1.1547 | 1 | sqrt(3)/2 = 0.8660 |

Welche der beiden Groessen springt am BKT-Uebergang auf 2 T_BKT / pi?

## 2. Theorie

Der universelle Sprung (Nelson & Kosterlitz, PRL 39, 1201 (1977)) ist fuer die
Steifigkeit K im Kontinuums-Funktional F = (K/2) int d^2x (grad theta)^2
formuliert - in NK als Flaechendichte angegeben. Die Twist-Definition
Upsilon = (1/A) d^2F/dk^2 (Fisher-Barber-Jasnow) liefert genau diese K; sie
ist invariant unter Umskalierung der Gitterkonstante. Dieselbe K bestimmt den
Tieftemperatur-Exponenten eta = T / (2 pi K); eta(T_BKT) = 1/4 ist aequivalent
zu K(T_BKT) = 2 T_BKT / pi. Die per-Site-Groesse ist nur auf Gittern mit
a_s = 1 (square) dieselbe.

## 3. Nachweis (PHY045, `results/260926 PHY045 helicity normalization O1 report.*`)

**A - exakt (MC-frei), harmonisches Modell auf den Repo-Gittern (L ~ 96):**
Die Steifigkeit, die <(theta_0 - theta_R)^2> ~ (T / (pi K)) ln R regiert, ist
auf allen vier Gittern Upsilon_area(0) auf <= 0.13 % (K/Upsilon_area:
1.0013 / 1.0004 / 1.0004 / 1.0007) und weicht von Upsilon_site(0) um genau
1/a_s ab (1.155 / 0.770 / 0.867).

**B - Monte Carlo, echtes XY-Modell, T ~ 0.5 T_BKT:** eta aus <m^2> ~ L^-eta
(L bis 128, Korrektur-Term 1/L^2) und Upsilon_site im selben Lauf;
R = 2 pi eta Upsilon_site / T ist 1 bei per-Site-, a_s bei per-Flaechen-
Korrektheit:

| Gitter | R | per Site (R=1) | per Flaeche (R=a_s) |
|---|---:|---:|---:|
| square (Kontrolle) | 0.989 +/- 0.003 | - | - |
| triangular | 0.865 +/- 0.005 | z = -30 | z = -0.25 |
| honeycomb | 1.289 +/- 0.005 | z = +54 | z = -1.8 |
| kagome | 1.152 +/- 0.007 | z = +22 | z = -0.4 |

Die square-Kontrolle zeigt eine Methoden-Systematik von ~1 %; innerhalb dieser
stimmt per Flaeche auf allen drei Gittern, per Site ist mit 22-54 sigma
falsifiziert.

**Literatur (nur Such-Snippets, Egress-Grenze der Arbeitsumgebung):**
- NK 1977 geben den Sprung als Flaechendichte an.
- Eine QMC-Gruppe normiert die Honeycomb-Steifigkeit ausdruecklich auf die
  Einheitszellen-Flaeche (Caci, Weber, Wessel, PRB 104, 155139 (2021)).
- Normierungsfreie Referenz triangular: Hochtemperatur-Reihe Butera & Pernici
  (arXiv:0806.1496) beta_c = 0.3412(4) in ihrer Konvention. Deren square-Wert
  0.5599(7) belegt den Faktor 2, also J/T_c = 0.6824(8), T_BKT = 1.465(1).
  Zusaetzlich nennt ein Snippet zu arXiv:1010.3075 (PRE 83, 011124) den Wert
  J_c = 0.6833(6), also T = 1.4635 - nur Snippet-Ebene.

## 4. Entscheid

1. **Das NK-Kriterium und die Weber-Minnhagen-Form verwenden Upsilon_area =
   Upsilon_site / a_s.** Upsilon_site bleibt als Rohgroesse und fuer die
   T=0-Orakel (Gate A) bestehen.
2. **Befund S2 (2026-06-04) ist falsifiziert.** Sein Argument - per Site
   treffe die Referenz 1.418, per Flaeche liege man ~15 % zu hoch - verglich
   stark driftende per-L-Crossings bei L = 9..19 mit einer Helicity-Referenz,
   deren eigene Normierung nicht belegt ist. Die normierungsfreie Referenz
   liegt bei ~1.465.
3. **Alle per-Site-T_BKT-Werte auf triangular, honeycomb und kagome sind
   konventions-verzerrt** (Vorzeichen wie beobachtet: triangular zu tief,
   honeycomb und kagome zu hoch). Sie bleiben als Lineage stehen und werden
   nicht still ersetzt. Die Neuauswertung mit Upsilon_area folgt als
   separates, datiertes Artefakt.
4. **Die triangular-Referenz 1.418** (Helicity-Crossing) gilt bis zum Beleg
   ihrer Normierung als konventions-abhaengig. Normierungsfreie Vergleichsbasis
   ist 1.465 (Reihe; MC-Snippet 1.4635).
5. **W4:** die Vorregistrierung v01 (per-Site-Upsilon) ist vor jeder W4-
   Datennahme durch v02 zu ersetzen
   (`spec/260926 PHI HEX w4 honeycomb preregistration v02.md`).

## 5. Offen

- Primaertext-Abgleich der Normierung in den zitierten Helicity-Arbeiten
  (Sorokin; Jiang; de Andrade et al.; Okabe & Otsuka).
- triangular per-Flaeche neu messen: das PHY030-T-Gitter endet bei 1.46 und
  enthaelt das per-Flaeche-Crossing nicht.
