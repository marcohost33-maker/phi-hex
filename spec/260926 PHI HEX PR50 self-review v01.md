# Selbstpruefung PR #50 (Issue #45): validiert, falsifiziert, korrigiert

Coworker Research / Coworkerz, 2026-09-26. Gegenstand: die drei Commits
ae96cd4 / 768fb19 / 001070b und die zugehoerige Berichterstattung. Methode:
jede Behauptung gegen Code, Daten, Tests, Primaerquellen (soweit erreichbar)
und unabhaengige Nachrechnung geprueft. Web-Recherche: nur Such-Snippets (die
Egress-Policy der Arbeitsumgebung sperrte arxiv/iop/oup/openalex/wayback).
Vorsicht: das oeffentliche phi-hex-Repo ist selbst indexiert, Snippets koennen
also zirkulaer sein.

## 1. Validiert (haelt)

| Behauptung | Pruefung |
|---|---|
| Einzel-Walker-Domaene war Nicht-Messung; Erratum reproduziert L32/L48 und alle Paar-Urteile | Tests + Erratum-Drift-Guard, CI gruen |
| Latenter Fail-open `min(0.6, nan)` in der Paar-Grenze | Regressionstest beider Reihenfolgen |
| Numba-WL-Kernel bit-identisch zum Python-Kernel (gleiche Maschine) | L=24-Produktionsjob, `results_identical` |
| PHY042-Trajektorien exakt reproduziert (7/7 wl_sweeps, Belegung) | VAL-BIT-P |
| np.dot bitweise thread-abhaengig ab n = 12288 | direkt gemessen (1/2/4 Threads) |
| 1/t-Phase greift bei L >= 48 mit lnf 1e-5 nicht | gemessen; in der Literatur bekannt (arXiv:2412.00809) |
| Walker-Spread produktions-dominiert | K3 qualitativ robust (lnf 1e-6 wirkungslos, prod 4x wirksam) |
| Vorregistrierung vor den Kalibrierlaeufen committet | Git-Historie 768fb19 vor 001070b |
| CI gruen auf py3.12/3.13/3.14 inkl. kompiliertem Numba-Gate | GitHub Checks |

## 2. Falsifiziert oder korrigiert

**F1 (kritisch): W4 v01 misst mit der falschen Normierung.**
- v01 legte Upsilon pro Site als Primaer-Observable fest und verglich deren
  2T/pi-Crossings mit Literaturwerten.
- Audit O1 ist jetzt entschieden (PHY045: exakt + MC, per Site mit 22-54 sigma
  falsifiziert): das NK-Kriterium verlangt Upsilon pro Flaeche.
- v01 haette damit einen um den Faktor 1.299 verzerrten Schaetzer gegen ein
  Literaturband getestet.
- Korrektur: W4 v02 vor jeder Datennahme.

**F2: Methodenwahl nicht kalibriert.**
- PHY044 verglich nur WL-Budget-Knoepfe, nicht WL gegen kanonischen Wolff.
- K4 (blind) zeigt: Wolff erreicht SE 0.002 je T-Punkt bei L = 128 in ~43 s.
  WL braucht fuer dieselbe Genauigkeit bei L = 64 ~6 CPU-h.
- Korrektur: v02 wechselt den Primaer-Sampler.

**F3: Primaer-Schaetzer ohne FSS-Extrapolation.**
- v01 nahm das groesste Paar ohne Extrapolation.
- Best Practice ist die HKS-Paar-Extrapolation mit ln^2-Korrektur.
- Korrektur: v02 §4.

**F4: Begruendung des Bandes B inkonsistent.**
- v01 nahm den Binder-beta-Kanal als obere Grenze und schloss die anderen
  beta-Kanaele aus.
- Korrektur: v02 §1 (B unveraendert, Begruendung ehrlich).

**F5: "PHY042 in 231 s statt ~3500 s" (README/CHANGELOG/PR) war irrefuehrend.**
- Der Vergleich ueberspannt zwei Maschinen.
- Der Speedup auf gleicher Maschine ist 5.6x.
- Korrektur im README.

**F6: "94 % produktions-dominiert" ist ueberpraezise.**
- Der Wert ist ein Quotient zweier Spannweiten aus je 3 Stichproben.
- Qualitativ haelt die Aussage, die Zahl traegt eine grosse Unsicherheit.

**F7: W4-GO (v01) war fragil.**
- T_max = T_req ohne Marge, eine Realisierung, Gitter-Quantisierung 0.005.
- Durch v02 gegenstandslos.

**F8: Blind-Leck (gering).**
- Der PHY044-Report zaehlt Walker-Kombinationen "mit Crossing im Fenster".
- Das verraet, ob ein Crossing in [0.52, 0.67] existiert.
- Dokumentiert; kein Einfluss auf v02 (andere Observable, anderer Sampler).

**F9: Provenienz-Nachtrag zu optimistisch.**
- Die beta-Werte 1.687(3)/1.635(11) waren als "v3/v4" gefuehrt.
- Eine zweite Recherche konnte sie nicht bestaetigen und fand eine
  Kontaminationsgefahr durch das eigene Repo.
- Belegt sind v1: T = 0.576(1) (SA); v2: beta 1.696(3)/1.67(1)/1.724(2);
  Journal: 0.575(8)/0.576(3).
- Korrektur: Provenienz-Nachtrag 2.

**F10: Triangular-Referenz 1.418 ungeprueft.**
- Der Wert stammt aus einem Helicity-Crossing mit unbelegter Normierung.
- Die normierungsfreie Hochtemperatur-Reihe gibt 1.465 (Butera & Pernici).
- Konsequenz: siehe O1-Entscheid.

## 3. Staerken (bleiben)

- Fail-closed-Denken (Walker-Plan, Paar-Grenze, Blind-Assert).
- Lineage-Disziplin: Errata statt Umschreiben, SHA-Pins.
- Vorregistrierung vor den Daten, mit Git-Zeitstempel belegt.
- Bit-Identitaets-Gates fuer Beschleuniger.
- Offene Benennung von Grenzen (search_corroborated, Egress).

Gerade diese Disziplin hat F1 bis F10 vor jeder W4-Datennahme auffindbar
gemacht.

## 4. Konsequenzen

1. `spec/260926 PHI HEX O1 helicity normalization decision v01.md` (Entscheid
   per Flaeche; Befund S2 vom 2026-06-04 falsifiziert).
2. `spec/260926 PHI HEX w4 honeycomb preregistration v02.md` (Upsilon pro
   Flaeche, Wolff, HKS; vor jeder W4-Datennahme).
3. PHY045: Normierungs-Nachweis und Wolff-Effizienz als Gate-Log.
4. Neuauswertung der committeten honeycomb/kagome-Daten pro Flaeche als
   separates, datiertes Artefakt, NACH den Punkten 1 und 2.
5. README/CHANGELOG/Provenienz-Nachtrag 2; Core-Docstring zur Normierung
   korrigiert (Lineage-Vermerk, keine stille Aenderung der Rohgroesse).
