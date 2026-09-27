# Honeycomb-Referenzband: Provenienz-Nachtrag (Issue #45 §3)

Status: Nachtrag zur Vertragsquelle `spec/260703 PHI HEX honeycomb reference
conventions audit v01.md`, gueltig ab 2026-09-26. Die Regel aus 260703
("beta_BKT und T_BKT duerfen nicht vermischt werden") bleibt unveraendert und
wird hier verschaerft: jeder Band-Kanal fuehrt Quelle, **Fassung**, berichtete
Groesse und Beleg-Status.

Coworker Research / Coworkerz, 2026-09-26.

## 1. Anlass

Issue #45 §3 (Machbarkeitsstudie W4, 2026-09-11) meldete drei Provenienz-Luecken:

1. Der Zielwert `0.576(3)` stand nicht im Band (`REF_BAND`).
2. `binder_beta` 1.724(2) liess sich in arXiv:2406.12076 v4 per automatischer
   Extraktion nicht finden.
3. (implizit) Es war nicht festgehalten, **welche Fassung** einer Quelle ein
   Wert stammt.

## 2. Befund (Web-Recherche 2026-09-26)

**Beleg-Grenze vorab:** Die Egress-Policy dieser Umgebung blockt arxiv.org,
export.arxiv.org, iopscience, academic.oup.com, crossref und ADS. Alle Befunde
unten stammen aus Such-Snapshots der Abstract-Seiten (mindestens zwei
unabhaengige Suchlaeufe je Wert), **nicht** aus dem Primaertext. Status daher
`search_corroborated` (Confidence-Tier: UNBEKANNT -> STABIL erst nach
Primaertext-Abgleich, §4).

### 2.1 arXiv:2406.12076 (de Andrade, Jorge, DaSilva)

- Fassungen: v1 2024-06-17, v2 2024-11-27, v3 2024-12-06, v4 2025-04-11;
  Journal: Physica Scripta 100(6) 065953 (2025), DOI 10.1088/1402-4896/add62a.
- Das aktuelle Abstract berichtet **beide** Wertesorten nebeneinander:
  - beta_BKT = 1.687(3) (Upsilon), 1.635(11) (Upsilon_4), 1.724(2) (Binder);
  - T_BKT = 0.575(8) (simulated annealing) und **0.576(3) (Wang-Landau)**.
- Fassungs-Drift: Snapshots der v2 zeigen beta = 1.696(3) / 1.67(1) / 1.724(2);
  die Werte 1.687(3) / 1.635(11) gehoeren zu v3/v4. `binder_beta` 1.724(2) ist
  ueber v2..v4 im Abstract stabil -> Luecke 2 ist auf Abstract-Ebene geschlossen
  (Ort im Haupttext weiter offen).
- Upsilon_4-T-Werte 0.568(1) (WL) und 0.551(11) (SA): aus dem Haupttext (W4-
  Studie las v4; Such-Snapshots ordnen sie v1 zu) -> `version_unclear`.
- **Innere Spannung der Quelle:** 1/1.687 = 0.5928, aber die direkt berichtete
  WL-T ist 0.576 (entspraeche beta = 1.736). beta- und T-Werte stammen offenbar
  aus verschiedenen Auswertungen. Sie werden deshalb als **getrennte Kanaele**
  gefuehrt und nie ineinander umgerechnet.

### 2.2 arXiv:2406.14812 (F.-J. Jiang)

- arXiv v1 (2024-06-21): NN 0.560(9), Helicity 0.571(8) (= Vertrag 260703).
- Journal-Fassung PTEP 2024(10) 103A02 (DOI 10.1093/ptep/ptae147, 2024-09-27):
  **NN 0.572(3), Helicity 0.576(4)**.
- Die Journal-Fassung ersetzt die arXiv-v1-Werte **nicht still**: beide stehen
  getrennt im Band (`helicity_direct`/`nn_mc` = v1; `helicity_ptep`/`nn_ptep`).

### 2.3 arXiv:2501.07388 (Okabe, Otsuka)

- J. Phys. A 58 065003 (2025), DOI 10.1088/1751-8121/ada988. Honeycomb 0.573
  als grobe Schaetzung ohne Fehlerbalken (unveraendert).

## 3. Korrektur einer frueheren Einordnung

`spec/260703` fuehrte `0.576 +/- 0.003` als "Legacy-Anker, fruehere interne
Kurznotation", und `tests/test_honeycomb_reference_conventions.py` nannte ihn
`legacy_dedicated_anchor_unattributed`. Nach §2.1 ist 0.576(3) **der direkt
berichtete WL-T-Wert von arXiv:2406.12076** - die urspruengliche Attribution
(PHY035, PHY041: "0.576(3) arXiv:2406.12076") war richtig. Richtig bleibt auch
der Kern von 260703: 0.576 ist **nicht** die Umrechnung eines beta-Kanals dieser
Quelle; die beta-Kanaele liegen in T-Form bei 0.5928 / 0.6116 / 0.5800.

## 4. Referenzband (T-Form), Stand 2026-09-26

| Kanal (`REF_BAND`) | Quelle / Fassung | berichtet | T-Form | Status |
|---|---|---|---:|---|
| multi_lattice | arXiv:2501.07388 v1 = J.Phys.A 58 065003 | T | 0.573 | contract_260703 |
| helicity_direct | arXiv:2406.14812 v1 | T, Helicity | 0.571 +/- 0.008 | contract_260703 |
| nn_mc | arXiv:2406.14812 v1 | T, NN | 0.560 +/- 0.009 | contract_260703 |
| upsilon_beta | arXiv:2406.12076 v3/v4 | beta = 1.687(3) | 0.5928 +/- 0.0011 | contract_260703 |
| upsilon4_beta | arXiv:2406.12076 v3/v4 | beta = 1.635(11) | 0.6116 +/- 0.0041 | contract_260703 |
| binder_beta | arXiv:2406.12076 v2..v4 | beta = 1.724(2) | 0.5800 +/- 0.0007 | search_corroborated |
| upsilon_wl_T | arXiv:2406.12076 v4 / Phys.Scr. | T, Upsilon, WL | 0.576 +/- 0.003 | search_corroborated |
| upsilon_sa_T | arXiv:2406.12076 v4 / Phys.Scr. | T, Upsilon, SA | 0.575 +/- 0.008 | search_corroborated |
| upsilon4_wl_T | arXiv:2406.12076 Haupttext | T, Upsilon_4, WL | 0.568 +/- 0.001 | version_unclear |
| upsilon4_sa_T | arXiv:2406.12076 Haupttext | T, Upsilon_4, SA | 0.551 +/- 0.011 | version_unclear |
| helicity_ptep | PTEP 2024 103A02 | T, Helicity | 0.576 +/- 0.004 | search_corroborated |
| nn_ptep | PTEP 2024 103A02 | T, NN | 0.572 +/- 0.003 | search_corroborated |

Lesart:

- **Direkt berichtete T-Werte** (Upsilon_2-/NN-/Multi-Lattice-Kanaele) liegen
  in [0.560, 0.576]; die Journal-Fassungen allein in [0.572, 0.576].
- **beta-Kanaele** (in T-Form 0.580 .. 0.612) sind eine davon getrennte,
  quellen-interne Spannung. Sie werden gefuehrt, aber nicht zur Band-Definition
  eines Konsistenz-Ziels benutzt (siehe W4-Vorregistrierung, `spec/260926 PHI
  HEX w4 honeycomb preregistration v01.md`).

**Offene Pflicht (nicht in dieser Umgebung erfuellbar):** Primaertext-Abgleich
fuer jede `search_corroborated`/`version_unclear`-Zeile (Fassung, Ort: Abstract/
Tabelle/Abschnitt, verwendetes L, Upsilon-Normierung). Bis dahin gilt: kein
Kanal dieser Status-Klassen darf allein einen Claim tragen.

## 5. Maschinelle Bindung

- `REF_BAND` und `REF_PROVENANCE` in `src/260706 PHY042 honeycomb wl fss v01.py`
  (jeder Band-Kanal hat genau einen Provenienz-Eintrag).
- `tests/test_phy042_wl_fss.py` und `tests/test_honeycomb_reference_conventions.py`
  pruefen: jeder Kanal-Wert steht zeichengenau in 260703 **oder** in diesem
  Nachtrag; beta-Kanaele sind korrekt konvertiert; 0.576 erscheint nur in
  explizit als direkte T-Kanaele attribuierten Eintraegen.

## 6. Nachtrag 2 (2026-09-26, Selbstpruefung F9/F10)

**Zirkularitaets-Warnung:** das oeffentliche phi-hex-Repo ist selbst von
Suchmaschinen indexiert. Die Werte 1.687/1.635/0.5928 stehen im README und in
den Specs; Such-Zusammenfassungen koennen sie aus dem Repo statt aus der Quelle
beziehen. Eine zweite, unabhaengige Recherche (~70 Suchlaeufe) fand:

- arXiv:2406.12076 **v1** (Abstract): T_BKT = 0.576 +/- 0.001 (SA), 2mu = 4.0(5).
- **v2** (2024-11-27): beta_BKT = 1.696(3) (Upsilon), 1.67(1) (Upsilon_4),
  1.724(2) (Binder).
- **Journal-Fassung** (Phys. Scr. 100 065953): 0.575(8) (SA), 0.576(3) (WL),
  2mu = 5.80(12).
- Die Werte **1.687(3) / 1.635(11)** (Vertrag 260703) konnten **nicht
  unabhaengig bestaetigt** werden. Ihre Herkunft (v3/v4?) ist offen, eine
  Kontamination durch das Repo ist nicht ausgeschlossen. Status der Kanaele
  `upsilon_beta` / `upsilon4_beta` damit: **ungeprueft**. Sie tragen keinen
  Claim und sind nicht Teil von Band B.
- Auffaellig, aber unbelegt: 1/1.696 = 0.5896 liegt beim internen
  per-Site-Crossing (0.588-0.592). Ob beta-Kanaele der Quelle per-Site-
  normierte Crossings sind, ist offen (Normierung im Primaertext pruefen).

**Normierung (Audit O1 entschieden, per Flaeche):** fuer jeden Helicity-
basierten Literaturwert ist zu klaeren, ob per Site oder per Flaeche normiert
wurde. Normierungsfreie Kanaele (Correlation-Ratio 0.573, NN 0.572) sind davon
unberuehrt.

**triangular:** die normierungsfreie Hochtemperatur-Reihe (Butera & Pernici,
arXiv:0806.1496) gibt beta_c = 0.3412(4) in ihrer Konvention. Deren square-Wert
0.5599(7) = 1.1199/2 belegt den Faktor 2, also T_BKT = 1/0.6824 = 1.465(1).
Die bisherige Repo-Referenz 1.418 (Helicity-Crossing, Normierung unbelegt)
gilt als konventions-abhaengig.


## 7. Nachtrag 3 - Primaertext-Abgleich 2026-09-27

Die bis 2026-09-26 offene Normierungsfrage ist fuer den publizierten
Jiang-Kanal primaertextlich geklaert und fuer die beiden anderen zentralen
Referenzfamilien enger eingegrenzt.

### Jiang, PTEP 2024 103A02

Die publizierte Fassung definiert den Honeycomb-Helicity-Modulus mit dem
expliziten Vorfaktor `4/(3 sqrt(3))`. Der Text erklaert diesen Faktor als
Verhaeltnis der Spindichten von Honeycomb- und Quadratgitter. Genau dieser
korrigierte `Gamma(L)` wird anschliessend mit der Nelson-Kosterlitz-Linie
`2T/pi` verglichen. Der publizierte Helicity-Wert `T_BKT = 0.576(4)`
ist deshalb **nicht** als roher per-Site-Wert zu klassifizieren.

Status ab 2026-09-27:
- `helicity_ptep`: **primary_text_verified** fuer Wert, Wolff-Sampler und
  Normierungsfaktor;
- `nn_ptep`: **primary_text_verified** fuer den publizierten Wert
  `0.572(3)`.

Damit ist die bisherige Hypothese in
`results/260926 PHY046 W4 v02 interpretation note.md`, die Uebereinstimmung
der internen per-Site-Linie mit Jiang koenne durch eine per-Site-
Literaturnormierung erklaert werden, fuer die PTEP-Fassung **falsifiziert**.

### de Andrade, Jorge, DaSilva, arXiv:2406.12076 v4 / Phys. Scr. 100 065953

Die v4-Metadaten/Abstract-Evidenz bestaetigt die direkten T-Werte
`0.575(8)` (SA) und `0.576(3)` (WL). Ein belastbarer maschinenlesbarer
Primaertextbeleg fuer den konkreten Honeycomb-Normierungsfaktor wurde in
dieser Sitzung nicht reproduzierbar extrahiert. Daher bleibt der
Normierungsstatus dieser Kanaele **primary_text_pending**, statt aus der
Jiang-Konvention uebertragen zu werden.

Die historische beta-Provenienz `1.687(3)` / `1.635(11)` bleibt ebenfalls
offen und traegt weiterhin keinen Claim.

### Okabe/Otsuka, J. Phys. A 58 065003 / arXiv:2501.07388

Die Quelle beschreibt als Monte-Carlo-Methode explizit die Groessenabhaengigkeit
des Verhaeltnisses von Korrelationsfunktionen bei zwei verschiedenen
Distanzen. Der Honeycomb-Wert `0.573` wird als Schaetzung berichtet. Dieser
Kanal ist damit fuer die PHI-Hex-Frage besonders wertvoll, weil er keinen
Helicity-Normierungsfaktor verwendet.

Konsequenz fuer W4:
- die ~1.3-%-Spannung zwischen PHY046 (~0.565) und 0.572-0.573 bleibt eine
  **echte offene Diskrepanz**;
- sie darf nicht mehr durch die Jiang-per-Site-Hypothese weginterpretiert
  werden;
- der naechste diskriminierende Test ist der vorregistrierte
  Correlation-Ratio-Quercheck W4-v03/PHY049.
