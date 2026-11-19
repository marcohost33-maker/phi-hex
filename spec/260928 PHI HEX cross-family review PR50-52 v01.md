# 260928 PHI HEX Cross-Family-Review PR #50 / #51 / #52 + Literatur-Provenienz v01

Coworker Research / Coworkerz, 2026-09-28. Unabhaengige Review-Session
(andere Session als PR #50; PR #52 stammt aus der Codex-Familie). Methode:
Code, CI-Logs, Git-Historie und Repo-Daten direkt geprueft; Literatur nur
ueber Web-Such-Snippets (Egress-Grenze, Abschnitt 3.0). Jede Zahl unten traegt
ihre Quelle; nichts wurde still umgeschrieben.

Gegenstand des Reviews ist der offene PR-Stapel zu Issue #45 (W4) und die
darin getroffene Entscheidung Audit O1 (Upsilon pro Flaeche). Der PR #54
(PHY051) setzt den daraus folgenden, noch niemandem zugeordneten Messschritt
um; sein Ergebnis wird NICHT in diesem Dokument interpretiert, sondern im
PHY051-Report und README (Trennung von Review und Messung).

## 1. Stand des PR-Stapels (Snapshot 2026-09-28, UTC)

| PR | Branch (Head) | Base | CI | Befund |
|---|---|---|---|---|
| #50 (Draft) | `claude/phihex-projekt-next-steps-trt4me` (`21b51d1`) | main | **gruen** (Run 36281005868: ruff, py3.12/3.13/3.14, CodeQL) | Integritaets-Reparatur Issue #45, W4-Vorregistrierung v01/v02, PHY044-048, O1-Entscheid. Lokal reproduziert: ruff clean, 213 Tests PASS. |
| #51 (offen, kein Draft) | `fix/issue-45-w4-evidence-contract` (`f1c3ced`) | main | **rot** (alle 3 Python) | Einziger Fehler: `test_sources_integrity` - `src/260927 PHY044 honeycomb w4 calibration v01.py` ohne passende SOURCES-Zeile. Autor-Kommentar 2026-09-27: nicht separat mergen; Inhalte selektiv nach #52 portiert. |
| #52 (Draft, agent:codex) | `codex/w4-v03-correlation-ratio` (`c92d322`, 2026-09-28 01:55 +0200) | #50-Branch | **rot** (alle 3 Python), lint + zizmor gruen | Einziger Fehler: `test_sources_integrity` - `tests/test_phy050_correlation_fss.py` (Abschnitt 2). 294/295 Tests PASS. |
| #43/#46/#48/#53 (Dependabot) | scipy>=1.18.1, zizmor 0.6.3, numpy>=2.5.3, ruff>=0.16.8 | main | - | unabhaengig; Floor-Bumps, keine Code-Beruehrung. |
| #54 (Draft, agent:claude) | `claude/phi-hex-next-steps-21z9lw` | #50-Branch | siehe PR | PHY051 (vorregistriert) + dieses Review. |

**Ueberlappung #51 vs #52** (beide gegen ihre Base): 10 gemeinsame Dateien
(README, SOURCES, `spec/260703`, `spec/260808`, PHY041/042/043 + deren Tests).
#51 und #50 fuehren zudem je ein eigenes PHY044-Modul mit verschiedenem
Dateinamen (`260927 ...` vs `260926 ...`) und verschiedener Semantik. Ein
Merge von #51 nach #50/#52 waere ein Doppel-PHY044 mit Konflikten in allen
zehn Dateien.

**Wichtig fuer gestapelte PRs:** `ci.yml` auf main triggert `pull_request`
nur mit Base `main`. #52 und #54 bekaemen ohne Aenderung KEINE CI. #52 hat
das im eigenen `ci.yml` korrigiert; #54 traegt denselben Hunk byte-identisch
(Commit `96be3cc`), damit beide Reihenfolgen konfliktfrei mergen.

### Empfohlene Merge-Reihenfolge (AGENTS.md: sequenziell, Auto-Merge)

1. Dependabot #43/#46/#48/#53 nach main (unabhaengig, harmlos).
2. **#50 -> main.** Gruen, in sich geschlossen; Draft -> ready nach Sichtung
   des O1-Entscheids (Abschnitt 4 dieses Reviews stuetzt ihn). Offene
   Checklisten-Punkte des PR-Bodys sind durch #52 (v03, Primaertext-Teil)
   und #54 (triangular/kagome) adressiert, nicht durch #50 selbst.
3. **#52** nach Merge von #50 (GitHub retargetet die Base automatisch auf
   main, wenn der Base-Branch nach dem Merge geloescht wird; sonst Base
   manuell auf main setzen): fehlende SOURCES-Zeile (Abschnitt 2) nachziehen,
   CI gruen abwarten, dann mergen.
4. **#54** analog; Konfliktflaeche zu #52: nur append-only `SOURCES.md`,
   `README.md`- und `CHANGELOG.md`-Abschnitte (rebase trivial).
   Squash-Merge-Hinweis (AGENTS.md: `--squash`, "zweiter paralleler PR
   rebased"): nach dem Squash von #50 enthalten #52 und #54 dessen
   Original-Commits weiter als Vorfahren; beide sind dann mit
   `git rebase --onto main 21b51d1 <branch>` auf main zu setzen
   (kein Force-Push auf fremde Branches; auf eigenen Agent-Branches ist
   `--force-with-lease` nach Rebase die Konvention).
5. **#51 schliessen** (nicht mergen) mit Verweis auf den Autor-Kommentar
   vom 2026-09-27; die noch nicht portierten Teile (PHY044 v01-Kalibrier-
   Runner, PHY041-1/t-Metadaten) sind in #52 als Lineage vermerkt.
   Entscheid liegt beim Autor.

## 2. CI-Rot-Ursachen (exakt, lokal reproduziert)

**#52** (Run 36360323565, Job 108736083312, py3.12; identisch auf 3.13/3.14):

```
FAILED tests/test_sources_integrity.py::test_every_scoped_file_has_matching_sources_row
  Dateien ohne byte-genau passende, pfadgebundene SOURCES.md-Zeile:
  tests/test_phy050_correlation_fss.py
1 failed, 294 passed, 4 deselected in 303.88s
```

Ursache: die Datei wurde nach der letzten Pin-Zeile (`4235FD75CE39C3C7`)
in `72987f9`/`c92d322` erneut geaendert; der committete Blob hat den
SHA-256-Praefix `0B6282E86E53462A`. Fix = eine append-only Zeile, lokal
verifiziert (Worktree auf `c92d322`, danach `test_sources_integrity` 4/4
PASS):

```
| 0B6282E86E53462A | 2026-09-28 | tests | repo-native: `tests/test_phy050_correlation_fss.py` (stand: G4-Gate-Log versionsuebergreifend semantisch verglichen - Hash nach letzter Test-Aenderung nachgepinnt; supersediert vorherige Zeile fuer diesen Pfad) |
```

Hinweis zur Klasse des Fehlers: #51 und #52 sind an derselben Stelle rot
(SHA-First-Gate nach Nach-Edit). Das Gate arbeitet wie vorgesehen; die
Lehre fuer Agenten-Workflows ist "SOURCES-Zeile im selben Commit wie die
letzte Aenderung", nicht ein weicheres Gate.

**#51** (Run 36302575984): gleiche Klasse, `src/260927 PHY044 honeycomb w4
calibration v01.py` ohne Zeile (1 failed, 161 passed). Wegen Empfehlung
"schliessen" kein Fix noetig.

## 3. Literatur-Provenienz (Stand 2026-09-28)

### 3.0 Beleg-Grenze

Die Egress-Policy dieser Umgebung blockt (am 2026-09-28 einzeln verifiziert):
arxiv.org, academic.oup.com, iopscience.iop.org, inspirehep.net,
alphaxiv.org, ui.adsabs.harvard.edu, semanticscholar.org, huggingface.co,
pubmed.ncbi.nlm.nih.gov, researchgate.net. Alle Befunde unten sind
Such-Snippets (`search_corroborated`), keine Primaertexte - derselbe Status
wie in PR #50/#52. Zirkularitaets-Warnung aus #50 gilt: das Repo ist
indexiert. Eine Drive-Suche (Volltext nach den arXiv-IDs, PDF-Typ) fand
KEINE abgelegten Primaertexte. **Konkreter Weg zum Abschluss:** die fuenf
PDFs (2406.12076v4, PTEP 103A02 / 2406.14812, 2501.07388, 2305.00651,
0806.1496) lokal laden und in `007acc2/00_EINGANG/` ablegen; Drive-PDFs
sind ueber die Drive-API lesbar, dann ist der Abgleich in einer Session
erledigt.

### 3.1 Befunde je Referenz

| Referenz | Was es ist | Wert(e) | Normierungs-Status |
|---|---|---|---|
| arXiv:2501.07388 = Okabe & Otsuka, J. Phys. A 58 (2025) 065003 | MC **Correlation-Ratio** (Verhaeltnis der Korrelationsfunktion bei zwei Distanzen) + ML-Phasenklassifikation; honeycomb, kagome, diced; triangular aus "previously calculated data" | honeycomb 0.573, kagome 0.825 (beide "rough estimate") | normierungsfrei (kein Helicity-Faktor) |
| arXiv:2305.00651 = Otsuka, Shiina, Okabe, J. Phys. A 56 (2023) 235001 | Hauptgegenstand laut Abstract: AF-3-State-Potts (square, NNN) und triangular AF-Ising (anisotrope NNN) - ML-Phasenklassifikation, Correlation-Ratio-MC, Level-Spectroscopy. Ein zweiter Such-Snapshot ordnet der Arbeit zusaetzlich XY + 6-state clock auf **triangular** (L = 48/72/96/144/192, Helicity + Correlation-Ratio) zu und nennt: **zitiert Sorokin, T_BKT = 1.418(2) aus dem Helicity-Modul**. Die beiden Snapshots widersprechen sich nicht, sind aber nicht am Primaertext geprueft. | 1.418(2) (Sorokin, Helicity) | **Helicity-Wert; Normierung unbelegt; Sorokin-Primaerquelle per Websuche nicht lokalisierbar** (2026-09-28). Ein eigener Okabe-Otsuka-Correlation-Ratio-Wert fuer triangular ist offen. |
| Butera & Pernici, arXiv:0806.1496 (2008) | Hochtemperatur-Reihe der Korrelationsfunktion, triangular bis Ordnung 20 | J_c = 0.6824(8) -> **T_BKT = 1.4654(17)** (Snippet, zusaetzlich als Zitat in arXiv:1010.3075 bestaetigt) | normierungsfrei (Reihe) |
| arXiv:1010.3075 (PRE 83, 011124) | MC, BKT-artige Perkolation, square + triangular | J_c = 0.6833(6) -> 1.4635(13) (Snippet-Ebene, aus PR #50 uebernommen) | normierungsfrei |
| Butera & Comi, PRB 50, 3052 (1994); cond-mat/9902326 | HT-Reihe triangular bis Ordnung 14 | (Vorgaenger von 0806.1496) | normierungsfrei |
| arXiv:2406.12076 (de Andrade, Jorge, DaSilva; Phys. Scr. 100 065953) | honeycomb, SA + WL, Upsilon_2 / Upsilon_4 / Binder | Snippets 2026-09-28: beta = 1.696(3) / 1.67(1) / 1.724(2) (v2-Stand), T_BKT = 0.576(1); Journal 0.575(8) SA / 0.576(3) WL (aus #50/#52) | #52 Nachtrag 4: Gl. (16)/(19) mit Faktor 4/(3 sqrt 3) = per Flaeche. **Hier nicht unabhaengig verifizierbar.** |
| arXiv:2406.14812 (Jiang) / PTEP 2024 103A02 | honeycomb, NN + Helicity | v1: NN 0.560(9), Helicity 0.571(8); PTEP: NN 0.572(3), Helicity 0.576(4) (Snippets 2026-09-28 bestaetigt) | #52 Nachtrag 3: Faktor 4/(3 sqrt 3) im PTEP-Text = per Flaeche. **Hier nicht unabhaengig verifizierbar.** |
| Caci, Weber, Wessel, PRB 104, 155139 (2021) | Spin-1-Heisenberg-AF honeycomb, QMC (Titel per Snippet bestaetigt) | - | Die in PR #50 zitierte Normierung auf die Einheitszellen-Flaeche ist per Snippet nicht belegbar (nur Titel/Abstract sichtbar). |

### 3.2 Konsequenzen

1. **triangular:** Die Repo-Referenz 1.418 (Core, PHY029/030/037/043) ist die
   **Sorokin-Helicity-Zahl**, keine normierungsfreie Groesse. Sie ist damit
   genau die Klasse von Wert, deren Konvention O1 in Frage stellt. Die
   einzige normierungsfreie triangular-Referenz MIT Fehlerbalken ist die
   HT-Reihe 1.4654(17); sie bildet die Mitte des PHY051-Bandes B_tri.
   Auffaellig: per-Site-Pipeline 1.4007 und Sorokin 1.418 liegen beide
   ~3 % unter 1.465 - dasselbe Muster wie auf honeycomb (per-Site-Pipeline
   trifft die Helicity-Literatur). Ob Sorokin per Site normiert hat, ist
   eine Primaertext-Frage; PHY051 misst, welche Konvention die
   normierungsfreie Reihe trifft.
2. **honeycomb:** Falls #52 Nachtrag 3/4 zutrifft (Faktor 4/(3 sqrt 3) in
   Jiang und de Andrade), sind die Helicity-Literaturwerte 0.576 bereits
   per Flaeche, und die Spannung PHY046 0.565(4) vs 0.572-0.576 ist echt
   (~2 sigma). PHY051 liefert dazu den Quercheck, ob die Pipeline auf einem
   zweiten nicht-quadratischen Gitter mit normierungsfreier Referenz
   ebenfalls tief liegt (Pipeline-Systematik) oder nicht (honeycomb-
   spezifisch).
3. **kagome:** 0.825 ist "rough estimate" ohne Fehlerbalken -> B_kag
   bewusst +/- 2 %.

## 4. Unabhaengige Nachrechnung des O1-Entscheids (Spin-Wellen-Ebene)

Per-Bond-Twist-Schaetzer (Repo, PHY045): Upsilon_site = (1/N) [<sum_b J cos(dtheta_b)(e_b.x)^2> - (1/T) <(sum_b J sin(dtheta_b)(e_b.x))^2>].
Bei T -> 0 mit uniformem Gradienten theta_j = g x_j:
E(g) = E0 + (J g^2 / 2) sum_b (e_b.x)^2 = E0 + (g^2 / 2) N Upsilon_site(0).
Kontinuum: E(g) = E0 + (rho_s / 2) g^2 A. Gleichsetzen:
**Upsilon_site = rho_s A/N = a_s rho_s.** Das NK-Kriterium
rho_s(T_c) = 2 T_c / pi ist fuer rho_s (Flaechendichte) formuliert; auf
Gittern mit a_s != 1 ist Upsilon_site die falsche Groesse. Zahlen (NN = 1):

| Gitter | Bonds/Site | <(e.x)^2> | Upsilon_site(0) | a_s | rho_s(0) = Upsilon_area(0) |
|---|---:|---:|---:|---:|---:|
| square | 2 | 1/2 | 1 | 1 | 1 |
| triangular | 3 | 1/2 | 3/2 | sqrt(3)/2 | sqrt(3) = 1.732 |
| honeycomb | 3/2 | 1/2 | 3/4 | 3 sqrt(3)/4 | 1/sqrt(3) = 0.577 |
| kagome | 2 | 1/2 | 1 | 2/sqrt(3) | sqrt(3)/2 = 0.866 |

Das deckt sich mit PHY045 Teil A (K_fit/Upsilon_area(0) = 1.0004..1.0013) und
mit der Faktor-Angabe 4/(3 sqrt 3) = 1/a_s(honeycomb), die #52 den beiden
honeycomb-Arbeiten zuschreibt. Der staerkste Beleg bleibt PHY045 Teil B:
R = 2 pi eta Upsilon_site / T = a_s auf allen drei Gittern (eta aus <m^2>,
konventionsfrei) - das ist unabhaengig von jeder Literatur-Konvention.

**Bewertung:** Der O1-Entscheid "per Flaeche" ist physikalisch korrekt
und hinreichend belegt. Die Restfrage ist ausschliesslich die Konvention
der zitierten Helicity-Literaturwerte (Sorokin, Jiang, de Andrade) und die
damit verbundene Lage-Spannung, nicht die Theorie.

Zwei Feinheiten, die der Entscheid nicht beruehrt, aber PHY051 beachtet:
(a) Auf dem Rhombus-Torus (tau = exp(i pi/3)) sind Upsilon_x und Upsilon_y
bei endlichem L verschieden; das x/y-Mittel des Repo-Schaetzers ist fuer
3- und 6-zaehlige Gitter im Limes isotrop. (b) Die Torus-Form beeinflusst
die finite-size-Korrekturen (nicht-universelles C), nicht den Grenzwert;
der C-eliminierte Paar-Schaetzer bleibt anwendbar.

## 5. Kleinere Befunde in #50 (nicht blockierend)

- `spec/260926 ... O1 decision v01.md` Abschnitt 3 zitiert Caci/Weber/
  Wessel als Beleg fuer Einheitszellen-Normierung; per Snippet nicht
  belegbar (Abschnitt 3.1). Empfehlung: als `search_corroborated` markieren.
- PHY046-Note Abschnitt 3 nennt die Hypothese "Literatur-Helicity per Site";
  #52 Nachtrag 3/4 widerspricht fuer Jiang/de Andrade. Sollte #52 mergen,
  ist die Note als supersedet zu markieren (nicht umschreiben).
- Der PHY048-Befund "3-Parameter-HKS rauschverstaerkend" ist in PHY051
  als Protokoll-Lehre umgesetzt (Paar-Mittel als Primaer; Diagnostik der
  sigma_FSS-Spannweite ohne v3, nicht entscheidend).
- PHY046 `estimate()`: Variante v3 "groesstes Paar" ist als `Ts[-1]`
  implementiert = groesstes Paar MIT Crossing. Fuer den committeten Lauf
  ohne Folge (alle 5 Paare kreuzen); bei fehlendem Crossing des groessten
  Paares wuerde still das naechstkleinere gezaehlt. PHY051 bindet v2 an das
  groesste Paar allein (None, falls kein Crossing).

- **PHY046 `produce()` (und PHY045 `part_b`/`part_c`) nutzen
  `ProcessPoolExecutor` ohne `mp_context`.** Die Module sind per importlib
  unter Kunstnamen geladen (Dateinamen mit Leerzeichen); spawn-/forkserver-
  Kinder koennen sie nicht per Modulname importieren. Python 3.14 macht
  `forkserver` zum Linux-Default (What's new in Python 3.14, docs.python.org,
  abgerufen 2026-09-28), Windows nutzt `spawn`: Produktionslaeufe mit
  `max_workers > 1` brechen dort im Worker (`ModuleNotFoundError`). Die CI
  (py3.14-Matrix) sieht das nicht, weil alle Tests `max_workers=1` fahren.
  Fix wie in PHY051 (`db90cd2`): expliziter `fork`-Kontext, sonst sequenziell.
  Nicht blockierend fuer #50 (committete Laeufe liefen unter fork); fuer
  jede Neu-Produktion auf py3.14/Windows relevant.

## 6. Offen nach diesem Review

1. Primaertext-Abgleich (Abschnitt 3.0, konkreter Weg genannt) - einziger
   Punkt, der in dieser Umgebung nicht abschliessbar ist.
2. #52: SOURCES-Zeile (Abschnitt 2), dann Codex-Review-Runde abschliessen.
3. PHY049-Produktion (W4-v03, mehrtaegig, 24-h-Budget-Stop) auf einer
   Maschine mit Laufzeit - nicht in einer Agenten-Session.
4. PHY051-Ergebnis (PR #54) nach Lauf in README/CHANGELOG; Interpretation
   nur ueber die vorregistrierte Wahrheitstafel.
5. Cross-Family-Review vor jeder externen Aussage bleibt bindend.
