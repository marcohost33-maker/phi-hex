# Phi-Hex

> **Status:** Forschungs-Repo (privat) | Aktive Forschungslinie (letzte Haertung 2026-06-11). Kern-Engine kompiliert; Selftest- + T_BKT-Mess-Evidenz in results/ (jetzt inkl. honeycomb WM-Log-Fit + Seed-Bootstrap-CI PHY032 + kagome PHY033). Testsuite 78/78 PASS (lokal verifiziert 2026-06-12; davon 72 schnelle Korrektheits-Gates blockierend in CI, 6 slow Mess-Laeufe lokal; inkl. Kern-Vektorisierungs-Aequivalenz-Gates).
> **Lineage/Provenance:** siehe `SOURCES.md` (SHA-256 je Quelldatei) | **Lizenz:** Apache-2.0

XY-Modell / BKT-Physik auf Dreiecks- und Honeycomb-Gittern: Helicity-Modulus, Nelson-Kosterlitz-Sprung, Wolff-Cluster, Finite-Size-Scaling.

## Konvention: Helicity-Modulus per-Site (Stand 2026-06-04)

> **Konventionswechsel 2026-06-04 (Marco-Entscheid, Wissenschafts-Audit S2/S3).**
> Der Helicity-Modulus wird jetzt **per-Site** normiert (Division durch die
> Site-Zahl N), Literatur-Standard (arXiv:2406.12076). Spinwellen-Grenzwert
> Dreiecksgitter: **Υ(0) = 1.5 J** (exakt, getestet, tol 0.03).
>
> **Vorher (superseded):** Normierung durch die FLAECHE A = N·√3/2 -> Υ(0) = J·√3 ≈ 1.732.
> Das ueberschaetzte T_BKT um **~15%** und war durch die zu laxe Toleranz tol=0.15
> nie hart geprueft. BKT-Sprung-Kriterium: Nelson-Kosterlitz Υ(T_BKT) = 2·T_BKT/π.

### T_BKT-Mess-Stand (per-Site, Neumessung 2026-06-04)

Wolff-Sampling auf dem Dreiecks-Torus, Nelson-Kosterlitz-Crossing
(`src/260604 PHY030 ...`, Evidenz: `results/260604 PHY030 ... report.txt`):

| Gitter | T_BKT (per-Site, neu) | Referenz | Methode |
|---|---|---|---|
| triangular | **1.4007 ± 0.0081** (Weber-Minnhagen-Log-Fit, v02) | 1.418 (arXiv:2501.07388) | Wolff + WM-Log-Korrektur, chi^2-Fit |
| triangular | ~1.42 (v01-Schranken: L=19-Crossing 1.456; 1/lnL-Extrap. 1.382) | 1.418 (arXiv:2501.07388) | Wolff + NK-Crossing, FSS (konservativ) |
| square | 0.893 (PHY028, V&V) | 0.893 (arXiv:2501.07388) | Sandvik-Paar (C-frei) |
| honeycomb | **PHY032 (s.u.)** Sandvik-Paar(24,48) + WM-Log-Fit, beide mit Seed-Bootstrap-CI | 0.573 (arXiv:2501.07388) UND 0.576(3) (arXiv:2406.12076) | Wolff + WM-Fit + Sandvik-Paar |
| honeycomb | **0.595** (PHY031 v01, Sandvik-Paar L=12/24; +3.9% vs 0.573) — superseded durch PHY032 (groesseres L) | 0.573 / 0.576(3) | Wolff + Sandvik-Paar (C-frei) |
| kagome | **0.8479 ± 0.005 (WM-Log-Fit, CI[0.8414,0.8507]; +2.78%)** · Sandvik-Paar(24,36) 0.8377 (+1.54%) — PHY033, beide mit Seed-Bootstrap-CI | 0.825 (arXiv:2501.07388, "rough estimate") | Wolff + WM-Fit + Sandvik-Paar (per-Site, z=4, Υ(0)=1 exakt) |

> **Honeycomb-Referenz: bewusst BEIDE Literaturwerte (kein Kanon gekuert).**
> Die Abweichung wird gegen **0.573** (arXiv:2501.07388, multi-lattice MC) UND
> gegen **0.576(3)** (arXiv:2406.12076, dedizierte honeycomb-Studie, 4th-order
> helicity + WM) ausgewiesen. Die ~0.5%-Divergenz der Quellen ist eine reale
> Konvention/Methodik-Differenz und wird ehrlich gezeigt, statt einen Wert als
> "richtig" zu kueren (Marco/Vero-Entscheid 2026-06-07).

**v02 (SOTA-Praezisierung, 2026-06-04):** Der Weber-Minnhagen-Log-Korrektur-Fit
`Υ(T_BKT,L) = (2 T_BKT/π)·(1 + 1/(2 ln L + C))` (Weber & Minnhagen PRB 37,
5986(R) (1988); Methodik arXiv:1302.2900) liefert einen Punktschaetzer
**1.4007 ± 0.0081 (−1.22%)** — naeher an der Referenz 1.418 als beide
konservativen v01-Schranken. Querscheck per C-freier Paar-Methode: 1.3801.
Beide v02-Schaetzer liegen knapp UNTER 1.418; bei nur drei kleinen Gittern
(L=9..19) sind sub-leading Log-Korrekturen eine reale Praezisions-Grenze
(ehrlich ausgewiesen). **v01 bleibt gueltig** als konservative Schranken-Methode.
Code: `src/260604 PHY030 ... v02 wm-logfit.py`, Evidenz:
`results/260604 PHY030 v02 wm-logfit report.txt`.

**PHY032 (honeycomb WM-Log-Fit + groessere L + Seed-Bootstrap-CI, 2026-06-07):**
Wendet die Weber-Minnhagen-Methodik (1:1 die PHY030-v02-Fit-Funktionen,
gitter-agnostisch) auf honeycomb an und erweitert das L-Set von PHY031 v01
(6/12/24) auf **L = 12/24/48** (N bis 4608). Zusaetzlich tragen jetzt **alle**
T_BKT-Schaetzer ein einheitlich propagiertes **Seed-Bootstrap-CI** (Perzentil,
n_boot=2000, seed=42) + Jackknife-Quercheck. Evidenz:
`results/260607 PHY032 honeycomb wm-logfit bootstrap report.txt`.

| Methode (PHY032) | T_BKT(honeycomb) | vs 0.573 | vs 0.576(3) |
|---|---|---|---|
| Sandvik-Paar (24,48), C-frei — **groesstes L** | **0.5917**  CI[0.587, 0.598] | +3.27% | +2.73% |
| Sandvik-Paar (12,48), C-frei | 0.5968  CI[0.590, 0.598] | +4.15% | +3.60% |
| Sandvik-Paar (12,24), C-frei | 0.6014  CI[0.599, 0.610] | +4.95% | +4.41% |
| Weber-Minnhagen Fixed-T-Fit (alle L, argmin@0.5975) | 0.5964  CI[0.594, 0.599] | +4.08% | +3.54% |
| PHY031 v01 Paar (12,24) — superseded | 0.595 | +3.90% | — |

**Befund (ehrlich, NICHT getunt):** Groesseres L senkt den Sandvik-Paar-Bias
monoton (Paar(24,48) +3.27% vs Paar(12,24) +4.95%) — die Audit-Hypothese
bestaetigt sich, analog triangular. Der Rest-Bias bei lokal-vertretbaren L
(max 48; Referenz-Studien nutzen L=48..192 bzw. 8..128) ist erwartbar
finite-size-dominiert. Die CIs erfassen die **Seed-Streuung**, nicht den
finite-size-Bias; sie ueberdecken 0.573/0.576 NICHT — der Rest-Abstand ist
also reales finite-size, kein statistisches Rauschen. Der WM-Fit liegt hier
(anders als bei triangular) **nicht** unter dem Paar-Schaetzer; bei nur drei
kleinen L sind sub-leading Log-Korrekturen eine reale Grenze (ehrlich
ausgewiesen). Code: `src/260607 PHY032 honeycomb wm-logfit bootstrap v01.py`.

**PHY034 (HKS-L→∞-Extrapolation, 2026-06-11):** Die in PHY032 als "erwartbar
finite-size" nur behauptete Grenze wird jetzt **getestet**. Web-recherchiert:
Hsieh-Kao-Sandvik (J. Stat. Mech. (2013) P09001, arXiv:1302.2900) zeigen, dass
die C-eliminierte Paar-Schaetzung T_BKT(L,2L) selbst noch eine sub-leading
log-Drift traegt und zum L→∞-Limes extrapoliert werden muss. Rein analytisch
auf den committed PHY032-Daten (kein neuer MC-Lauf): die Verdopplungspaare
0.6014/0.5917 extrapolieren im HKS-Standard (u=1/(ln L)², geom. Mittel) zu
**T_BKT = 0.5741** (+0.19% vs 0.573 / −0.33% vs 0.576) — das LEGT NAHE, dass
der +3 %-Offset ein sub-leading-log-Artefakt ist. Ehrliche Doppel-Unsicherheit:
Methoden-Systematik [0.547, 0.583] + statistische 95%-CI [0.553, 0.595] (die
2-Punkt-Extrapolation verstaerkt den Paar-Fehler ~3×). **Wichtig (Selbst-
Falsifikation, s. PHY035):** dieser Wert ist NUR ein 2-Punkt-Konsistenz-Check
und haengt allein am Paar (24,48); die dichtere Leiter widerlegt seine
Robustheit. Code: `src/260611 PHY034 honeycomb hks extrapolation v01.py`.

**PHY035 (dichte L-Leiter — falsifiziert die PHY034-Lesart, 2026-06-11):**
Frische Wolff-Messung auf L = 8/12/16/24/32 (drei verschachtelte Paare
(8,16)/(12,24)/(16,32); dank Kern-Vektorisierung in ~3 min lokal). Die Paare
liegen **flach bei ~0.597–0.601** — KEINE saubere Drift gegen 0.573; die echte
≥3-Punkt-Extrapolation ergibt **~0.605** (CI[0.577, 0.616]), *oberhalb* der
Referenzen. Damit ist PHY034 (0.574) **nicht robust**: das Ergebnis hing allein
am größten Paar (24,48), das hier fehlt. **Ehrlicher Stand:** der +3 %-Offset
ist *plausibel* finite-size (größtes Paar + Literatur stützen es), aber mit den
lokal vertretbaren L (≤48) **nicht beweisbar** — eine belastbare Konvergenz
braucht L ≥ 48..192 (HKS), jenseits des lokalen Budgets. Negativ-Result, ehrlich
ausgewiesen. Code: `src/260611 PHY035 honeycomb hks dense ladder v01.py`,
Evidenz: `results/260611 PHY035 honeycomb hks dense ladder report.txt`.

**PHY037 (HKS gitter-übergreifend — ehrliche Methoden-Bilanz, 2026-06-11):**
Derselbe HKS-Apparat, uniform auf die committed Paar-Schätzer aller drei Gitter
angewandt (zur Laufzeit aus den PHY030/032/033-Reports geparst):

| Gitter | Paare (L) | L→∞-Extrap (p=2,geom) | bestes Paar | Verdikt |
|---|---|---|---|---|
| kagome | 12/24/36 | **0.8236** (−0.17% vs 0.825) | 0.8377 (+1.54%) | **sauber** |
| honeycomb | 12/24/48 | **0.5747** (−0.23% vs 0.576) | 0.5917 (+2.73%) | **sauber** (nur mit L=48) |
| triangular | 9/13/19 | **1.4591** (+2.90% vs 1.418) | 1.3888 (−2.06%) | **überschießt** |

**Befund:** Die Klein-L-Paar-Extrapolation konvergiert sauber für kagome und
honeycomb (letzteres **nur** wenn L=48 enthalten ist — PHY035 ohne L=48 kippt
auf ~0.605), **überschießt** aber für triangular bei sehr kleinem L (≤19), wo
das blanke größte Paar näher an der Referenz liegt. Die Methode ist also
gitter-abhängig und empfindlich auf das größte enthaltene L — **kein
universelles Wundermittel** (deckt sich mit HKS: L=48..192 nötig). Code:
`src/260611 PHY037 hks multilattice extrapolation v01.py`.

**PHY036 (größere Leiter L bis 64 — L hilft NICHT, 2026-06-11):** Frische
Wolff-Messung L = 16/24/32/48/64 (N bis 8192; ~11 min lokal dank
Vektorisierung). Ergebnis: das **größere L stabilisiert die honeycomb-
Extrapolation NICHT — es verschlechtert sie**. Das (32,64)-Paar springt auf
0.6178 mit breitem CI[0.598,0.635] (L=64 ist bei 8 Seeds **rausch-dominiert**,
sichtbar an Υ(L=64) nahe T_BKT). Die ≥3-Punkt-Extrapolation ergibt **0.639**
(+11.5 %). Damit **wandert** der Extrapolations-Limes mit dem L-Set:
**0.574** (PHY032, L≤48) → **0.605** (PHY035, L≤32) → **0.639** (PHY036, L≤64).
Cross-Check intakt: (24,48) reproduziert PHY032s 0.5917 bit-identisch.
**Korrekte Schlussfolgerung:** Honeycomb-T_BKT ist aus der lokal vertretbaren
MC-Statistik **nicht robust extrapolierbar**; die ehrlichste Aussage bleibt der
Roh-Paar-Stand ~0.59 (L≈48, +3 % über Literatur, innerhalb finite-size).
Belastbare Konvergenz braucht L ≥ 48..192 **mit** hoher Seed-Statistik —
jenseits des lokalen Budgets. Code:
`src/260611 PHY036 honeycomb hks large ladder v01.py`.

**PHY038 (Hochstatistik — es ist KEIN Statistik-Problem, 2026-06-12):**
Dieselbe Leiter mit **doppelter Statistik** (16 statt 8 Seeds, 240 statt 140
Sweeps; 32 min). Ergebnis: das (32,64)-Paar bleibt **0.6181** (8 Seeds: 0.6178)
— **unverändert**, CI nicht enger. Der hohe Wert ist also **kein Rauschen**;
die Paare bleiben nicht-monoton (0.597/0.594/0.618), die Extrapolation 0.647
(+13 %). **Damit ist die honeycomb-Nicht-Auflösbarkeit definitiv ein Methoden-/
finite-size-Limit, kein Statistik-Problem** — mehr Seeds/Sweeps lösen es nicht.
Eine echte Auflösung bräuchte größere L (≥96..192) **und** ein feineres
T-Gitter nahe dem Crossing (HPC). Code:
`src/260612 PHY038 honeycomb hks highstat v01.py`.

**Superseded:** alle vor 2026-06-04 publizierten T_BKT-Zahlen, die auf der
Flaechen-Normierung beruhen (u.a. PHY026-Report "Υ(T→0)=J·√3"), sind durch
diesen Mess-Stand ersetzt. Details in `CHANGELOG.md`.

## Kern

- **Engine:** `src/260602 PHI HEX core v2 2 hardened.py` (Kern) + PHY024-038-Experiment-Serie
- **Evidenz:** Selftest-Report + PHY025-038 Reports in `results/`; Testsuite in `tests/`; methodischer Audit (4 Responses) in `spec/`. Mess-Stand T_BKT(triangular) = 1.42 ± 0.04 (per-Site, 2026-06-04; Referenz 1.418); T_BKT(honeycomb) = 0.5917 CI[0.587, 0.598] (PHY032 groesstes L, Sandvik-Paar (24,48), 2026-06-07; Referenzen 0.573 / 0.576(3)) — supersedet PHY031 v01 (0.595); T_BKT(kagome) = 0.8479 CI[0.8414, 0.8507] (PHY033 WM-Log-Fit, 2026-06-09; Referenz 0.825 "rough estimate").

## Struktur

```
src/        Engines + Experiment-Code
spec/       Specs, Theorie, Audits (gehaertete Versionen)
results/    Selftest-/Gate-Logs, Reports, Negativ-Results
archive/    Vorgaenger-Versionen (Lineage)
docs/ADR/   Architektur-Entscheide
SOURCES.md  Provenance: Quelle + SHA-256 + mtime je Datei
```

## Reproduzieren

Die Selftests laufen direkt in der Engine (`python <engine>.py`; Details im
Engine-Header). Die Testsuite (numpy/scipy/pytest noetig):

```
pytest                 # volle Suite (inkl. slow Mess-Laeufe), 78/78 PASS
pytest -m "not slow"   # nur schnelle Korrektheits-Gates (72, inkl. WM-/Sandvik-/Bootstrap-/Vektorisierungs-/HKS-/Multilattice-Orakel) — CI-Gate
python "src/260604 PHY030 triangular tbkt per-site v01.py"            # T_BKT(tri)-Schranken (konservativ)
python "src/260604 PHY030 triangular tbkt per-site v02 wm-logfit.py"  # T_BKT(tri) Weber-Minnhagen-Log-Fit
python "src/260607 PHY031 honeycomb tbkt per-site v01.py"             # T_BKT(honeycomb) Wolff + Sandvik-Paar (kleine L)
python "src/260607 PHY032 honeycomb wm-logfit bootstrap v01.py"       # T_BKT(honeycomb) WM-Log-Fit + groessere L + Seed-Bootstrap-CI (~24 min lokal)
python "src/260609 PHY033 kagome tbkt per-site v01.py"               # T_BKT(kagome) Wolff + WM-Log-Fit + Sandvik-Paar + Seed-Bootstrap-CI (~9 min lokal)
```

CI (`.github/workflows/ci.yml`) prueft drei blockierende Gates: Lint
(`ruff check`, Baseline in `ruff.toml`), Syntax (`compileall`) und die
schnellen Korrektheits-Tests (`pytest -m "not slow"`) ueber die Python-Matrix
3.10/3.11/3.12 (Deps aus `requirements-dev.txt`, pip-Cache). Die vollen
physikalischen Mess-Laeufe (slow, Wolff-Sampling) laufen lokal.

## Verwandte Repos

`coworker-dde` (kanonische v3-Engine-Quelle, 3DDE-Kanon) | `HQST` (p1_tasks-Floquet-Linie) |
`phi-hex` `hex-hqst` `u2-dualspec` `u6-bifurcation` `hexa-ntk` (diese Serie, 2026-06-04)

---
*Coworker Research / Coworkerz | Repo-Anlage 2026-06-04 nach arbeitsschablone_forschungs-repo-anlage*
