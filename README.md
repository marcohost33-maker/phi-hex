# Phi-Hex

> **Status:** Forschungs-Repo (privat) | Aktive Forschungslinie (letzte Haertung 2026-06-11). Kern-Engine kompiliert; Selftest- + T_BKT-Mess-Evidenz in results/ (jetzt inkl. honeycomb WM-Log-Fit + Seed-Bootstrap-CI PHY032 + kagome PHY033). Testsuite 53/53 PASS (lokal verifiziert 2026-06-11; davon 48 schnelle Korrektheits-Gates blockierend in CI, 5 slow Mess-Laeufe lokal; inkl. Kern-Vektorisierungs-Aequivalenz-Gates).
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

**Superseded:** alle vor 2026-06-04 publizierten T_BKT-Zahlen, die auf der
Flaechen-Normierung beruhen (u.a. PHY026-Report "Υ(T→0)=J·√3"), sind durch
diesen Mess-Stand ersetzt. Details in `CHANGELOG.md`.

## Kern

- **Engine:** `src/260602 PHI HEX core v2 2 hardened.py` (Kern) + PHY024-033-Experiment-Serie
- **Evidenz:** Selftest-Report + PHY025-033 Reports in `results/`; Testsuite in `tests/`; methodischer Audit (4 Responses) in `spec/`. Mess-Stand T_BKT(triangular) = 1.42 ± 0.04 (per-Site, 2026-06-04; Referenz 1.418); T_BKT(honeycomb) = 0.5917 CI[0.587, 0.598] (PHY032 groesstes L, Sandvik-Paar (24,48), 2026-06-07; Referenzen 0.573 / 0.576(3)) — supersedet PHY031 v01 (0.595); T_BKT(kagome) = 0.8479 CI[0.8414, 0.8507] (PHY033 WM-Log-Fit, 2026-06-09; Referenz 0.825 "rough estimate").

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
pytest                 # volle Suite (inkl. slow Mess-Laeufe), 53/53 PASS
pytest -m "not slow"   # nur schnelle Korrektheits-Gates (48, inkl. WM-/Sandvik-/Bootstrap-/Vektorisierungs-Orakel) — CI-Gate
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
