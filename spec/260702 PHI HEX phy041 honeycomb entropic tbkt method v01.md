# Methoden-Spec: PHY041 — WL-entropischer Helicity-Modul auf honeycomb

> **Status:** Vertrags-Quelle (spec/). Repo-native, 2026-07-02. Setzt den in
> der WL-Methoden-Spec (`260616 … wang-landau … method v01.md`, §6) benannten
> **honeycomb-Follow-up** um. Code folgt dieser Spec (`src/260702 PHY041 …`).

## 1. Motivation

PHY040 hat den entropischen Sampler (Wang-Landau + 1/t) am square-Goldstandard
validiert und Statistik- vom finite-size-Limit getrennt. Das eigentliche
offene Ziel der Serie PHY031–039 ist **honeycomb** (Referenzen 0.573
arXiv:2501.07388 / 0.576(3) arXiv:2406.12076 — letztere bestimmt T_BKT mit
exakt dieser Methode: WL + 2.-/4.-Ordnungs-Helicity-Modul). PHY041 traegt den
Sampler dorthin und haertet ihn dabei dreifach (web-recherchiert).

## 2. Methode (Delta zu PHY040)

Analyse-Mathematik (kanonische Rueckgewichtung, Υ₂/Υ₄, C-eliminierte
WM-Paare) **1:1 aus PHY040 importiert** (single source of truth). Neu:

1. **Gitter-agnostischer Kernel.** Der WL-/Produktions-Walk arbeitet auf
   (Adjazenz, Kanten, Bond-Projektionen) statt square-Hardcode; skalar-
   optimiert (Block-RNG je Sweep + `math.cos` statt numpy auf 3-Element-
   Slices): **~18× schneller** als der PHY040-Innerloop, physik-identische
   Akzeptanz `log U ≤ ln g(b) − ln g(b')`. Aequivalenz-Gates:
   `tests/test_phy041_honeycomb_wl.py` (lokale ΔE == volle Energie-Differenz;
   Aggregat-Orakel Υ₂(0) = 3/4 J exakt).
2. **Echte 1/t-Phase (Belardinelli & Pereyra).** Die lnf-Halbierungen nutzen
   das B&P-Kriterium „alle Bins besucht" (H_min ≥ 1) statt strikter
   Flachheit; der Wechsel auf lnf = 1/t wird **am Halbierungs-Ereignis**
   geprueft. Befund (Selbst-Debug, honeycomb L=12): im PHY040-Kernel waren
   die strikten Flachheits-Epochen so langsam, dass lnf(t) nie 1/t erreichte
   — die 1/t-Politur griff **nie**, und der Standard-WL-Saettigungsfehler
   blieb (Walker-abhaengiger g(E)-Bias, auf L=12 ~0.02 in ⟨E⟩/N bzw. bis
   0.06 in Υ₂(0.65) messbar). Mit echter 1/t-Politur (lnf_final = 1e-5,
   d.h. ≥ 1e5 Sweeps Politur): ⟨E⟩/N-Abweichung < 0.006. Genau die
   Saettigung, die B&P (JCP 127, 184105 (2007)) beheben.
3. **Auto-Energiefenster + hartes Leak-Gate.** Das BKT-Fenster wird je L aus
   kurzen Wolff-Ankern bei T_lo=0.50 / T_hi=0.70 bestimmt (Mittel ± k·std,
   Minimal-Margen, Grundzustands-Floor-Guard bei −1.5 J + 0.03 gegen
   unreachable Bins / dynamical traps, arXiv:1508.01888). PHY040 hatte das
   Fenster hartkodiert und das Rand-Leck nur im Kommentar dokumentiert;
   PHY041 prueft HART: kanonisches Gewicht in den aeussersten 3 Bins je Rand
   < 1e-3 fuer **jedes** Analyse-T (Gate `PASS_NO_CANONICAL_EDGE_LEAK`).

**Lauf-Konfiguration:** L ∈ {12, 16, 24} (N = 288/512/1152), T-Gitter
0.52…0.67 (Δ=0.005), lnf_final = 1e-5, prod_sweeps = 30000, master_seed = 42.
RNG-Streams: Anker 600+s+1000L, WL-Walker 650+L, Wolff-Referenz 660+s+1000L
(kollisionsfrei zu PHY031/032/033/040).

## 3. Korrektheit / Validierung (arXiv-unabhaengig, zweifach)

- **VAL-A (frische Wolff-Referenz, L=12):** WL-⟨E⟩(T)/N und WL-Υ₂(T) muessen
  direktes Wolff treffen (Toleranzen wie PHY040: 0.03 / 0.04); Υ₂ dabei wie
  in der WL-Pipeline ueber beide Twist-Richtungen gemittelt.
- **VAL-B (Drift-Guard gegen committed Evidenz):** die WL-Kurven muessen das
  im Repo festgeschriebene PHY032-Wolff-Messgitter Υ(T,L) (L=12/24, 8 T-
  Punkte je L, `results/260607 PHY032 …`) innerhalb max(4·SEM, 0.025)
  reproduzieren — deterministisch, zur Laufzeit geparst (Muster PHY037).
- Physik-Hypothesen (Paar-Schaetzer nahe Referenz) sind bewusst **kein**
  Build-Breaker (Negativ-Results sind Buerger erster Klasse).

## 4. Ergebnis (2026-07-02)

> Zahlen aus `results/260702 PHY041 … report.txt` (deterministisch, seed=42).

Siehe Gate-Log; Kurzfassung im README (Mess-Stand). Kern-Punkte:

- Pipeline-Integritaet: alle Gates PASS (Orakel exakt, VAL-A/B, kein Leak,
  Υ₂ glatt); Υ₂/Υ₄(T) fuer alle drei L **rauschfrei** aus je EINEM Lauf.
- T_BKT (C-eliminierte WM-Paare auf glatten WL-Kurven): siehe Report —
  ehrlicher Vergleich gegen **beide** Referenzen (0.573 / 0.576(3)) und
  gegen den Wolff-Paar-Stand PHY032 (0.6014 bei (12,24)).
- Der 4.-Ordnungs-Dip (PHY039-Observable) ist unter dem entropischen Sampler
  deterministisch glatt lokalisierbar (Lage je L im Report).

## 5. Wiederverwendbarer API

| Funktion | Zweck |
|---|---|
| `wl_entropic_lattice` | gitter-agnostischer WL-Kern (B&P-1/t) + Produktion |
| `wang_landau_honeycomb` | Auto-Fenster (Wolff-Anker) + WL fuer ein L |
| `window_from_anchors` | Fenster-Logik (pur, getestet) |
| `canonical_edge_leak` | Leak-Mass fuer das Rand-Gate |
| `parse_phy032_grid` | Drift-Guard-Parser der committed PHY032-Evidenz |
| `wolff_reference` | frische Wolff-Quervalidierung (E, Υ₂ beide Richtungen) |

Triangular/kagome sind jetzt reine Adjazenz+Projektions-Frage (PHY033-
Geometrie vorhanden) — der Kern ist wiederverwendbar.

## 6. Naechste Stufe (Vertrag)

- **Groessere L (32/48) + FSS wie arXiv:2406.12076** (2.-+4.-Ordnung
  kombiniert): mit dem ~18×-Kernel und B&P-1/t lokal in Reichweite
  (L=24 lief in ~2 min); fuer L≥48 Politur-Budget (lnf_final) und
  Produktions-Statistik neu bemessen.
- **Walker-Mittelung** (mehrere unabhaengige g(E)-Walker) als optionale
  weitere Varianz-Reduktion, falls die Paar-Schaetzer bei groesseren L
  Seed-Systematik zeigen.

## 7. Provenance / Referenzen

- Wang & Landau, PRL 86, 2050 (2001); Belardinelli & Pereyra, JCP 127,
  184105 (2007) (cond-mat/0702414); kontinuierliche Modelle: PRE 89, 013311
  (2014); dynamical traps: arXiv:1508.01888; 2D-XY-DOS: cond-mat/0611039;
  Methodik-Vorbild honeycomb: arXiv:2406.12076; WM-Paare: PRB 37, 5986(R)
  (1988), arXiv:1302.2900; Referenz-Uebersicht: arXiv:2501.07388.
- Code/Evidenz: `src/260702 PHY041 …`; Gate-Log `results/260702 PHY041 …`;
  Tests `tests/test_phy041_honeycomb_wl.py`.

## 8. Definition of Done

- [x] Gitter-agnostischer WL-Kern (Adjazenz + Bond-Projektionen), skalar-
      optimiert, Aequivalenz-/Orakel-Gates.
- [x] Echte B&P-1/t-Phase (Saettigungs-Bug des PHY040-Phasenschemas behoben,
      Befund dokumentiert).
- [x] Auto-Energiefenster aus Wolff-Ankern + hartes Leak-Gate.
- [x] Zweifache arXiv-unabhaengige Quervalidierung (frisches Wolff +
      committed PHY032-Gitter als Drift-Guard).
- [x] T_BKT(honeycomb) via C-eliminierte Paare auf glatten WL-Kurven,
      ehrlich gegen beide Referenzen; Gate-Log in results/.
