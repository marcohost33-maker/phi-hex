# Methoden-Spec: PHY042 - honeycomb WL-FSS L=24/32/48 + getrennte Y2/Y4-Kanaele + Multi-Walker-Systematik

> **Status:** AKTIV. Setzt die in PHY041 (Spec + README "Naechste Stufe nach
> PR #18") vertraglich benannten Punkte 2-4 um: L=32/48 mit demselben Kernel,
> Upsilon_2-/Upsilon_4-FSS getrennt, mehrere g(E)-Walker als Systematik-Check.
> Referenzvergleiche laufen ausschliesslich ueber das Band aus
> `spec/260703 PHI HEX honeycomb reference conventions audit v01.md`.

## 1. Zweck

PHY041 hat den gitter-agnostischen Wang-Landau/1-t-Apparat auf honeycomb bei
L=12/16/24 validiert (OVERALL=PASS) und als FINDING festgehalten: die
C-eliminierten WM-Paare steigen mit L (0.5951/0.6029/0.6087) - der
finite-size-Bias der kleinen Gitter dominiert. PHY042 beantwortet die daraus
vertraglich folgende Frage: **Was machen die Paar-Schaetzer bei L=32/48, und
liefern Y2- und Y4-Kanal getrennt konsistente Aussagen?**

PASS = Pipeline-Integritaet (Orakel, Quervalidierung, Leak, Glaette,
Walker-Reproduzierbarkeit). Der Physik-Befund ist FINDING, kein Build-Breaker
(AGENTS.md: Negativ-Results sind Buerger erster Klasse).

## 2. Methode (Delta zu PHY041)

Kernel, Analyse-Primitiven und Fenster-Logik werden 1:1 aus PHY041/PHY040
importiert (single source of truth, kein Reinvent):

- WL-g(E) mit echter Belardinelli-Pereyra-1/t-Phase (JCP 127, 184105 (2007)),
  lnf_final=1e-5; Produktion mit mikrokanonischen Aggregaten; kanonische
  Rueckgewichtung -> glatte Y2(T)/Y4(T); C-eliminierte WM-Paare.
- Auto-Energiefenster je L aus Wolff-Ankern (T_lo=0.50, T_hi=0.70),
  Grundzustands-Floor-Guard. **Neu:** die Fenster werden im Parent EINMAL je
  L bestimmt und an alle Walker desselben L uebergeben (identisches Binning,
  Walker-Kurven direkt vergleichbar).
- **Skalierte Produktionsphase:** prod_sweeps = max(30000, 60 * nbins),
  damit die Bin-Belegungsdichte (~60 Samples/Bin, PHY041-Niveau bei L=24)
  bei L=32/48 nicht ausduennt. Ausnahme L=24: exakt 30000 (siehe §3 Bruecke).
- **Multi-Walker:** je 3 unabhaengige g(E)-Walker bei L=32 UND L=48
  (Walker 0: Standard-Stream 650+L; Walker w>=1: Stream 90000+1000*w+L,
  kollisionsfrei zu allen bisherigen Stream-Vertraegen). Die Streuung der
  Y2-Kurven ueber Walker ist der Systematik-Schaetzer des Samplers und wird
  auf die Paar-T_BKT propagiert (Paar-Neuberechnung ueber alle
  Walker-Kombinationen).
- Unabhaengige Jobs (L, Walker) duerfen parallel laufen (Prozess-Pool);
  jede Job-RNG ist vollstaendig durch (master_seed, stream) bestimmt -
  Determinismus unabhaengig vom Scheduling.

## 3. Deterministische Bruecke zu PHY041

Der L=24-Lauf verwendet exakt die PHY041-Parameter (master_seed=42,
Streams 600/650, lnf_final=1e-5, prod_sweeps=30000, identisches T-Gitter
0.52..0.67 Schritt 0.005). Damit ist die L=24-Kurve die bitgleiche
Reproduktion des committed PHY041-Stands - verifizierbar ueber:

- [VAL-B] PHY032-Drift-Guard (committed Messgitter L=24, Toleranz
  max(4*SEM, 0.025) wie PHY041);
- [VAL-C] Y4-Dip L=24 == 0.6500 (committed PHY041-Reportwert, exakt auf
  dem T-Gitter).

## 3b. Haertung nach Lauf 1 (2026-07-06) — Lineage

Lauf 1 (Einzel-Walker bei L=48, Walker nur bei L=32) hat ein reales
Negativ-Result geliefert (Gate-Log in `results/`): die Y2-Kurven sind bei
L>=32 oberhalb T~0.60 sampler-limitiert — Walker-Spread bis 0.14 im
Kernfenster, FSS-Ordnungsverletzung (Y2(L=32) < Y2(L=48) bei T>=0.62),
Y4-Dip bei L=32 nicht walker-robust (0.615/0.64/0.66). Ursache:
Dekorrelation der vortex-reichen Konfigurationen in der Produktionsphase,
kein Bug (Leak-Gate und <E>-Validierung blieben unauffaellig).

Konsequenz-Haertungen (v01b, in-place vor Merge, PR-Draft):

1. Multi-Walker an BEIDEN grossen L (32 und 48), Hauptanalyse auf
   Walker-Mittel-Kurven.
2. **Validitaets-Domaene** je L aus dem Walker-Spread: zusammenhaengend ab
   dem unteren T-Rand, solange Spread < 0.04 (= VAL-A-Y2-Toleranz);
   0.02/0.01 als strengere Stufen ausgewiesen.
3. **Paar-Quotierbarkeit:** ein Paar ist nur quotierbar, wenn sein
   Crossing innerhalb beider Domaenen liegt; sonst explizites NR.
4. **Coverage-Massen-Gate:** kanonische Masse (volle lng) auf
   produktions-unbesetzten Bins < 1e-3 fuer jedes in-Domaene-T
   (verallgemeinert das Rand-Leak-Gate; Lauf 1: L=48 nur 1043/1641 Bins
   in der Produktion besetzt).
5. VAL-A-Y2-Gate gilt nur fuer in-Domaene-Punkte; ausserhalb liegende
   Punkte werden als Evidenz ausgewiesen, nicht als Gate.

Negativ-Result-IDs: **NR-PHY042-01** (Y2-Hoch-T-Sampler-Limit L>=32),
**NR-PHY042-02** (Paare mit Crossing ausserhalb der Domaene nicht
belastbar), **NR-PHY042-03** (Y4-Dip L>=32 nicht walker-robust bei diesem
Budget).

## 4. Gates (PASS = Pipeline-Integritaet)

| Gate | Kriterium |
|---|---|
| A aligned-Orakel | Y2(0) = 3/4 J exakt (Aggregat-Maschinerie) |
| LEAK | kanonisches Randgewicht < 1e-3 fuer jedes (L, Walker, T) |
| COVER | kanonische Masse auf unbesetzten Bins < 1e-3 (in-Domaene-T) |
| VAL-A Wolff | frisches Wolff L=32 (T=0.55/0.60/0.65): dE<0.03 immer; dY2<0.04 fuer in-Domaene-T |
| VAL-B PHY032 | L=24-Kurve trifft committed PHY032-Gitter (Toleranz s.o.) |
| VAL-C PHY041 | Y4-Dip L=24 bei T=0.6500 (committed Reportwert) |
| Y2-Glaette | Y2(T) L=48 (Mittel-Kurve) im Kernfenster [0.54,0.66] streng monoton fallend |
| PAIR-2432 | Crossing des Paars (24,32) liegt in beiden Validitaets-Domaenen |

## 5. Auswertungs-Vertrag (getrennte Kanaele)

- **Y2-Kanal:** C-eliminierte WM-Paare (24,32), (24,48), (32,48) auf den
  glatten Kurven; Einordnung NUR relativ zum Referenzband; explizite
  Trendaussage gegen die PHY041-Paare (steigt der Paar-Schaetzer weiter
  mit L, ist das ein ehrliches Negativ-Result fuer die
  kleine-L-Paar-Extrapolation).
- **Y4-Kanal:** Dip-Lage von N*Upsilon_4 je L (24/32/48), Rand-Flag,
  Trendaussage (erwartet: Dip wandert mit L nach unten Richtung T_BKT,
  vgl. Y4-beta-Kanal des Referenzbands).
- **Systematik:** Walker-Spread in Y2 (L=32) und propagiert als
  Delta-T_BKT auf die Paare (24,32)/(32,48). Statistik- vs.
  finite-size-Trennung im Stil PHY040.

## 6. Ehrlichkeit / Grenzen

- L=48 ist weiterhin endlich; auch ein monotoner Paar-Trend ersetzt keine
  unendliche-L-Extrapolation. Kein "neuer T_BKT-Bestwert" ohne
  Cross-Family-Review.
- Die Walker-Systematik misst NUR den Sampler (g(E)-Bias), nicht den
  finite-size-Bias.
- Referenz-Einordnung ausschliesslich als Band; die alte
  Einzelanker-Lesart bleibt nicht mergefaehig (Doc-Drift-Gate).

## 7. Naechste Stufe (nach PHY042)

1. Falls Paar-Trend weiter steigend: WM-Log-Fit/HKS-Formen direkt auf den
   glatten WL-Kurven (statt Paaren) pruefen, bevor groessere L bezahlt werden.
2. Falls Paar-Trend dreht: L=64 als Bestaetigungslauf.
3. Kagome-Uebertrag des WL-Kernels (PHY033-Nachfolge) erst nach
   Honeycomb-Abschluss.
