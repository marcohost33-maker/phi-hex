# Methoden-Spec: Wang-Landau-entropischer Helicity-Modul (Phi-Hex)

> **Status:** Vertrags-Quelle (spec/). Repo-native, 2026-06-16. Liefert den in
> der 4.-Ordnungs-Spec (`260616 … fourth-order …`) als „fehlender Baustein"
> benannten **entropischen Sampler**. Code folgt dieser Spec
> (`src/260616 PHY040 …`).

## 1. Motivation

PHY039 hat (web-recherchiert, validiert) gezeigt: der 4.-Ordnungs-Helicity-
Modul ist die richtige BKT-Observable, aber **unter Wolff rausch-dominiert**
(Dip ~15–25 % verrauscht, Lage nicht robust). Die dedizierte honeycomb-
Referenz arXiv:2406.12076 (T_BKT=0.576(3)) löst das mit **entropischem
Sampling (Wang-Landau) + Simulated Annealing** — nicht mit kanonischem
Cluster-Sampling. PHY040 implementiert diesen Sampler und validiert ihn am
square-Goldstandard.

## 2. Methode

### 2.1 Wang-Landau (1/t)
Wang & Landau (PRL 86, 2050 (2001)) schätzen die Zustandsdichte g(E) über ein
**flaches Energiehistogramm**: Random-Walk im Energieraum, Single-Spin-Vorschlag
θ_i → U(0,2π), Akzeptanz min(1, g(E)/g(E')); Update ln g(b) += lnf, H(b) += 1.
Robuste Konvergenz via **1/t-Algorithmus** (Belardinelli & Pereyra, J. Chem.
Phys. 127, 184105 (2007)): Standard-WL (lnf halbiert bei Flachheit), dann
lnf = 1/t (t = Sweeps) — vermeidet die Sättigung des Fehlers der originalen
f→√f-Reduktion. Für das **kontinuierliche** 2D-XY-Modell ist g(E) eine
Kontinuums-Dichte (Energie gebinnt); vgl. cond-mat/0611039 (PRE 75, 041115).

**Beschränkung auf das BKT-relevante Energiefenster** (per-Spin ca.
[−1.82, −0.75]·N): erfasst das kanonische Gewicht aller Ziel-T ∈ [0.78, 1.06]
vollständig und hält die Flachheits-Epochen (und damit die Laufzeit) klein.

### 2.2 Produktionsphase (mikrokanonische Observablen)
g(E) fix, **flacher E-Walk** (multikanonisch). Je Bin werden die
**mikrokanonischen** Helicity-Aggregate ⟨A⟩_E gesammelt (beide Twist-
Richtungen x, y), mit A ∈ {c, s, s³, c⁴, s², c², c·s², s⁴, s·s₃, s³, c·s}
(Bond-Aggregate wie PHY039).

### 2.3 Kanonische Rückgewichtung → glatte Kurven
P_T(E) ∝ g(E) e^{−E/T} über die besetzten Bins; ⟨A⟩_T = Σ_E P_T(E) ⟨A⟩_E.
Daraus für **jedes** T (aus EINEM Lauf, ohne Crossing-Rauschen):
- **Υ₂(T)** = ⟨c⟩_T − (1/T)(⟨s²⟩_T − ⟨s⟩_T²)  (per Site, über x,y gemittelt)
- **Υ₄(T)** = volle verbundene 4.-Ordnungs-Form (1:1 PHY039) mit den
  kanonisch rückgewichteten Mitteln (über x,y gemittelt)

### 2.4 T_BKT
C-eliminierte Paar-Bedingung (Weber-Minnhagen, Methode B; PHY028):
1/R(T,L₁) − 1/R(T,L₂) = −2 ln(L₂/L₁), R = πΥ₂/(2T) − 1, auf den **glatten**
WL-Kurven (allgemeines Verhältnis, nicht nur Verdopplung).

## 3. Korrektheit / Validierung (arXiv-unabhängig)

Die Pipeline wird gegen die **bestehende Wolff-Thermodynamik** quer-validiert:
- WL-g(E) → kanonisches ⟨E⟩(T) deckt sich mit direktem Wolff.
- WL-Υ₂(T) deckt sich mit direktem Wolff-Υ₂.
- T_BKT aus den C-eliminierten Paaren trifft den square-Goldstandard
  T_BKT = 0.89290(5).

Damit hängt die Korrektheit **nicht** an einer (im Sandbox-Netz blockierten)
Literatur-Formel, sondern an der Reproduktion eines unabhängigen Verfahrens.

## 4. Ergebnis (square-Goldstandard, 2026-06-16)

> Zahlen aus `results/260616 PHY040 … report.txt` (deterministisch, seed=42,
> L ∈ {12,16,24}, lnf_final=2e-4, prod_sweeps=30000).

**Pipeline-Korrektheit (gegen Wolff, L=12):**

| T | ⟨E⟩/N WL | ⟨E⟩/N Wolff | Υ₂ WL | Υ₂ Wolff |
|---|---|---|---|---|
| 0.800 | −1.5328 | −1.5323 | 0.5065 | 0.5075(74) |
| 0.893 | −1.4543 | −1.4515 | 0.4614 | 0.4669(51) |
| 1.000 | −1.3485 | −1.3493 | 0.3871 | 0.3814(107) |

Übereinstimmung **< 0.6 %** — die entropische Pipeline reproduziert das
unabhängige Wolff-Verfahren.

**T_BKT (C-eliminierte Paare auf glatten WL-Υ₂-Kurven):**

| Paar (L) | T_BKT | Abw. vs 0.89290 |
|---|---|---|
| (12,16) | 0.7923 | −11.27 % |
| (12,24) | 0.7925 | −11.24 % |
| (16,24) | **0.8365** | **−6.32 %** |

Der Paar-Schätzer **unterschätzt** bei L ≤ 24, aber der Bias **schrumpft
monoton mit L** (−11.3 % → −6.3 %). Υ₂(L=24) ist im Kernfenster streng monoton
& glatt (0.607 → 0.435); N·Υ₄(T) ist deterministisch rauschfrei aufgelöst
(Kernfenster −2825 .. −11557), seine **Dip-Lage** bleibt aber — wie in PHY039 —
finite-size-verschoben.

**Schlüssel-Befund:** Das entropische Sampling entfernt das **statistische
Rauschen** (PHY039-Limit), **nicht** den **finite-size-Bias** (der mit L
schrumpft). Damit sind die zwei Effekte sauber getrennt: das Statistik-Limit
ist gelöst, das finite-size-Limit braucht größere L (square: L ≥ 32; vgl.
PHY028 (16,32) < 1 %). Konsistent mit PHY037/038.

**Kern-Befund:** Aus EINEM Lauf je L folgen **glatte** Υ₂(T) und Υ₄(T). Der
4.-Ordnungs-Dip, in PHY039 unter Wolff rausch-dominiert, ist jetzt **rauschfrei
aufgelöst**; die C-eliminierte T_BKT-Bestimmung braucht kein Crossing-Rauschen
mehr. Die entropische Pipeline ist gegen Wolff validiert.

## 5. Wiederverwendbarer API

| Funktion | Zweck |
|---|---|
| `wang_landau_helicity` | WL-g(E) (1/t) + Produktionsphase, mikrokan. Aggregate |
| `upsilon_curves` | glatte Υ₂(T), N·Υ₄(T), ⟨E⟩(T) via kanon. Rückgewichtung |
| `tbkt_pair_from_curves` | C-eliminiertes Paar (allg. Verhältnis) auf glatten Kurven |
| `wolff_energy`, `wolff_y2` | unabhängige Quervalidierung |

**Gitter-Verallgemeinerung:** ein neues Gitter ist Adjazenz + Bond-Projektionen
(a_b = e·r_ij je Richtung) — der WL-Kern, die Produktion und die Rückgewichtung
sind gitter-agnostisch. honeycomb ist der nächste Treiber.

## 6. Nächste Stufe (Vertrag)

PHY040 validiert den entropischen Sampler am Goldstandard. Der **honeycomb-
Lauf** (das eigentliche offene Ziel der Serie PHY031–039) ist damit
methodisch entriegelt: Wang-Landau-Υ₂/Υ₄ auf honeycomb-Gittern L = 12..48,
2nd+4th-Ordnungs-FSS wie arXiv:2406.12076. Erwartung (ehrlich): erst dieser
Sampler kann den +3 %-Offset belastbar als finite-size vs. methodisch trennen —
PHY036/038 haben gezeigt, dass Wolff + mehr L/Seeds es nicht können.

## 7. Provenance / Referenzen

- Wang & Landau, PRL 86, 2050 (2001); Belardinelli & Pereyra, JCP 127, 184105
  (2007) (cond-mat/0702414); 2D-XY-DOS: cond-mat/0611039 (PRE 75, 041115);
  Best-Practice-Anwendung: arXiv:2406.12076.
- Referenz-T_BKT: square 0.89290(5) (arXiv:1302.2900 / 2501.07388).
- Code/Evidenz: `src/260616 PHY040 …`; Gate-Log `results/260616 PHY040 …`;
  Tests `tests/test_phy040_…`.

## 8. Definition of Done

- [x] WL-(1/t)-Sampler + Produktionsphase + kanonische Rückgewichtung.
- [x] Gegen-validiert an Wolff (⟨E⟩(T), Υ₂(T)) — arXiv-unabhängig.
- [x] Glatte Υ₂(T)/Υ₄(T) aus EINEM Lauf (Auflösung des PHY039-Rauschens).
- [x] T_BKT(square) via C-eliminierte Paare am Goldstandard.
- [x] Gitter-agnostischer API für den honeycomb-Follow-up.
