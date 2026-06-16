# Methoden-Spec: Vierte-Ordnung-Helicity-Modul als BKT-Ordnungsparameter (Phi-Hex)

> **Status:** Vertrags-Quelle (spec/). Repo-native, 2026-06-16. Ergaenzt die
> HKS-Extrapolations-Spec (`260612 … hks …`) um die EINZIGE im Repo bisher
> ungenutzte Best-Practice-Observable: den 4.-Ordnungs-Helicity-Modul. Code
> folgt dieser Spec (`src/260616 PHY039 …`).

## 1. Motivation (web-recherchiert + gegen-recherchiert, 2026-06-16)

PHY034–038 haben die 2.-Ordnungs-(Weber-Minnhagen-/Paar-)Methode fuer
honeycomb ausgereizt und ehrlich befunden: der L→∞-Limes **wandert** mit dem
L-Set (0.574→0.605→0.639→0.647), und mehr Statistik heilt das nicht (PHY038).
Befund: Methoden-/finite-size-Limit.

**Websuche** (Best-Practice für honeycomb-T_BKT) und **Gegen-Recherche**
(Validierung der Treffer an Primärquellen) ergeben übereinstimmend:

- Die dedizierte honeycomb-Referenz mit T_BKT = **0.576(3)** / 0.575(8)
  (arXiv:2406.12076 = Phys. Scr., IOP `1402-4896/add62a`; sowie PTEP 2024
  `103A02`) bestimmt T_BKT **mit der zweiten UND vierten Ordnung** des
  Helicity-Modulus als Ordnungsparameter — und mit **Wang-Landau-entropischem
  Sampling + Simulated Annealing**, *nicht* mit blossem Wolff.
- Der 4.-Ordnungs-Modul (Minnhagen & Kim, **Phys. Rev. B 67, 172509 (2003)**,
  cond-mat/0304226) zeigt am BKT-Punkt eine **Diskontinuität** (negativer, mit
  L vertiefender „Dip") — das direkte Diagnostikum des universellen Sprungs,
  ohne a-priori-Annahme über die Natur des Übergangs.

Damit ist der 4.-Ordnungs-Modul der konkrete, literaturgestützte Hebel, den
die Serie noch nicht probiert hatte.

## 2. Estimator (Erst-Prinzipien-Herleitung — arXiv-unabhängig)

Twist Δ entlang Richtung **e**: θ_i → θ_i + Δ·(**e**·**r**_i). Pro Bond b mit
a_b = **e**·**r**_ij, φ_b = θ_i − θ_j gilt E(Δ) = −J Σ_b cos(φ_b + a_b Δ).
Ableitungen bei Δ=0 (pro Konfiguration):

| | Wert | Bond-Summe |
|---|---|---|
| E₁ = dE/dΔ | J·s | s = Σ a sin φ |
| E₂ | J·c | c = Σ a² cos φ |
| E₃ | −J·s₃ | s₃ = Σ a³ sin φ |
| E₄ | −J·c₄ | c₄ = Σ a⁴ cos φ |

Freie Energie F(Δ) = −T ln Z(Δ); ψ = ln Z erzeugt die Kumulanten von X = −βE.
Faà-di-Bruno-Kumulanten-Faltung bei Δ=0 (ungerade Terme verschwinden im
thermischen Ensemble durch θ→−θ):

```
ψ⁽²⁾ = −β⟨E₂⟩ + β² κ₂(E₁)
ψ⁽⁴⁾ = −β⟨E₄⟩ + β²[4 κ(E₃,E₁) + 3 κ(E₂,E₂)] − 6 β³ κ(E₂,E₁,E₁) + β⁴ κ₄(E₁)
```

mit F⁽ⁿ⁾ = −T·ψ⁽ⁿ⁾. Pro Site: **Υ₂ = F⁽²⁾/N**, **Υ₄ = F⁽⁴⁾/N**.
Diagnostik: **N·Υ₄ = F⁽⁴⁾** (der mit L vertiefende negative Dip markiert T_BKT
— Minnhagen-Kim-Signatur). Υ₂ = F⁽²⁾/N reproduziert exakt das bestehende
per-Site-Helicity-Schema (Konvention 2026-06-04).

**Korrektheits-Garantie (statt Literatur-Formel-Abschrift):** ein numerisches
**Orakel** vergleicht die geschlossene Form F⁽⁴⁾ gegen die numerische 4.
Ableitung der *exakten* freien Energie F(Δ) = −T ln mean_k exp(−βE(Δ;θ_k)) auf
**einer fixen Stichprobe** — die Übereinstimmung ist exakt (Maschinengenauigkeit,
rel-Fehler ~1e-8) und unabhängig von jedem (im Sandbox-Netz blockierten)
arXiv-Volltext. Gate: `f4_estimator_oracle()` / `test_phy039_…`.

## 3. Ergebnis — Quadratgitter-Goldstandard (Selbst-Falsifikation)

Validierung am Gitter mit hochpräzise bekanntem T_BKT = **0.89290(5)**
(arXiv:1302.2900 / 2501.07388). Wolff, n_seeds=12, n_measure=1500, L∈{8,16,24}.

| Gate | Ergebnis |
|---|---|
| F4-Estimator-Orakel (closed vs numerisch) | rel-Fehler ~1e-8 — **PASS** |
| Ausgerichtet T→0: Υ₂ | 1.000000 (square z=4) — **PASS** |
| N·Υ₄-Dip vorhanden, mit L vertiefend | ja — **PASS** |

**Befund (ehrlich):** Der negative Dip von N·Υ₄ EXISTIERT und vertieft sich mit
L wie erwartet. ABER seine **Lage konvergiert unter Wolff nicht sauber**: L=16
trifft 0.8901 (−0.31 %), doch L=8 **und** L=24 schieben das argmin an den
Gitterrand (>0.95). Der L=16-Treffer ist also **nicht robust**, sondern
rausch-/finite-size-dominiert (Dip ~15–25 % verrauscht). Der 4.-Ordnungs-
Kumulant gewinnt **unter Wolff keine Präzision** gegenüber der 2.-Ordnungs-
Paar-Methode.

Das deckt sich exakt mit PHY036/038 (Limes wandert mit L; mehr Statistik heilt
nicht) **und** erklärt, warum die Referenz 2406.12076 T_BKT gerade **mit
Wang-Landau-entropischem Sampling + SA** bestimmt — der 4.-Ordnungs-Kumulant
braucht eine Sampling-Methode mit flacher Energiehistogramm-Abdeckung, nicht
blosses Cluster-Sampling.

## 4. Korrekte nächste Stufe (Vertrag)

Der einzige im Repo noch ungenutzte Hebel ist **nicht „mehr L / mehr Seeds"**
(PHY036/038 falsifiziert), sondern ein **anderer Sampler**:

1. **Wang-Landau / multikanonisch** (flaches Energiehistogramm) — wie in
   2406.12076. Erlaubt belastbare 4.-Ordnungs-Statistik + freie-Energie-
   basierte FSS bei lokal vertretbaren L.
2. Optional **Simulated Annealing** als unabhängiger Quercheck (2406.12076).

Bis dieser Sampler implementiert ist, bleibt die ehrlichste honeycomb-Aussage
der Roh-Paar-Stand **~0.59** (L≈48, +3 % über Literatur, innerhalb finite-size;
PHY032). PHY039 liefert die **validierte Observable** (F4-Estimator + Orakel)
für diesen Lauf — der MC-Sampler ist der fehlende Baustein.

## 5. Wiederverwendbarer API (single source of truth)

| Funktion | Zweck |
|---|---|
| `bond_aggregates` | c, s, s₃, c₄ pro Konfiguration (gitter-agnostisch über a_b) |
| `fourth_order_F` | F⁽⁴⁾ = N·Υ₄, volle verbundene Kumulanten-Form |
| `second_order_F` | F⁽²⁾ = N·Υ₂ (reduced = bestehende Helicity-Konvention) |
| `f4_estimator_oracle` | numerisches Korrektheits-Orakel (Maschinengenauigkeit) |
| `measure_square_y2y4`, `dip_minimum` | Mess-Treiber + Dip-Lokalisierung |

Ein neues Gitter ist ein dünner Treiber über `bond_aggregates` +
`fourth_order_F` (a_b aus den Bond-Projektionen), **kein** Reimplement.

## 6. Provenance / Referenzen

- **Methode:** Minnhagen & Kim, Phys. Rev. B 67, 172509 (2003)
  (cond-mat/0304226); Best-Practice-Anwendung honeycomb: arXiv:2406.12076
  (Phys. Scr., IOP 1402-4896/add62a; PTEP 2024 103A02).
- **Referenz-T_BKT:** square 0.89290(5), honeycomb 0.573 / 0.576(3)
  (arXiv:2501.07388, arXiv:2406.12076).
- **Code/Evidenz:** `src/260616 PHY039 …`; Gate-Log
  `results/260616 PHY039 … report.txt`; Tests `tests/test_phy039_…`.

## 7. Definition of Done

- [x] Estimator erst-prinzipiell hergeleitet **und** numerisch (Orakel) auf
      Maschinengenauigkeit bestätigt — arXiv-unabhängig.
- [x] Selbst-Falsifikation am Goldstandard (square 0.8929) ausgeführt.
- [x] Negativ-Befund (Dip-Lage nicht robust unter Wolff) als Bürger erster
      Klasse dokumentiert — kein erzwungener Pass.
- [x] Korrekte nächste Stufe (Wang-Landau) als Vertrag benannt.
