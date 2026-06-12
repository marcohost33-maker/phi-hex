# Methoden-Spec: HKS-L→∞-Extrapolation des BKT-T_BKT (Phi-Hex)

> **Status:** Vertrags-Quelle (spec/). Repo-native, 2026-06-12. Fasst den
> wiederverwendbaren Apparat **und** die ehrliche Methoden-Bilanz der Serie
> PHY034–038 zusammen. Code folgt dieser Spec.

## 1. Zweck

Ein gitter-agnostischer, getesteter Apparat, um aus dem Helicity-Modulus
Υ(T, L) den BKT-Übergang T_BKT im thermodynamischen Limes zu schätzen — und
die **Grenzen** dieser Schätzung ehrlich auszuweisen.

## 2. Methode (drei Stufen)

1. **Weber–Minnhagen-Finite-Size-Form** am BKT-Punkt
   (Weber & Minnhagen, PRB 37, 5986(R) (1988)):

   `Υ(T_BKT, L) = (2 T_BKT/π) · (1 + 1/(2 ln L + C))`

   mit unbekanntem `C`. Per-Site-Normierung (Υ(0) = z/2·J, Konvention
   2026-06-04).

2. **C-eliminierte Paar-Schätzung** (Methode B, Sandvik):
   `1/R(T,L) = 2 ln L + C` mit `R = πΥ/(2T) − 1` ⟹
   `1/R(T,L1) − 1/R(T,L2) = −2 ln(L2/L1)` ist **C-frei**. Der Nulldurchgang
   über das T-Gitter liefert `T_BKT(L1,L2)`.

3. **HKS-L→∞-Extrapolation** (Hsieh, Kao & Sandvik,
   J. Stat. Mech. (2013) P09001, arXiv:1302.2900): `T_BKT(L1,L2)` trägt selbst
   noch eine sub-leading log-Drift und wird linear in
   `u = 1/(ln Lc)^p` (charakteristische Größe `Lc`, HKS-Standard `p = 2`,
   `Lc =` geom. Mittel) zum Intercept `u → 0` extrapoliert (OLS über ≥2 Paare).

**Unsicherheiten (immer beide ausweisen):**
- *Methoden-Systematik* = Spannweite über `(Lc-Modus × p)`-Varianten.
- *Statistische CI* = Seed-Bootstrap des Extrapolations-Werts (Perzentil) bzw.
  MC-Propagation der per-Paar Bootstrap-σ. Achtung: die Extrapolation
  **verstärkt** den Paar-Fehler (Hebel > 1).

## 3. Wiederverwendbarer API (single source of truth)

| Funktion | Modul (src) |
|---|---|
| `wm_form`, `tbkt_pair_C_eliminated`, `fit_tbkt_fixedT` | PHY030 v02 wm-logfit |
| `bootstrap_tbkt_over_seeds`, `jackknife_…`, `make_sandvik_pair_estimator` | PHY032 |
| `ols_intercept`, `char_size`, `log_variable`, `extrapolate_pairs`, `extrapolation_ci` | PHY034 |
| `measure_cube`, `pair_estimates`, `make_hks_extrap_estimator` (param. `ladder`/`pairs`/`t_grid`) | PHY035 |

Ein neues Gitter/Experiment ist ein dünner Treiber über diese Bausteine
(vgl. PHY036/037/038), **kein** Reimplement.

## 4. Ergebnisse (Stand 2026-06-12)

### 4.1 Gitter-übergreifend (PHY037, committed Paare)

| Gitter | L | L→∞ (p=2,geom) | vs Ref | bestes Paar | Verdikt |
|---|---|---|---|---|---|
| kagome | 12/24/36 | 0.8236 | −0.17 % (0.825) | 0.8377 | **sauber** |
| honeycomb | 12/24/48 | 0.5747 | −0.23 % (0.576) | 0.5917 | **sauber, nur mit L=48** |
| triangular | 9/13/19 | 1.4591 | +2.90 % (1.418) | 1.3888 | **überschießt** |

### 4.2 Honeycomb — der Limes wandert mit dem L-Set

| Lauf | L-Set | T_BKT(L→∞) |
|---|---|---|
| PHY034 (2-Punkt) | 12/24/48 | 0.574 |
| PHY035 | 8…32 (kein 48) | 0.605 |
| PHY036 | 16…64 (8 Seeds) | 0.639 |
| PHY038 | 16…64 (16 Seeds, 240 Sweeps) | 0.647 |

**Schlüssel-Befund PHY038:** Verdoppelte Statistik (16 statt 8 Seeds, 240 statt
140 Sweeps) ändert das (32,64)-Paar nicht (0.6181 vs 0.6178) und verengt sein
CI nicht. Der hohe (32,64)-Wert ist also **kein Rausch-Artefakt** — die
honeycomb-Nicht-Auflösbarkeit ist ein **Methoden-/finite-size-Limit, kein
Statistik-Problem**. Mehr Seeds/Sweeps lösen es nicht; nötig wären größere L
(≥96..192) **und** ein feineres T-Gitter nahe dem Crossing.

## 5. Ehrliche Bilanz (Vertrag für künftige Nutzung)

- Die Klein-L-Paar-Extrapolation ist **gitter-abhängig** und **empfindlich auf
  das größte enthaltene L**. Sie ist **kein universelles Wundermittel**.
- Sie liefert SAUBERE Werte, wenn die Paare monoton fallen (kagome; honeycomb
  *mit* L=48). Sie **überschießt** bei sehr kleinem L (triangular ≤19) oder
  wenn das größte L rausch-dominiert ist (honeycomb L=64 bei 8 Seeds).
- Das **blanke größte Paar** (`T_BKT(L_max, 2·L_max)`) ist eine
  konkurrenzfähige, konservative Schranke und im Zweifel vorzuziehen.
- **Honeycomb-T_BKT** ist aus der lokal vertretbaren MC-Statistik **nicht
  robust extrapolierbar** (Limes wandert 0.574→0.605→0.639). Ehrlichster
  Stand: Roh-Paar ~0.59 (L≈48), +3 % über Literatur, innerhalb finite-size.
- Belastbare Konvergenz braucht **L ≥ 48..192 mit hoher Seed-Statistik**
  (HKS-Erfahrung), jenseits des lokalen Budgets.

## 6. Provenance / Referenzen

- Methodik: arXiv:1302.2900 (HKS, J. Stat. Mech. (2013) P09001);
  Weber & Minnhagen, PRB 37, 5986(R) (1988).
- Referenz-T_BKT: arXiv:2501.07388 (multi-lattice MC: square 0.89290(5),
  triangular 1.418, honeycomb 0.573, kagome 0.825 „rough estimate");
  arXiv:2406.12076 (dedizierte honeycomb-Studie, 0.576(3)).
- Code/Evidenz: `src/260611 PHY034…`, `…PHY035…`, `…PHY036…`,
  `…PHY037…`, `src/260612 PHY038…`; Gate-Logs in `results/`.

## 7. Definition of Done

- [x] Apparat gitter-agnostisch + getestet (Drift-Guards, Synthetik-Orakel,
      Wiring, Geometrie-Orakel).
- [x] Ergebnisse + Unsicherheiten **doppelt** (Systematik + statistisch).
- [x] Negativ-Befunde (PHY035/036, wandernder Limes) als Bürger erster Klasse.
- [x] Keine erzwungenen Passes; OVERALL=PASS = Analyse-Integrität.
