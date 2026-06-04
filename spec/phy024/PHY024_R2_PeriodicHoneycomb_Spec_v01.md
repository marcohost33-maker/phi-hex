# PHY024-R2 — Periodic Honeycomb XY-PhiHex Validation v0.1

Datum: 2026-06-01  
Status: R2-METHOD-HARDENING, nicht finaler Ergebnisstand  
Scope: periodische Honeycomb-/Torus-Geometrie, Hex-vs-Nullmodell, Twist-/Vortex-/g1-Diagnostik

## 1. Zweck

R0/R1mini zeigte:
- J_theta repariert lokale Phasenausrichtung.
- Shuffled-null synchronisiert ebenfalls.
- Daher ist noch nicht belegt, dass die Hex-Geometrie eine eigene Ordnungsphysik erzeugt.

R2 verschiebt deshalb den Fokus von offenen Mini-Patches auf periodische Honeycomb-Geometrie.

## 2. Harte Korrektur

In R0/R1 war `helicity_analog = mean cos(theta_i-theta_j)` nur ein lokaler Kohärenzproxy.
R2 führt zusätzlich einen Twist-Energy-Curvature-Proxy ein. Das ist noch nicht der vollständige
thermodynamische Helicity Modulus, aber methodisch näher an der Steifigkeitsfrage.

Echte Helicity-Modulus-Auswertung erfordert später:
- equilibrium Monte Carlo oder sauber kalibrierte stationäre Langevin-Dynamik,
- Temperatur-/Beta-Definition,
- freie Energie statt nur Energiekrümmung,
- grössere L und Finite-Size-Scaling.

## 3. Modell

theta_i dynamics:

    dtheta_i = [omega_i + J_theta * sum_j sin(theta_j - theta_i)] dt
               + sigma * sqrt(dt) dW_i

R2 setzt lambda=0. Backreaction bleibt deaktiviert, bis XY-only robust ist.

## 4. Periodische Honeycomb-Geometrie

Unit cell: (x,y,A/B), periodic in x/y.
Nodes: 2 * Lx * Ly.
Degree: 3.
Hex-plaquettes: one ordered six-site ring per unit cell.

## 5. Diagnostik

Primary:
- g1(r): graph-distance phase correlation
- vortex_density: hex-plaquette winding density
- vortex_pair_fraction: Anteil von Vortices mit Gegenladung im Abstand <= threshold
- twist_energy_curvature_proxy: finite-difference energy curvature under boundary twist
- hex-vs-degree-preserving shuffled null comparison

Secondary:
- local_edge_coherence
- magnetization magnitude
- algebraic/exponential fit diagnostic

## 6. Gates

PASS darf erst erwogen werden, wenn:
1. J_theta=0 bleibt ungeordnet.
2. J_theta>J_c erzeugt g1-Struktur.
3. Effekt unterscheidet sich von shuffled null.
4. Vortex-Dichte/Pairing reagiert konsistent.
5. Twist-energy-curvature-proxy steigt mit J_theta.
6. Trend bleibt für mehrere L erhalten.
7. FDR/Bootstrap-Auswertung wird in R3 ergänzt.

R2 kann maximal REVIEW-POSITIVE erzeugen, kein BKT-PASS.

## 7. Forbidden Claims

- "BKT bewiesen."
- "HQST beweist Raumzeit."
- "Twist-energy-curvature-proxy ist bereits echter Helicity Modulus."
- "Hex-Spezifität bewiesen", solange shuffled null ähnlich reagiert.
