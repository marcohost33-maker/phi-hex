# PHY024-R0 — XY-Phi-Hex Reparatur- und Validierungssprint v0.1

Datum: 2026-06-01  
Status: R0-DESIGN, nicht Ergebnis  
Warnstufe: ORANGE -> Reparaturpfad  
Ziel: Testen, ob explizite XY-Phasenkopplung in Phi-Hex v2.0 BKT-kompatible Ordnungsindikatoren erzeugt.

## 1. Ausgangslage

PHY023 zeigte für Phi-Hex v1.x einen fundamentalen Negativbefund:
- keine robuste räumliche Phasenordnung,
- g1(r) verschwindet in den getesteten Noise-Bereichen,
- Defekt-Lifetime ist als Hauptmetrik ungeeignet,
- Backreaction/PID darf nicht mehr als validierter Steuermechanismus gelten.

## 2. Forschungsfrage

Erzeugt eine explizite Nachbar-Phasenkopplung

    dtheta_i/dt = omega_i
                  + J_theta * sum_j sin(theta_j - theta_i)
                  + lambda * Q_i
                  + sigma * noise_i

auf einem Hex-Patch einen reproduzierbaren, finite-size-stabilen, BKT-kompatiblen Ordnungsindikator?

## 3. Nicht-Ziele

- Kein Beweis realer Raumzeitstruktur.
- Kein abgeschlossener BKT-Nachweis.
- Keine Rehabilitation von PID/Lifetime als Hauptmetrik.
- Keine Vermischung von Topologie-Stack PHY009 mit dynamischer Phi-Hex-Ordnungsbildung.

## 4. Primäre Nullhypothesen

H0-1: Auch bei J_theta > 0 bleibt g1(r) statistisch nahe Null.  
H0-2: J_theta erzeugt nur lokale Synchronisation, aber kein skalenstabiles räumliches Regime.  
H0-3: Es gibt keinen robusten Dosis-Wirkungs-Effekt über J_theta.  
H0-4: Backreaction lambda verbessert nur Artefaktmetriken, nicht g1(r), Pairing oder finite-size-Kurven.  
H0-5: Ein degree-erhaltendes Random-Graph-Nullmodell zeigt gleiche Effekte wie das Hex-Gitter.

## 5. Primäre Metriken

1. g1(r) = mean cos(theta_i - theta_j) nach Graphdistanz r.
2. Fit-Vergleich: algebraischer Zerfall vs. exponentieller Zerfall.
3. Vortex-Dichte über Plaquette-Winding.
4. Vortex-Antivortex-Pairing.
5. Helicity-Analog: mean cos(theta_i - theta_j) über Kanten.
6. Binder-artiger Winkel-/Magnetisationsindikator.
7. Bootstrap-Konfidenzintervalle.
8. FDR-korrigierte p-Werte bei Parametergitter-Auswertung.

## 6. Minimal-Gates

PASS nur wenn:
- J_theta = 0 reproduziert PHY023-Baseline.
- J_theta > J_c erzeugt reproduzierbare g1(r)-Struktur.
- Effekt bleibt bei wachsender Patchgrösse sichtbar.
- Algebraischer Fit schlägt exponentiellen Fit in einem stabilen Fenster.
- Vortex-Antivortex-Pairing steigt im selben Fenster.
- Random-Graph-Nullmodell schwächt oder zerstört den Effekt.
- Bootstrap-CI und FDR-Kontrolle bestehen.

REVIEW wenn:
- Effekt nur auf kleinen Patches sichtbar ist.
- Nur g1(r), aber kein Vortex-Pairing reagiert.
- Backreaction hilft nur über Lifetime/Artefaktmetriken.
- Fits widersprüchlich sind.

FAIL wenn:
- g1(r) überall verschwindet.
- J_theta nur lokale Glättung erzeugt.
- kein Dosis-Wirkungs-Zusammenhang entsteht.
- Hex-Geometrie irrelevant bleibt.

## 7. Run-Matrix R0

Smoke:
- radius: 2, 3
- J_theta: 0, 0.1, 0.4, 0.8
- sigma: 0.01, 0.05, 0.1
- lambda: 0
- seeds: 8

Main:
- radius: 2, 3, 5, 8
- J_theta: 0, 0.025, 0.05, 0.1, 0.2, 0.4, 0.8, 1.2
- sigma: 0.005, 0.01, 0.02, 0.05, 0.1, 0.2
- lambda: 0, 0.1, 0.5, 1.0
- seeds: 64

Confirm:
- selected transition windows only
- seeds: 128
- step-halving test required

## 8. Claim-Ledger

Allowed:
- "XY-Phi-Hex v2.0 is a repair candidate."
- "R0 tests BKT-compatible indicators."
- "PHY023 remains the baseline negative result."

Forbidden:
- "BKT transition proven."
- "HQST proves spacetime structure."
- "PID stabilizes Phi-Hex."
- "Internal verification equals physical proof."
