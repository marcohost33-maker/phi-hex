# PHY024-R1mini — Kontrollierter Mini-Main-Run

Status: exploratory mini-main, nicht finaler PASS/FAIL.

Parameter: sigma=0.05, lambda=0, seeds=4, steps=400, burn_in=150, dt=0.05.

Wichtig: helicity_proxy ist nur mean cos(theta_i-theta_j) ueber Kanten, nicht der echte Helicity Modulus.

| radius | nodes | J_theta | null | helicity_proxy | vortex_density | g1(d=1) | g1(last) | fit |
|---:|---:|---:|---|---:|---:|---:|---:|---|
| 2 | 19 | 0.0 | hex | -0.02691 | 0.3571 | -0.02691 | -0.04261 | none |
| 2 | 19 | 0.0 | shuffled | 0.09498 | nan | 0.09498 | -0.04419 | none |
| 2 | 19 | 0.1 | hex | 0.9169 | 0.01 | 0.9169 | 0.6857 | exponential |
| 2 | 19 | 0.1 | shuffled | 0.8034 | nan | 0.8034 | 0.4447 | exponential |
| 2 | 19 | 0.4 | hex | 0.9991 | 0 | 0.9991 | 0.9963 | exponential |
| 2 | 19 | 0.4 | shuffled | 0.9985 | nan | 0.9985 | 0.9949 | exponential |
| 2 | 19 | 0.8 | hex | 0.9996 | 0 | 0.9996 | 0.9992 | exponential |
| 2 | 19 | 0.8 | shuffled | 0.9996 | nan | 0.9996 | 0.9994 | exponential |
| 3 | 37 | 0.0 | hex | 0.0487 | 0.3279 | 0.0487 | -0.03633 | none |
| 3 | 37 | 0.0 | shuffled | -0.05054 | nan | -0.05054 | 0.02148 | none |
| 3 | 37 | 0.1 | hex | 0.8007 | 0.1495 | 0.8007 | -0.164 | exponential |
| 3 | 37 | 0.1 | shuffled | 0.7846 | nan | 0.7846 | 0.4081 | exponential |
| 3 | 37 | 0.4 | hex | 0.9185 | 0.08316 | 0.9185 | 0.1196 | exponential |
| 3 | 37 | 0.4 | shuffled | 0.9977 | nan | 0.9977 | 0.9932 | exponential |
| 3 | 37 | 0.8 | hex | 0.9371 | 0.06842 | 0.9371 | 0.3807 | exponential |
| 3 | 37 | 0.8 | shuffled | 0.9997 | nan | 0.9997 | 0.9994 | exponential |

## Erste Auswertung

1. J_theta=0 reproduziert die PHY023-Richtung: g1 ist nahe Null oder instabil.
2. J_theta>0 erzeugt starke lokale Phasenausrichtung.
3. Shuffled-Nullmodelle zeigen ebenfalls starke Synchronisation, besonders bei J_theta=0.4/0.8.
4. Deshalb ist die Reparatur mit XY-Kopplung noch kein Hex-spezifischer Nachweis.
5. Radius 2/3 mit 4 Seeds ist nur Mini-Smoke; ein echter Crossover-Claim braucht radius 5/8/12 und mehr Seeds.

## Konsequenz

PHY024 darf nach diesem Mini-Run maximal als REVIEW-POSITIVE fuer lokale XY-Synchronisation gelten, nicht als PASS fuer BKT-kompatible Hex-Ordnungsphysik.
