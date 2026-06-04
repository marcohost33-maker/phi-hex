# PHY024-R2 Runbook v0.1

## Ziel

R2 prüft, ob XY-PhiHex auf periodischer Honeycomb-Geometrie mehr zeigt als generische Graph-Synchronisation.

## Minimal-Sequenz

1. Baseline:
   python PHY024_R2_PeriodicHoneycomb_XY_PhiHex_v01.py --lx 4 --ly 4 --jtheta 0 --sigma 0.05 --seeds 8 --out r2_baseline.json

2. XY honeycomb:
   python PHY024_R2_PeriodicHoneycomb_XY_PhiHex_v01.py --lx 4 --ly 4 --jtheta 0.4 --sigma 0.05 --seeds 8 --out r2_xy_honeycomb.json

3. XY shuffled:
   python PHY024_R2_PeriodicHoneycomb_XY_PhiHex_v01.py --lx 4 --ly 4 --jtheta 0.4 --sigma 0.05 --seeds 8 --null shuffled --out r2_xy_shuffled.json

## R2 Review-Regeln

- Wenn honeycomb und shuffled beide stark synchronisieren, ist der Effekt generisch.
- Wenn honeycomb klare Plaquette-/Vortex-Struktur zeigt, shuffled aber nicht, ist das REVIEW-POSITIVE.
- Wenn J_theta=0 bereits Ordnung zeigt, ist Setup ungültig.
- Wenn twist_energy_curvature_proxy steigt, aber g1/vortex nicht reagieren, bleibt es REVIEW.
- Kein BKT-PASS in R2.
