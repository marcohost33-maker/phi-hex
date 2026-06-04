# PHY024-R0 Runbook v0.1

## Reihenfolge

1. Baseline reproduzieren:
   - J_theta=0, lambda=0
   - sigma in [0.01, 0.05, 0.1]
   - Erwartung: g1(r) flach / nahe Null.

2. XY-only testen:
   - lambda=0
   - J_theta sweep: 0.025 bis 1.2
   - Primär: g1(r), Helicity-Analog, Vortex-Dichte.

3. Nullmodell prüfen:
   - gleiche Parameter auf degree-erhaltendem shuffled graph.
   - PASS nur, wenn Hex-Geometrie stärker/anders reagiert.

4. Backreaction sekundär:
   - Nur wenn XY-only Signal erzeugt.
   - lambda sweep: 0.1, 0.5, 1.0.
   - Kein PASS, wenn lambda nur Lifetime/Artefakte beeinflusst.

5. Confirm:
   - Transition-Fenster mit Seeds=128.
   - Step-Halving: dt halbieren und Signal prüfen.

## Minimalbefehle

Smoke baseline:

```bash
python PHY024_R0_XY_PhiHex_Reference_v01.py --radius 3 --jtheta 0 --sigma 0.05 --seeds 16 --out baseline.json
```

XY candidate:

```bash
python PHY024_R0_XY_PhiHex_Reference_v01.py --radius 3 --jtheta 0.4 --sigma 0.05 --seeds 16 --out xy_candidate.json
```

Shuffled null:

```bash
python PHY024_R0_XY_PhiHex_Reference_v01.py --radius 3 --jtheta 0.4 --sigma 0.05 --seeds 16 --null shuffled --out shuffled_null.json
```

## Review-Kriterien

- Wenn g1 nur auf radius=2/3 sichtbar ist: REVIEW.
- Wenn shuffled graph gleich gut ist: FAIL/REVIEW.
- Wenn power-fit besser ist, aber Vortex-Pairing nicht reagiert: REVIEW.
- Wenn J_theta=0 nicht den Negativbefund reproduziert: Setup ungültig.
