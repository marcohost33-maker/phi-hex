# Erratum zu `results/260808 PHY043 triangular convention-free crossing report.txt`

Coworker Research / Coworkerz, 2026-09-26 — Issue #45 §3.

**Befund.** Der Grenzen-Abschnitt des Reports (Zeile 165) lautet

    - n_seeds=4, L<=19: keine 1%-Diskriminierung erwartbar (Spec §6).

Die Kopfzeile desselben Reports (Zeile 12) und die Spec §4 nennen fuer den
committeten Finallauf `n_seeds=8`, und die Gitterleiter reicht bis `L=25`
(`Gitter L = [9, 13, 19, 25]`). Zeile 165 war im Report-Generator
(`write_report`) hart kodiert und stammte aus dem Pilot-/PHY030-v02-Budget
(Spec §4, "Budget-Herleitung": radii 4/6/9, n_seeds=4), das vor dem
Evidenz-Pin superseded wurde.

**Korrekte Lesart.**

    - n_seeds=8, L<=25: keine 1%-Diskriminierung erwartbar (Spec §6).

**Reichweite.** Reiner Text-Drift. Messwerte, Gates (OVERALL PASS 6/6), die
Splay-/Crossing-Analyse und die Aussage NR-PHY043-01 (keine
1%-Diskriminierung zwischen 1.4007 und 1.418) sind unberuehrt; die
Grenz-Aussage selbst gilt bei 8 Seeds und L<=25 unveraendert (der Report
belegt sie mit Splay-Lagen 1.50..1.64).

**Behandlung.** Der Report bleibt byte-unveraendert (SHA-gepinnt in
`SOURCES.md`, AGENTS.md "Lineage ehrlich"). Der Generator leitet die Zeile
jetzt aus den Laufparametern ab (`n_seeds`, `max(lattices_L)`); ein Test
(`tests/test_phy043_convention_crossing.py`) bindet sowohl den Fix als auch
dieses Erratum an die Kopfzeile des committeten Reports.
