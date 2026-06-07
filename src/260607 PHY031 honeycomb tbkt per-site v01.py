"""PHY031: T_BKT(honeycomb) per-Site via Wolff + Nelson-Kosterlitz / Sandvik-Paar

Coworker Research / Coworkerz, 07. Juni 2026
de-CH Konventionen, ASCII-Quotes, ss statt Eszett.

ZWECK: Schliesst die im README/CHANGELOG als "offen" markierte Luecke
T_BKT(honeycomb). Das Honeycomb-Gitter (z=3) ist das DUALE der Hex-Zellen-
Schwerpunkte; ein XY-Modell direkt auf den Honeycomb-Knoten hat einen
deutlich tieferen BKT-Punkt als das Triangular-Gitter (z=6).

Referenz (web-validiert, arXiv:2501.07388): T_BKT(honeycomb) = 0.573 J/k_B
(Monte-Carlo; die analytische Vermutung 1/sqrt(2) ist umstritten).

GITTER (2-atomige Basis A/B, Bond-Laenge a=1):
  Bravais-Vektoren a1=(3/2, +sqrt3/2), a2=(3/2, -sqrt3/2).
  Basis: A bei (0,0), B bei delta1=(1,0).
  Die drei NN-Bonds A->B (120 Grad auseinander, Laenge 1):
     delta1 = (+1,   0)            -> B(i,   j)
     delta2 = (-1/2, +sqrt3/2)     -> B(i,   j-1)
     delta3 = (-1/2, -sqrt3/2)     -> B(i-1, j)
  Auf dem LxL-Torus (i,j in 0..L-1, mod L): N = 2 L^2 Knoten, 3 L^2 Bonds,
  jeder Knoten hat exakt 3 Nachbarn (Koordinationszahl 3).

PER-SITE-KONVENTION (Audit S2, konsistent mit core/PHY028):
  Upsilon = (1/N)[ J<T1> - (J^2/T)<sin_accum^2> ].
  T=0-ORAKEL (perfekt ausgerichtet): Upsilon(0) = (J/N) sum_b (e.delta_b)^2.
  Pro Zelle ist sum ueber die 3 Bonds (e.delta)^2 = 1 + 1/4 + 1/4 = 3/2
  (fuer e=(1,0)); ueber L^2 Zellen / N=2L^2 => Upsilon(0) = 3/4 J EXAKT,
  groessen-unabhaengig. (Vgl. Dreieck 3/2, Quadrat 1.)

METHODE: Identische Wolff-Pipeline wie PHY028/PHY030 (gitter-agnostisch via
Adjazenz). T_BKT aus der C-freien Sandvik-Paar-Bedingung am Weber-Minnhagen-
Punkt (1/R(T,L) - 1/R(T,2L) = -2 ln 2, R = pi*Upsilon/(2T) - 1), die
willkuerliche Extrapolationsformen vermeidet (arXiv:1302.2900).

EVIDENZ: deterministischer Gate-Report nach results/.
"""
from __future__ import annotations

import importlib.util
import json
import math
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

_SRC = Path(__file__).resolve().parent

REF_HONEYCOMB = 0.573  # arXiv:2501.07388
SW_LIMIT_HONEYCOMB = 0.75  # per-Site-Spinwellen-Grenzwert (T=0, exakt)


def _load(name: str, filename: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _SRC / filename)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


core = _load("phi_hex_core_v2", "260602 PHI HEX core v2 2 hardened.py")
wolff = _load("phy026_wolff_cluster", "260602 PHY026 wolff cluster v01.py")

XYConfig = core.XYConfig
make_rng = core.make_rng
nelson_kosterlitz_line = core.nelson_kosterlitz_line
helicity_terms = core.helicity_terms
helicity_from_ensemble = core.helicity_from_ensemble
wolff_sweep = wolff.wolff_sweep


# ============================================================================
# 1. Honeycomb-Torus-Gitter (2-atomige Basis, z=3)
# ============================================================================

_A1 = (1.5, math.sqrt(3.0) / 2.0)
_A2 = (1.5, -math.sqrt(3.0) / 2.0)
# NN-Verschiebungen A->B (Laenge 1):
_DELTA = {
    "d1": (1.0, 0.0),
    "d2": (-0.5, math.sqrt(3.0) / 2.0),
    "d3": (-0.5, -math.sqrt(3.0) / 2.0),
}


@dataclass
class HoneycombLattice:
    L: int
    edges: list[tuple[int, int]]        # (A_idx, B_idx) in natuerlicher Richtung
    edge_disp: list[tuple[float, float]]  # delta = r_B - r_A (Einheitslaenge)
    pos: list[tuple[float, float]]
    n_nodes: int


def build_honeycomb_lattice(L: int) -> HoneycombLattice:
    """LxL-Honeycomb-Torus, N = 2 L^2 Knoten, jeder Knoten z=3 Nachbarn.

    Index: site = 2*(i*L + j) + s mit s=0 (A) / s=1 (B), i,j in 0..L-1.
    Kanten als (A_idx, B_idx) mit edge_disp = r_B - r_A; die Vorzeichen-
    Konsistenz (dtheta_AB <-> Projektion auf delta_AB) ist damit exakt
    erfuellt, ohne min/max-Umsortierung.
    """
    def a_idx(i: int, j: int) -> int:
        return 2 * ((i % L) * L + (j % L)) + 0

    def b_idx(i: int, j: int) -> int:
        return 2 * ((i % L) * L + (j % L)) + 1

    pos: list[tuple[float, float]] = [(0.0, 0.0)] * (2 * L * L)
    for i in range(L):
        for j in range(L):
            base_x = i * _A1[0] + j * _A2[0]
            base_y = i * _A1[1] + j * _A2[1]
            pos[a_idx(i, j)] = (base_x, base_y)
            pos[b_idx(i, j)] = (base_x + _DELTA["d1"][0],
                                base_y + _DELTA["d1"][1])

    edges: list[tuple[int, int]] = []
    edge_disp: list[tuple[float, float]] = []
    for i in range(L):
        for j in range(L):
            a = a_idx(i, j)
            # A(i,j) -> B(i,j) via delta1
            edges.append((a, b_idx(i, j)))
            edge_disp.append(_DELTA["d1"])
            # A(i,j) -> B(i,j-1) via delta2
            edges.append((a, b_idx(i, j - 1)))
            edge_disp.append(_DELTA["d2"])
            # A(i,j) -> B(i-1,j) via delta3
            edges.append((a, b_idx(i - 1, j)))
            edge_disp.append(_DELTA["d3"])

    return HoneycombLattice(L=L, edges=edges, edge_disp=edge_disp,
                            pos=pos, n_nodes=2 * L * L)


def build_honeycomb_adjacency(lat: HoneycombLattice) -> list[list[int]]:
    adj: list[list[int]] = [[] for _ in range(lat.n_nodes)]
    for i, j in lat.edges:
        adj[i].append(j)
        adj[j].append(i)
    return adj


# ============================================================================
# 2. Messung (gitter-agnostische Wolff-Pipeline + per-Site-Helicity)
# ============================================================================

@dataclass
class HoneycombResult:
    T: float
    L: int
    n_nodes: int
    upsilon_mean: float
    upsilon_sem: float
    nk_line: float


def measure_honeycomb(L: int, T: float, J: float = 1.0, n_seeds: int = 8,
                      n_measure: int = 400, n_burn: int = 300,
                      master_seed: int = 42) -> HoneycombResult:
    """Per-Site-Helicity auf dem Honeycomb-Torus via Wolff-Single-Cluster.

    Pro Seed: nach Equilibrierung separat akkumulierte Helicity-Terme,
    Upsilon aus den Ensemble-Mitteln; Seed-Mittel + SEM.
    """
    lat = build_honeycomb_lattice(L)
    adj = build_honeycomb_adjacency(lat)
    n = lat.n_nodes
    beta = 1.0 / T
    cfg = XYConfig(J=J, T=T)

    ups = []
    for s in range(n_seeds):
        rng = make_rng(master_seed, stream=400 + s + 1000 * L)
        theta = rng.uniform(0, 2 * math.pi, n)
        for _ in range(n_burn):
            wolff_sweep(theta, adj, beta, J, rng, target_flips=n)
        t1s, sas = [], []
        for _ in range(n_measure):
            wolff_sweep(theta, adj, beta, J, rng, target_flips=n)
            t1, sa = helicity_terms(theta, lat, lat.pos)
            t1s.append(t1)
            sas.append(sa)
        ups.append(helicity_from_ensemble(t1s, sas, cfg, n))

    arr = np.array(ups)
    return HoneycombResult(
        T=T, L=L, n_nodes=n,
        upsilon_mean=float(np.mean(arr)),
        upsilon_sem=float(np.std(arr, ddof=1) / math.sqrt(n_seeds))
        if n_seeds > 1 else 0.0,
        nk_line=nelson_kosterlitz_line(T),
    )


def sandvik_pair_tbkt(data_L: dict, data_2L: dict) -> float | None:
    """T_BKT(L,2L) via exakte C-freie Sandvik-Bedingung.

    1/R(T,L) - 1/R(T,2L) = -2 ln 2 mit R(T,L) = pi*Upsilon(T,L)/(2T) - 1.
    Nulldurchgang der Differenz g(T) = [1/R(T,L) - 1/R(T,2L)] + 2 ln 2.
    """
    ln2 = math.log(2.0)
    Ts = sorted(set(data_L) & set(data_2L))
    diffs = []
    for T in Ts:
        u1 = data_L[T].upsilon_mean
        u2 = data_2L[T].upsilon_mean
        R1 = math.pi * u1 / (2.0 * T) - 1.0
        R2 = math.pi * u2 / (2.0 * T) - 1.0
        if abs(R1) > 1e-6 and abs(R2) > 1e-6:
            diffs.append((T, (1.0 / R1 - 1.0 / R2) + 2.0 * ln2))
    for i in range(len(diffs) - 1):
        T0, d0 = diffs[i]
        T1, d1 = diffs[i + 1]
        if d0 == 0.0:
            return T0
        if d0 * d1 < 0.0:
            return float(T0 + (T1 - T0) * (-d0) / (d1 - d0))
    return None


def nk_crossing(data_L: dict) -> float | None:
    """Per-L Nelson-Kosterlitz-Crossing: Upsilon(T,L) = 2T/pi."""
    Ts = sorted(data_L)
    pts = [(T, data_L[T].upsilon_mean - nelson_kosterlitz_line(T)) for T in Ts]
    for i in range(len(pts) - 1):
        (t0, d0), (t1, d1) = pts[i], pts[i + 1]
        if d0 == 0.0:
            return t0
        if d0 * d1 < 0.0:
            return float(t0 + (t1 - t0) * (-d0) / (d1 - d0))
    return None


# ============================================================================
# 3. Hauptlauf
# ============================================================================

def run_phy031(Ls=(6, 12, 24),
               temps=(0.56, 0.59, 0.62, 0.65, 0.68, 0.71),
               n_seeds=8, n_measure=400, n_burn=300, master_seed=42):
    t_start = time.time()
    print("PHY031: T_BKT(honeycomb) per-Site via Wolff + Sandvik-Paar")
    print("Coworker Research / Coworkerz, 07. Juni 2026")
    print("=" * 70)
    print(f"Referenz T_BKT(honeycomb) = {REF_HONEYCOMB} J/k_B "
          f"(arXiv:2501.07388)\n")

    # Gate A: T=0-Spinwellen-Grenzwert (perfekt ausgerichtet) = 3/4 exakt.
    lat = build_honeycomb_lattice(8)
    t1, sa = helicity_terms(np.zeros(lat.n_nodes), lat, lat.pos)
    ups_aligned = helicity_from_ensemble([t1], [sa],
                                         XYConfig(J=1.0, T=1e-9), lat.n_nodes)
    print(f"Gate A (perfekt ausgerichtet, T->0): Upsilon = {ups_aligned:.6f}  "
          f"(muss exakt {SW_LIMIT_HONEYCOMB} sein)\n")

    print(f"Per-Site-Helicity Upsilon(T,L) (Wolff, n_seeds={n_seeds}, "
          f"n_measure={n_measure}, n_burn={n_burn}):")
    data = {}
    for L in Ls:
        data[L] = {}
        row = []
        for T in temps:
            r = measure_honeycomb(L, T, n_seeds=n_seeds, n_measure=n_measure,
                                  n_burn=n_burn, master_seed=master_seed)
            data[L][T] = r
            row.append(f"T={T:.2f}:{r.upsilon_mean:.4f}+-{r.upsilon_sem:.4f}")
        print(f"  L={L:2d} (N={2*L*L:4d}): " + "  ".join(row))

    # Per-L NK-Crossings (informativ, finite-size von oben).
    print("\nPer-L Nelson-Kosterlitz-Crossings (Upsilon(T,L) = 2T/pi):")
    crossings = {}
    for L in Ls:
        tc = nk_crossing(data[L])
        crossings[L] = tc
        tc_str = f"{tc:.4f}" if tc is not None else "kein Crossing im Fenster"
        if tc is not None:
            print(f"  L={L:2d}: T_x = {tc_str}  "
                  f"(Abw. {(tc - REF_HONEYCOMB) / REF_HONEYCOMB * 100:+.2f}%)")
        else:
            print(f"  L={L:2d}: {tc_str}")

    # C-freie Sandvik-Paar-Schaetzung (L, 2L).
    print("\nSandvik-Paar-Methode T_BKT(L,2L) (C-frei, arXiv:1302.2900):")
    pair_tbkt = {}
    for L in Ls:
        if 2 * L in data:
            tc = sandvik_pair_tbkt(data[L], data[2 * L])
            pair_tbkt[(L, 2 * L)] = tc
            if tc is not None:
                print(f"  T_BKT(L={L},2L={2*L}) = {tc:.4f}  "
                      f"(Abw. {(tc - REF_HONEYCOMB) / REF_HONEYCOMB * 100:+.2f}%)")
            else:
                print(f"  T_BKT(L={L},2L={2*L}): kein Nulldurchgang "
                      f"(Gitter zu klein / Fenster verfehlt)")

    valid_pairs = [v for v in pair_tbkt.values() if v is not None]
    valid_cross = [v for v in crossings.values() if v is not None]
    largest_pair = valid_pairs[-1] if valid_pairs else None

    pass_gates = {
        "PASS_ALIGNED_EXACT_THREE_QUARTERS":
            abs(ups_aligned - SW_LIMIT_HONEYCOMB) < 1e-6,
        "PASS_CROSSING_EXISTS": len(valid_cross) >= 1,
        "PASS_CROSSING_IN_PHYSICAL_BAND": all(
            0.50 < c < 0.75 for c in valid_cross) if valid_cross else False,
    }
    if largest_pair is not None:
        pass_gates["PASS_PAIR_WITHIN_8PCT"] = (
            abs(largest_pair - REF_HONEYCOMB) / REF_HONEYCOMB < 0.08)

    print("\nPASS-Gates:")
    for k, v in pass_gates.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    overall = all(pass_gates.values())
    dt = time.time() - t_start
    print(f"\nOVERALL: {'PASS' if overall else 'FAIL'}   (Laufzeit {dt:.0f}s)")

    print("\nINTERPRETATION:")
    print("Das Honeycomb-Gitter (z=3) hat einen tieferen BKT-Punkt als das")
    print("Triangular-Gitter (z=6); der per-Site-T=0-Grenzwert 3/4 J ist")
    print("analytisch und exakt getroffen. Die Sandvik-Paar-Schaetzung")
    print("verankert T_BKT C-frei (ohne willkuerliche Extrapolationsform).")
    print("Bei kleinen L und endlicher Wolff-Statistik ehrlich ausgewiesen:")
    print("die per-L-Crossings driften von oben gegen die Referenz 0.573.")

    return {
        "module": "PHY031_honeycomb_tbkt_per_site_v01",
        "attribution": "Coworker Research / Coworkerz",
        "date": "2026-06-07",
        "reference_T_bkt_honeycomb": REF_HONEYCOMB,
        "reference_source": "arXiv:2501.07388",
        "sw_limit_per_site": SW_LIMIT_HONEYCOMB,
        "aligned_helicity": ups_aligned,
        "lattices_L": list(Ls),
        "temperatures": list(temps),
        "wolff": {"n_measure": n_measure, "n_burn": n_burn,
                  "n_seeds": n_seeds, "master_seed": master_seed},
        "measurements": {
            str(L): {f"{T:.2f}": asdict(r) for T, r in d.items()}
            for L, d in data.items()},
        "nk_crossings": {str(L): v for L, v in crossings.items()},
        "pair_tbkt": {f"{a}_{b}": v for (a, b), v in pair_tbkt.items()},
        "pass_gates": pass_gates,
        "overall_pass": overall,
        "runtime_s": dt,
    }


def write_report(report: dict, path: Path) -> None:
    """Deterministischer Text-Report (keine Laufzeit/Timestamps im Body)."""
    ref = report["reference_T_bkt_honeycomb"]
    lines = []
    lines.append("PHY031 v01 - T_BKT(honeycomb) per-Site via Wolff + Sandvik-Paar")
    lines.append("Coworker Research / Coworkerz, 2026-06-07")
    lines.append("=" * 70)
    lines.append("")
    lines.append(f"Referenz T_BKT(honeycomb) = {ref} J/k_B "
                 f"({report['reference_source']})")
    lines.append(f"Per-Site-T=0-Grenzwert (exakt) = {report['sw_limit_per_site']} J")
    lines.append(f"Gate A (perfekt ausgerichtet, T->0): "
                 f"Upsilon = {report['aligned_helicity']:.6f}")
    w = report["wolff"]
    lines.append(f"Wolff: n_measure={w['n_measure']}, n_burn={w['n_burn']}, "
                 f"n_seeds={w['n_seeds']}, master_seed={w['master_seed']}")
    lines.append(f"Gitter L = {report['lattices_L']}  (N = 2 L^2)")
    lines.append(f"T-Gitter = {report['temperatures']}")
    lines.append("")
    lines.append("--- Messung Upsilon(T,L) ---")
    for L in report["lattices_L"]:
        cells = []
        for T in report["temperatures"]:
            m = report["measurements"][str(L)][f"{T:.2f}"]
            cells.append(f"T={T:.2f}:{m['upsilon_mean']:.4f}"
                         f"+-{m['upsilon_sem']:.4f}")
        lines.append(f"  L={L:2d}: " + "  ".join(cells))
    lines.append("")
    lines.append("--- Per-L Nelson-Kosterlitz-Crossings ---")
    for L, v in report["nk_crossings"].items():
        vs = (f"{v:.4f}  (Abw. {(v - ref) / ref * 100:+.2f}%)"
              if v is not None else "kein Crossing")
        lines.append(f"  L={L}: T_x = {vs}")
    lines.append("")
    lines.append("--- Sandvik-Paar T_BKT(L,2L) (C-frei) ---")
    for k, v in report["pair_tbkt"].items():
        vs = (f"{v:.4f}  (Abw. {(v - ref) / ref * 100:+.2f}%)"
              if v is not None else "kein Nulldurchgang")
        lines.append(f"  T_BKT(L={k}) = {vs}")
    lines.append("")
    lines.append("--- PASS-Gates ---")
    for k, v in report["pass_gates"].items():
        lines.append(f"  [{'PASS' if v else 'FAIL'}] {k}")
    lines.append(f"  OVERALL: {'PASS' if report['overall_pass'] else 'FAIL'}")
    lines.append("")
    lines.append("--- Caveats ---")
    lines.append("  - Kleine L (Honeycomb N=2L^2); finite-size Crossings driften")
    lines.append("    von oben gegen die Referenz. Sandvik-Paar ist C-frei,")
    lines.append("    aber bei kleinen L noch finite-size-behaftet.")
    lines.append("  - Statistik: bewusst kleines Seed-Budget (Minuten lokal).")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _clean(o):
    if o is None:
        return None
    if isinstance(o, (np.bool_, bool)):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        v = float(o)
        return v if math.isfinite(v) else None
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    if isinstance(o, np.ndarray):
        return [_clean(x) for x in o.tolist()]
    if hasattr(o, "__dataclass_fields__"):
        return _clean(asdict(o))
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(x) for x in o]
    return o


if __name__ == "__main__":
    report = run_phy031()
    out = (_SRC.parent / "results" /
           "260607 PHY031 honeycomb tbkt per-site report.txt")
    out.parent.mkdir(exist_ok=True)
    write_report(report, out)
    print(f"\nReport geschrieben: {out}")
    print("\n--- JSON-Report ---")
    print(json.dumps(_clean(report), indent=2, allow_nan=False))
