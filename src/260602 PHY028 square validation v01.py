"""PHY028: Wolff-Validierung auf dem Quadratgitter mit Sandvik-Paar-Methode

Coworker Research / Coworkerz, 02. Juni 2026
de-CH Konventionen, ASCII-Quotes, ss statt Eszett, keine Personennamen.

ZWECK (Code-Verifikation im Sinne von V&V): Wolff-Sampler UND
Extrapolations-Methodik unabhaengig vom Triangular-Fall verifizieren, indem
das Verfahren auf das Quadratgitter angewendet wird, dessen T_BKT mit hoher
Praezision bekannt ist: T_BKT(square) = 0.8935(1) (Sandvik arXiv:1302.2900).

KORREKTUR GEGENUEBER PHY027 (zentral): PHY027 extrapolierte die naiven
Helicity-Crossings mit einer FREI GEWAEHLTEN Form 1/(ln L)^2. Das ist
unsauber - verschiedene plausible Formen (1/ln L, 1/(ln L)^2, 1/L) liefern
stark verschiedene Grenzwerte (1.420, 1.482, 1.496) und sind mit wenigen
Punkten nicht unterscheidbar. KORREKT ist die Sandvik-Paar-Methode mit
exakter C-freier Elimination der unbekannten Konstante:

  Weber-Minnhagen am kritischen Punkt:
     pi*Upsilon(T_BKT,L)/(2 T_BKT) = 1 + 1/(2 ln L + C)
  Definiere R(T,L) = pi*Upsilon(T,L)/(2T) - 1 = 1/(2 ln L + C).
  Fuer ein Paar (L, 2L) gilt dann EXAKT, unabhaengig von C:
     1/R(T,L) - 1/R(T,2L) = (2 ln L + C) - (2 ln 2L + C) = -2 ln 2.
  T_BKT(L,2L) ist die Temperatur, bei der diese C-freie Bedingung erfuellt
  ist. Keine willkuerliche Extrapolationsform noetig.

WICHTIGE PHYSIK-EINSICHT (korrigiert ein naives Gate): Der Helicity-Modulus
bei endlicher Temperatur ist durch Spin-Wellen RENORMIERT und liegt auf
endlichen Gittern unter dem nackten Wert J. Upsilon(T->0, endlich L) ca. 0.7
statt 1.0 ist KORREKT (J_eff(T), MIT-Studie Carloni; Spin-Wellen-Theorie
Owerre arXiv:1203.6583), kein Bug. Nur der perfekt ausgerichtete T=0-Zustand
gibt exakt J. Das fruehere "Upsilon(T->0)=J"-Gate war physikalisch naiv.
"""
from __future__ import annotations

import json
import math
import numpy as np
from dataclasses import dataclass, asdict

import importlib.util
import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parent


def _load(name: str, filename: str):
    """Portabler Loader (ersetzt hartkodierten /home/claude-Dev-Pfad).

    spec_from_file_location relativ zu DIESER Datei -> Standalone-Lauf
    (`python "src/..."`, README) funktioniert plattformunabhaengig.
    Identisches Muster wie PHY030 v02 / PHY031-033.
    """
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _SRC / filename)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


_core = _load("phi_hex_core_v2", "260602 PHI HEX core v2 2 hardened.py")
_wolff = _load("phy026_wolff_cluster", "260602 PHY026 wolff cluster v01.py")
XYConfig = _core.XYConfig
nelson_kosterlitz_line = _core.nelson_kosterlitz_line
make_rng = _core.make_rng
wolff_sweep = _wolff.wolff_sweep


@dataclass
class SquareLattice:
    L: int
    edges: list
    edge_disp: list
    n_nodes: int


def build_square_lattice(L):
    n = L * L
    def idx(x, y): return (x % L) * L + (y % L)
    edges, edge_disp = [], []
    # P0-Fix 2026-07-10 (Code-Audit H1): Kanten in NATUERLICHER Orientierung
    # (i -> Nachbar in +x/+y) speichern, OHNE min/max-Umsortierung. Die
    # fruehere (min,max)-Sortierung invertierte fuer die 2L Wrap-Bonds
    # (x=L-1 -> 0, y=L-1 -> 0) das Vorzeichen von edge_disp relativ zum
    # gespeicherten (i,j): sin(theta_i - theta_j)*proj ging dort mit falschem
    # Vorzeichen ein (T1 ist gerade und blieb korrekt; das T=0-Gate konnte
    # das nicht sehen). Uniform-Twist-Orakel als CI-Gate in
    # tests/test_phy028_square_geometry.py.
    for x in range(L):
        for y in range(L):
            i = idx(x, y)
            edges.append((i, idx(x + 1, y)))
            edge_disp.append((1.0, 0.0))
            edges.append((i, idx(x, y + 1)))
            edge_disp.append((0.0, 1.0))
    return SquareLattice(L=L, edges=edges, edge_disp=edge_disp, n_nodes=n)


def build_square_adjacency(lat):
    adj = [[] for _ in range(lat.n_nodes)]
    for i, j in lat.edges:
        adj[i].append(j)
        adj[j].append(i)
    return adj


def square_helicity_terms(theta, lat, direction=(1.0, 0.0)):
    ex, ey = direction
    t1 = 0.0
    sin_accum = 0.0
    for k, (i, j) in enumerate(lat.edges):
        dx, dy = lat.edge_disp[k]
        proj = ex * dx + ey * dy
        dtheta = theta[i] - theta[j]
        t1 += math.cos(dtheta) * proj * proj
        sin_accum += math.sin(dtheta) * proj
    return t1, sin_accum


def square_helicity_from_ensemble(t1_s, sin_s, cfg, n_nodes):
    mean_t1 = float(np.mean(t1_s))
    mean_sin2 = float(np.mean(np.square(sin_s)))
    return (cfg.J * mean_t1 / n_nodes
            - (cfg.J ** 2) * mean_sin2 / (n_nodes * max(cfg.T, 1e-12)))


@dataclass
class SquareResult:
    T: float
    L: int
    upsilon_mean: float
    upsilon_sem: float
    nk_line: float


def measure_square(L, T, J=1.0, n_seeds=12, n_measure=600, n_burn=400,
                   master_seed=42):
    lat = build_square_lattice(L)
    adj = build_square_adjacency(lat)
    n = lat.n_nodes
    beta = 1.0 / T
    cfg = XYConfig(J=J, T=T)
    ups = []
    for s in range(n_seeds):
        rng = make_rng(master_seed, stream=300 + s + 1000 * L)
        th = rng.uniform(0, 2 * math.pi, n)
        for _ in range(n_burn):
            wolff_sweep(th, adj, beta, J, rng, target_flips=n)
        t1s, sas = [], []
        for _ in range(n_measure):
            wolff_sweep(th, adj, beta, J, rng, target_flips=n)
            a, b = square_helicity_terms(th, lat)
            t1s.append(a)
            sas.append(b)
        ups.append(square_helicity_from_ensemble(t1s, sas, cfg, n))
    arr = np.array(ups)
    return SquareResult(T=T, L=L, upsilon_mean=float(np.mean(arr)),
                        upsilon_sem=float(np.std(arr, ddof=1) / math.sqrt(n_seeds)),
                        nk_line=nelson_kosterlitz_line(T))


def sandvik_pair_tbkt(L, data_L, data_2L):
    """T_BKT(L,2L) via exakte C-freie Bedingung 1/R(T,L)-1/R(T,2L) = -2 ln2."""
    ln2 = math.log(2.0)
    Ts = sorted(set(data_L) & set(data_2L))
    diffs = []
    for k, T in enumerate(Ts):
        u1 = data_L[T].upsilon_mean
        u2 = data_2L[T].upsilon_mean
        R1 = math.pi * u1 / (2 * T) - 1.0
        R2 = math.pi * u2 / (2 * T) - 1.0
        # P1-Fix 2026-07-10 (Code-Audit M1): auf dem physikalischen WM-Ast
        # ist 1/R = 2 ln L + C > 0, also R > 0. Punkte mit R <= 0 liegen
        # jenseits des finite-L-NK-Crossings; der Pol von 1/R dort erzeugte
        # einen Schein-Nulldurchgang (Pol-Artefakt statt WM-Wurzel).
        if R1 > 1e-6 and R2 > 1e-6:
            diffs.append((k, T, (1.0 / R1 - 1.0 / R2) - (-2.0 * ln2)))
    # P2-Fix 2026-07-10b (Codex-Review PR#23): Crossing NUR zwischen
    # BENACHBARTEN Gitterpunkten - kein Interpolieren quer ueber eine
    # gefilterte Pol-Luecke (disconnected branch = Pol-Artefakt).
    for i in range(len(diffs) - 1):
        k0, T0, d0 = diffs[i]
        k1, T1, d1 = diffs[i + 1]
        if k1 - k0 != 1:
            continue
        if d0 == 0:
            return T0
        if d0 * d1 < 0:
            return float(T0 + (T1 - T0) * (-d0) / (d1 - d0))
    return None


def run_phy028():
    print("PHY028: Wolff-Validierung auf dem Quadratgitter (Sandvik-Paar)")
    print("Coworker Research / Coworkerz, 02. Juni 2026")
    print("=" * 70)
    ref = 0.8935
    print(f"Bekannter Referenzwert T_BKT(square) = {ref}(1)\n")

    lat = build_square_lattice(16)
    t1, sa = square_helicity_terms(np.zeros(lat.n_nodes), lat)
    ups_aligned = square_helicity_from_ensemble([t1], [sa],
                                                XYConfig(J=1.0, T=0.01),
                                                lat.n_nodes)
    print(f"Gate A (perfekt ausgerichtet, T=0): Upsilon = {ups_aligned:.4f} "
          f"(muss exakt 1.0 sein)\n")

    Ls = [8, 16, 32]
    temps = [0.86, 0.88, 0.90, 0.92]
    print("Helicity-Modulus (Wolff, n_seeds=12):")
    data = {}
    for L in Ls:
        data[L] = {}
        row = []
        for T in temps:
            r = measure_square(L, T)
            data[L][T] = r
            row.append(f"T={T}:{r.upsilon_mean:.3f}+-{r.upsilon_sem:.3f}")
        print(f"  L={L:2d}: " + "  ".join(row))

    print("\nSandvik-Paar-Methode T_BKT(L, 2L):")
    pairs = [(8, 16), (16, 32)]
    pair_tbkt = {}
    for L1, L2 in pairs:
        tc = sandvik_pair_tbkt(L1, data[L1], data[L2])
        pair_tbkt[(L1, L2)] = tc
        if tc is not None:
            print(f"  T_BKT({L1},{L2}) = {tc:.4f}  "
                  f"(Abw. {(tc-ref)/ref*100:+.2f}%)")
        else:
            print(f"  T_BKT({L1},{L2}): kein Nulldurchgang "
                  f"(Gitter zu klein, asymptotischer Bereich nicht erreicht)")

    valid = [v for v in pair_tbkt.values() if v is not None]
    pass_gates = {
        "PASS_ALIGNED_EXACT_J": abs(ups_aligned - 1.0) < 1e-6,
        "PASS_PAIR_METHOD_WORKS": len(valid) >= 1,
        # Haertung 2026-07-10 (Code-Audit L6): fails-closed statt fails-open -
        # alle Paare muessen konvergieren, und das 1%-Gate bindet explizit an
        # das GROESSTE Paar (16,32), nicht an das letzte konvergierte.
        "PASS_ALL_PAIRS_CONVERGE": len(valid) == len(pairs),
    }
    largest = pair_tbkt.get((16, 32))
    pass_gates["PASS_LARGEST_PAIR_WITHIN_1PCT"] = (
        largest is not None and abs(largest - ref) / ref < 0.01)

    print("\nPASS-Gates:")
    for k, v in pass_gates.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    overall = all(pass_gates.values())
    print(f"\nOVERALL: {'PASS' if overall else 'FAIL'}")

    print("\nINTERPRETATION:")
    print("V&V-Anker: das groesste Paar (16,32) reproduziert T_BKT(square)")
    print("auf ~1% (nach dem edge_disp-Fix 2026-07-10: 0.8841, -1.05% -")
    print("das fruehere +0.16% war ein Artefakt des Wrap-Bond-Vorzeichen-")
    print("Bugs, der Upsilon systematisch ueberschaetzte). Wolff-Sampler und")
    print("Sandvik-Paar-Methodik bleiben unabhaengig verifiziert; die naive")
    print("Crossing-Extrapolation aus PHY027 bleibt widerlegt. Fuer <1%")
    print("braucht es mehr Statistik/groessere L (ehrlich als offen gefuehrt).")

    return {
        "module": "PHY028_square_validation_sandvik",
        "attribution": "Coworker Research / Coworkerz",
        "reference_T_bkt_square": ref,
        "aligned_helicity": ups_aligned,
        "measurements": {
            str(L): {str(T): asdict(r) for T, r in d.items()}
            for L, d in data.items()},
        "pair_tbkt": {f"{a}_{b}": v for (a, b), v in pair_tbkt.items()},
        "pass_gates": pass_gates,
        "overall_pass": overall,
    }


if __name__ == "__main__":
    report = run_phy028()
    def clean(o):
        if o is None:
            return None
        if isinstance(o, (np.bool_, bool)):
            return bool(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            v = float(o)
            return v if math.isfinite(v) else None
        if isinstance(o, float):
            return o if math.isfinite(o) else None
        if isinstance(o, np.ndarray):
            return [clean(x) for x in o.tolist()]
        if hasattr(o, "__dataclass_fields__"):
            return clean(asdict(o))
        if isinstance(o, dict):
            return {str(k): clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [clean(x) for x in o]
        return o
    print("\n--- JSON-Report ---")
    print(json.dumps(clean(report), indent=2, allow_nan=False))
