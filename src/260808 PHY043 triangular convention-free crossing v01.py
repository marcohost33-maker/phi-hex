"""PHY043: Konventionsfreier Quercheck triangular — Binder + xi_2/L

Coworker Research / Coworkerz, 08. August 2026
de-CH Konventionen, ASCII-Quotes, ss statt Eszett, keine Personennamen.

KONTEXT (Audit-Punkt O1, spec/260710 PHI HEX code audit v01.md §2):
Die per-Site-Normierung des Helicity-Modulus koennte auf triangular
(Flaeche/Site = sqrt(3)/2 < 1) das universelle Sprung-Crossing im Limes
nach UNTEN verschieben; der interne per-Site-Bestwert 1.4007 +/- 0.0081
(PHY030 v02) liegt -1.22 % unter der Referenz 1.418 (arXiv:2501.07388) —
konsistent mit genau so einem Restbias. O1 fordert einen KONVENTIONSFREIEN
Quercheck (Binder/Korrelations-Ratio) VOR jedem Code-Fix oder neuen
triangular-Claim. PHY043 ist dieser Quercheck.

METHODE (Vertragsquelle: spec/260808 PHI HEX phy043 triangular
convention-free crossing method v01.md):
Dimensionslose, normierungsfreie Observablen auf dem triangular-Torus,
gemessen mit der bestehenden Wolff-Pipeline (PHY026-Kern):

  U4(T, L)     = <|m|^4> / <|m|^2>^2,  m = (1/N) sum_j exp(i theta_j)
  xi_2/L(T, L) = sqrt(S(0)/S(k1) - 1) * sqrt(3) / (4 pi)

S(k) = (1/N) |sum_j exp(i theta_j - i k.r_j)|^2; k1 = Mittel ueber die
drei symmetrieaequivalenten Minimalmoden (n1,n2) in {(1,0),(0,1),(1,1)}
der torus-periodischen Wellen exp(2 pi i (n1 q + n2 r)/L)
(|k|^2 propto n1^2 + n2^2 - n1 n2 auf dem 60-Grad-Rhombus).

BKT-Signatur: Kurven mergen unterhalb T_BKT (kritische Phase) und splayen
oberhalb. Je Paar (L1 < L2) und Observable werden bestimmt:
  (1) Adjacent-Crossings von D(T) = O(T,L2) - O(T,L1) — Nulldurchgang NUR
      zwischen benachbarten Gitterpunkten (M1/P2-Vertrag des Root-Finder-
      Stacks), mit Flanken-Signifikanz max(|D|/sigma_D);
  (2) Splay-Temperatur T_splay — kleinstes Gitter-T, ab dem D fuer ALLE
      folgenden Punkte das Splay-Vorzeichen traegt und |D| > 2 sigma_D.
Unsicherheit: parametrischer Bootstrap (Punkt ~ Normal(mean, sem),
n_boot=300, stream=99043), Perzentil-CI + None-Quote (fails-closed).

BEWUSST NICHT: der universelle Wert (xi_2/L)* ~ 0.7507 (Hasenbusch 2005)
gilt fuer den quadratischen Torus (tau = i); unser Rhombus-Torus hat
tau = exp(i pi/3) — der Wert ist geometrieabhaengig und wird NICHT als
Anker benutzt.

RNG-VERTRAG (Audit O2/O3 umgesetzt): Neu-Produktion ohne Bit-Repro-
Altlast, daher T-Index im Stream: stream = 900 + s + 1000*L + 100000*t_idx
(Basis 900 frei; 800+ bleibt fuer die PHY039-Umbuchung reserviert).
T-Punkte eines Scans sind damit statistisch unabhaengig (kein CRN).

KEIN T_BKT-BESTWERT-CLAIM: Gates pruefen Pipeline-Integritaet, nicht
Physik-Wahrheit; der Physik-Befund ist ein FINDING (ggf. NR-PHY043-01).

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

REF_TRIANGULAR = 1.418          # arXiv:2501.07388
INTERNAL_PER_SITE = 1.4007      # PHY030 v02 WM-Fit (per-Site-Konvention)
SPLAY_Z = 2.0                   # Persistenz-Schwelle in sigma_D
N_BOOT = 300
BOOT_STREAM = 99043             # frei (99001 = PHY032-Bootstrap)
XI_PREFACTOR = math.sqrt(3.0) / (4.0 * math.pi)
MIN_MODES = ((1, 0), (0, 1), (1, 1))  # Norm-1-Moden; (1,-1) hat Norm 3


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

build_triangular_lattice = core.build_triangular_lattice
build_adjacency = core.build_adjacency
make_rng = core.make_rng
wolff_sweep = wolff.wolff_sweep


# ============================================================================
# Observablen-Kern (pure, testbar ohne Monte-Carlo)
# ============================================================================

def structure_factors(theta: np.ndarray, q_ax: np.ndarray, r_ax: np.ndarray,
                      L: int) -> tuple[float, float]:
    """(S(0), S(k1)) fuer eine Winkelkonfiguration auf dem LxL-Rhombus-Torus.

    S(k) = (1/N) |sum_j exp(i theta_j) exp(-2 pi i (n1 q_j + n2 r_j)/L)|^2.
    S(k1) ist das Mittel ueber die drei Minimalmoden MIN_MODES.
    """
    n = len(theta)
    spin = np.exp(1j * theta)
    s0 = float(np.abs(np.sum(spin)) ** 2) / n
    sk = 0.0
    for n1, n2 in MIN_MODES:
        phase = np.exp(-2j * math.pi * (n1 * q_ax + n2 * r_ax) / L)
        sk += float(np.abs(np.sum(spin * phase)) ** 2) / n
    return s0, sk / len(MIN_MODES)


def xi_over_L_from_S(s0: float, sk1: float) -> float | None:
    """xi_2/L = sqrt(S(0)/S(k1) - 1) * sqrt(3)/(4 pi); None fails-closed.

    None bei numerisch verschwindendem S(k1) (exakt 0 ist z.B. bei einer
    uniformen Konfiguration die Geometrie-Summe; Float-Roundoff liefert
    dann O(1e-30) statt 0 — ein Schwellwert RELATIV zu S(0) faengt beides)
    oder bei S(0) < S(k1) (Wurzel-Argument negativ — bei MC-Rauschen tief
    in der ungeordneten Phase moeglich; ein solcher Punkt ist kein valider
    xi-Schaetzer und wird NICHT auf 0 geclampt).
    """
    if sk1 <= 1e-12 * max(s0, 1.0) or s0 < sk1:
        return None
    return math.sqrt(s0 / sk1 - 1.0) * XI_PREFACTOR


def binder_u4(m2_samples: np.ndarray, m4_samples: np.ndarray) -> float:
    """U4 = <|m|^4> / <|m|^2>^2 aus Sample-Reihen der Momente."""
    m2 = float(np.mean(m2_samples))
    m4 = float(np.mean(m4_samples))
    return m4 / (m2 * m2)


# ============================================================================
# Crossing-/Splay-Finder (pure; Vertraege siehe Spec §3)
# ============================================================================

def adjacent_crossings(T_grid, d_vals, d_sigmas):
    """Nulldurchgaenge von D(T) NUR zwischen benachbarten Gitterpunkten.

    d_vals darf None-Punkte enthalten (fails-closed): ein Vorzeichenwechsel
    QUER ueber eine None-Luecke ist kein Crossing (M1/P2-Vertrag, gleiche
    Regel wie tbkt_pair_C_eliminated in PHY030 v02).

    Rueckgabe: Liste von dicts {T_cross, significance} mit
    significance = max(|D|/sigma) der beiden flankierenden Punkte
    (Crossings im Merge-Bereich koennen Rauschen sein — die Signifikanz
    macht das im Report sichtbar, statt sie zu verstecken).
    """
    out = []
    for i in range(len(T_grid) - 1):
        d0, d1 = d_vals[i], d_vals[i + 1]
        if d0 is None or d1 is None:
            continue
        t0, t1 = T_grid[i], T_grid[i + 1]
        if d0 == 0.0:
            sig = abs(d1) / max(d_sigmas[i + 1], 1e-300)
            out.append({"T_cross": float(t0), "significance": float(sig)})
            continue
        if d0 * d1 < 0.0:
            t_x = t0 + (t1 - t0) * (-d0) / (d1 - d0)
            sig = max(abs(d0) / max(d_sigmas[i], 1e-300),
                      abs(d1) / max(d_sigmas[i + 1], 1e-300))
            out.append({"T_cross": float(t_x), "significance": float(sig)})
    return out


def splay_temperature(T_grid, d_vals, d_sigmas, sign: float,
                      z: float = SPLAY_Z, min_points: int = 2):
    """Kleinstes Gitter-T, ab dem D persistent signifikant splayt.

    Vertrag (Spec §3): fuer ALLE Gitterpunkte T_j >= T_splay muss
    sign * D(T_j) > z * sigma_D(T_j) gelten, und der persistente Tail
    muss mindestens min_points = 2 Gitterpunkte umfassen — ein einzelner
    signifikanter Endpunkt ist kein Splay-Nachweis. None-Punkte brechen
    die Persistenz (fails-closed). Rueckgabe None, wenn kein solches T
    existiert.
    """
    n = len(T_grid)
    for k in range(n - min_points + 1):
        ok = True
        for j in range(k, n):
            d = d_vals[j]
            if d is None or sign * d <= z * d_sigmas[j]:
                ok = False
                break
        if ok:
            return float(T_grid[k])
    return None


def pair_analysis(T_grid, obs_L1, sem_L1, obs_L2, sem_L2, splay_sign):
    """Crossing + Splay fuer ein Paar (L1 < L2) einer Observable.

    obs_*: Kurven ueber T_grid (Eintraege duerfen None sein, fails-closed).
    splay_sign: erwartetes Vorzeichen von D = O(L2) - O(L1) OBERHALB des
    Uebergangs (-1 fuer xi_2/L, +1 fuer U4).
    """
    d_vals, d_sigmas = [], []
    for i in range(len(T_grid)):
        a, b = obs_L1[i], obs_L2[i]
        if a is None or b is None:
            d_vals.append(None)
            d_sigmas.append(float("inf"))
        else:
            d_vals.append(b - a)
            d_sigmas.append(math.hypot(sem_L1[i], sem_L2[i]))
    return {
        "crossings": adjacent_crossings(T_grid, d_vals, d_sigmas),
        "T_splay": splay_temperature(T_grid, d_vals, d_sigmas, splay_sign),
    }


def bootstrap_splay_ci(T_grid, obs_L1, sem_L1, obs_L2, sem_L2, splay_sign,
                       rng, n_boot: int = N_BOOT):
    """Parametrischer Bootstrap der Splay-Temperatur (Spec §3).

    Jeder Kurvenpunkt wird als Normal(mean, sem) resampelt; None-Punkte
    bleiben None. Rueckgabe: {ci_low, ci_high, none_frac} mit Perzentil-CI
    [2.5 %, 97.5 %] ueber die Nicht-None-Resamples; none_frac wird
    ausgewiesen statt verschwiegen (fails-closed-Transparenz).
    """
    vals = []
    n_none = 0
    for _ in range(n_boot):
        o1, o2 = [], []
        for i in range(len(T_grid)):
            a, b = obs_L1[i], obs_L2[i]
            o1.append(None if a is None
                      else a + sem_L1[i] * rng.standard_normal())
            o2.append(None if b is None
                      else b + sem_L2[i] * rng.standard_normal())
        res = pair_analysis(T_grid, o1, sem_L1, o2, sem_L2, splay_sign)
        if res["T_splay"] is None:
            n_none += 1
        else:
            vals.append(res["T_splay"])
    if vals:
        lo, hi = np.percentile(vals, [2.5, 97.5])
    else:
        lo = hi = float("nan")
    return {"ci_low": float(lo), "ci_high": float(hi),
            "none_frac": n_none / n_boot}


# ============================================================================
# Messung (Wolff-Pipeline, PHY026-Kern)
# ============================================================================

@dataclass
class RatioMeasurement:
    T: float
    L: int
    u4_mean: float
    u4_sem: float
    xi_ratio_mean: float | None
    xi_ratio_sem: float | None


def measure_ratios_wolff(radius: int, T: float, t_idx: int, J: float = 1.0,
                         n_measure: int = 400, n_burn: int = 300,
                         n_seeds: int = 4,
                         master_seed: int = 42) -> RatioMeasurement:
    """Misst U4 und xi_2/L per Seed; Seed-Mittel + SEM.

    Stream-Vertrag: 900 + s + 1000*L + 100000*t_idx (Header/Spec §4).
    xi_2/L wird je Seed aus den Ensemble-Mitteln <S(0)>, <S(k1)> gebildet
    (Verhaeltnis der Mittel, nicht Mittel der Verhaeltnisse — gleiche
    Logik wie helicity_from_ensemble). Liefert ein Seed keinen validen
    xi-Schaetzer (None), ist der Punkt insgesamt None (fails-closed).
    """
    lattice = build_triangular_lattice(radius=radius, periodic=True)
    adj = build_adjacency(lattice)
    n = lattice.n_nodes
    L = 2 * radius + 1
    q_ax = np.array([a[0] for a in lattice.ax_of], dtype=float)
    r_ax = np.array([a[1] for a in lattice.ax_of], dtype=float)
    beta = 1.0 / T

    u4_seeds, xi_seeds = [], []
    for s in range(n_seeds):
        rng = make_rng(master_seed,
                       stream=900 + s + 1000 * L + 100000 * t_idx)
        theta = rng.uniform(0, 2 * math.pi, n)
        for _ in range(n_burn):
            wolff_sweep(theta, adj, beta, J, rng, target_flips=n)
        m2_s, m4_s, s0_s, sk_s = [], [], [], []
        for _ in range(n_measure):
            wolff_sweep(theta, adj, beta, J, rng, target_flips=n)
            s0, sk1 = structure_factors(theta, q_ax, r_ax, L)
            m2 = s0 / n          # |m|^2 = S(0)/N
            m2_s.append(m2)
            m4_s.append(m2 * m2)
            s0_s.append(s0)
            sk_s.append(sk1)
        u4_seeds.append(binder_u4(np.array(m2_s), np.array(m4_s)))
        xi_seeds.append(xi_over_L_from_S(float(np.mean(s0_s)),
                                         float(np.mean(sk_s))))

    u4_arr = np.array(u4_seeds)
    u4_sem = (float(np.std(u4_arr, ddof=1) / math.sqrt(n_seeds))
              if n_seeds > 1 else 0.0)
    if any(x is None for x in xi_seeds):
        xi_mean, xi_sem = None, None
    else:
        xi_arr = np.array(xi_seeds, dtype=float)
        xi_mean = float(np.mean(xi_arr))
        xi_sem = (float(np.std(xi_arr, ddof=1) / math.sqrt(n_seeds))
                  if n_seeds > 1 else 0.0)
    return RatioMeasurement(T=T, L=L, u4_mean=float(np.mean(u4_arr)),
                            u4_sem=u4_sem, xi_ratio_mean=xi_mean,
                            xi_ratio_sem=xi_sem)


# ============================================================================
# Hauptlauf
# ============================================================================

def run_phy043(radii=(4, 6, 9, 12),
               temps=(1.36, 1.38, 1.40, 1.42, 1.44, 1.46, 1.48, 1.50,
                      1.52, 1.54, 1.56, 1.58, 1.60, 1.62, 1.64, 1.66,
                      1.68, 1.70),
               n_measure=800, n_burn=300, n_seeds=8, master_seed=42):
    t_start = time.time()
    print("PHY043: Konventionsfreier Quercheck triangular (U4 + xi_2/L)")
    print("Coworker Research / Coworkerz, 08. August 2026")
    print("=" * 70)
    print(f"Referenz T_BKT(triangular) = {REF_TRIANGULAR} (arXiv:2501.07388)")
    print(f"Interner per-Site-Wert (PHY030 v02) = {INTERNAL_PER_SITE}")
    print("Observablen: U4 = <|m|^4>/<|m|^2>^2;  "
          "xi_2/L = sqrt(S0/Sk1 - 1)*sqrt(3)/(4 pi)\n")

    Ls = [2 * r + 1 for r in radii]
    T_grid = list(temps)
    print(f"Gitter L = {Ls}  T-Gitter (Delta=0.02) = {T_grid}")
    print(f"Wolff: n_measure={n_measure}, n_burn={n_burn}, "
          f"n_seeds={n_seeds}, master_seed={master_seed}")
    print("Stream-Vertrag: 900 + s + 1000*L + 100000*t_idx "
          "(T-unabhaengige Streams, Audit O2)\n")

    meas: dict[int, dict[float, RatioMeasurement]] = {L: {} for L in Ls}
    for r in radii:
        L = 2 * r + 1
        cells = []
        for t_idx, T in enumerate(T_grid):
            m = measure_ratios_wolff(radius=r, T=T, t_idx=t_idx,
                                     n_measure=n_measure, n_burn=n_burn,
                                     n_seeds=n_seeds,
                                     master_seed=master_seed)
            meas[L][T] = m
            xi_txt = ("None" if m.xi_ratio_mean is None
                      else f"{m.xi_ratio_mean:.4f}")
            cells.append(f"T={T:.2f}:U4={m.u4_mean:.4f},xi/L={xi_txt}")
        print(f"  L={L:2d}: " + "  ".join(cells))

    def curve(L, field_mean, field_sem):
        vals = [getattr(meas[L][T], field_mean) for T in T_grid]
        sems = [getattr(meas[L][T], field_sem) for T in T_grid]
        sems = [0.0 if s is None else s for s in sems]
        return vals, sems

    pairs = [(Ls[i], Ls[j]) for i in range(len(Ls))
             for j in range(i + 1, len(Ls))]
    analysis = {"xi_ratio": {}, "u4": {}}
    boot_rng = make_rng(master_seed, stream=BOOT_STREAM)
    for L1, L2 in pairs:
        for obs, mean_f, sem_f, sign in (
                ("xi_ratio", "xi_ratio_mean", "xi_ratio_sem", -1.0),
                ("u4", "u4_mean", "u4_sem", +1.0)):
            o1, s1 = curve(L1, mean_f, sem_f)
            o2, s2 = curve(L2, mean_f, sem_f)
            res = pair_analysis(T_grid, o1, s1, o2, s2, sign)
            res["splay_ci"] = bootstrap_splay_ci(
                T_grid, o1, s1, o2, s2, sign, boot_rng)
            analysis[obs][f"({L1},{L2})"] = res

    print("\n--- Paar-Analyse (D = O(L2) - O(L1)) ---")
    for obs in ("xi_ratio", "u4"):
        print(f"  Observable {obs}:")
        for key, res in analysis[obs].items():
            cr = ", ".join(f"{c['T_cross']:.4f} (sig {c['significance']:.1f})"
                           for c in res["crossings"]) or "-"
            ts = res["T_splay"]
            ci = res["splay_ci"]
            ts_txt = "None" if ts is None else f"{ts:.2f}"
            print(f"    {key}: T_splay={ts_txt} "
                  f"CI[{ci['ci_low']:.3f},{ci['ci_high']:.3f}] "
                  f"none_frac={ci['none_frac']:.2f}  Crossings: {cr}")

    # --- Gates (Pipeline-Integritaet, Spec §5) ---
    T_min, T_max = T_grid[0], T_grid[-1]
    all_points = [meas[L][T] for L in Ls for T in T_grid]
    input_complete = all(
        m.xi_ratio_mean is not None and math.isfinite(m.u4_mean)
        and math.isfinite(m.xi_ratio_mean) for m in all_points)
    u4_range = all(0.9 <= m.u4_mean <= 2.1 for m in all_points)

    def xi_at(L, T):
        return meas[L][T].xi_ratio_mean

    high_t_xi = all(
        xi_at(Ls[i], T_max) is not None and xi_at(Ls[i + 1], T_max) is not None
        and xi_at(Ls[i], T_max) > xi_at(Ls[i + 1], T_max)
        for i in range(len(Ls) - 1))
    high_t_u4 = all(
        meas[Ls[i]][T_max].u4_mean < meas[Ls[i + 1]][T_max].u4_mean
        for i in range(len(Ls) - 1))

    low_t_merge = True
    for i in range(len(Ls)):
        for j in range(i + 1, len(Ls)):
            a, b = meas[Ls[i]][T_min], meas[Ls[j]][T_min]
            if a.xi_ratio_mean is None or b.xi_ratio_mean is None:
                low_t_merge = False
                continue
            diff = abs(a.xi_ratio_mean - b.xi_ratio_mean)
            sig3 = 3.0 * math.hypot(a.xi_ratio_sem or 0.0,
                                    b.xi_ratio_sem or 0.0)
            rel8 = 0.08 * max(a.xi_ratio_mean, b.xi_ratio_mean)
            if diff >= max(rel8, sig3):
                low_t_merge = False

    largest_pair = f"({Ls[-2]},{Ls[-1]})"
    splay_exists = analysis["xi_ratio"][largest_pair]["T_splay"] is not None

    pass_gates = {
        "PASS_INPUT_COMPLETE": bool(input_complete),
        "PASS_U4_RANGE": bool(u4_range),
        "PASS_HIGH_T_ORDERING_XI": bool(high_t_xi),
        "PASS_HIGH_T_ORDERING_U4": bool(high_t_u4),
        "PASS_LOW_T_MERGE_XI": bool(low_t_merge),
        "PASS_SPLAY_XI_LARGEST_PAIR": bool(splay_exists),
    }

    print("\n--- PASS-Gates (Pipeline-Integritaet) ---")
    for k, v in pass_gates.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    overall = all(pass_gates.values())
    dt = time.time() - t_start
    print(f"\nOVERALL: {'PASS' if overall else 'FAIL'}   (Laufzeit {dt:.0f}s)")

    # --- FINDING (Diagnostik, kein Gate; Interpretations-Vertrag Spec §6) ---
    ts_largest = analysis["xi_ratio"][largest_pair]["T_splay"]
    print("\nFINDING (kein Gate, kein Bestwert-Claim):")
    if ts_largest is None:
        print("  Kein persistenter xi_2/L-Splay im Fenster — nicht "
              "diskriminierend (NR-PHY043-01).")
    else:
        print(f"  xi_2/L-Splay {largest_pair} beginnt bei T={ts_largest:.2f}; "
              f"Einordnung relativ zu intern {INTERNAL_PER_SITE} und "
              f"Referenz {REF_TRIANGULAR} siehe Report.")

    return {
        "module": "PHY043_triangular_convention_free_crossing_v01",
        "attribution": "Coworker Research / Coworkerz",
        "date": "2026-08-08",
        "method": "Binder U4 + second-moment xi_2/L auf triangular-Torus; "
                  "adjacent-Crossings + persistente Splay-Temperatur "
                  "(z=2), parametrischer Bootstrap",
        "spec": "spec/260808 PHI HEX phy043 triangular convention-free "
                "crossing method v01.md",
        "reference_T_bkt_triangular": REF_TRIANGULAR,
        "internal_per_site_T_bkt": INTERNAL_PER_SITE,
        "lattices_L": Ls,
        "temperatures": T_grid,
        "wolff": {"n_measure": n_measure, "n_burn": n_burn,
                  "n_seeds": n_seeds, "master_seed": master_seed,
                  "stream_contract": "900 + s + 1000*L + 100000*t_idx"},
        "bootstrap": {"n_boot": N_BOOT, "stream": BOOT_STREAM},
        "measurements": {
            str(L): {f"{T:.2f}": asdict(meas[L][T]) for T in T_grid}
            for L in Ls},
        "pair_analysis": analysis,
        "pass_gates": pass_gates,
        "overall_pass": bool(overall),
        "runtime_s": dt,
    }


# ============================================================================
# Report
# ============================================================================

def _clean(o):
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, float) and not math.isfinite(o):
        return str(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.integer):
        return int(o)
    return o


def write_report(report: dict, path: Path) -> None:
    lines = []
    lines.append("PHY043 - Konventionsfreier Quercheck triangular "
                 "(U4 + xi_2/L)")
    lines.append("Coworker Research / Coworkerz, 2026-08-08")
    lines.append("=" * 70)
    lines.append("")
    lines.append("Zweck: Audit-Punkt O1 (spec/260710 §2) — liegt die "
                 "konventionsfreie")
    lines.append("BKT-Signatur beim internen per-Site-Wert 1.4007 oder bei "
                 "der Referenz 1.418?")
    lines.append(f"Spec (Vertragsquelle): {report['spec']}")
    lines.append(f"Methode: {report['method']}")
    lines.append("")
    lines.append(f"Gitter L = {report['lattices_L']}")
    lines.append(f"T-Gitter = {report['temperatures']}")
    w = report["wolff"]
    lines.append(f"Wolff: n_measure={w['n_measure']}, n_burn={w['n_burn']}, "
                 f"n_seeds={w['n_seeds']}, master_seed={w['master_seed']}")
    lines.append(f"Stream-Vertrag: {w['stream_contract']} (Audit O2: "
                 "T-unabhaengige Streams)")
    lines.append("")
    lines.append("--- Messung ---")
    for L, row in report["measurements"].items():
        cells = []
        for T, m in row.items():
            xi = m["xi_ratio_mean"]
            xi_txt = "None" if xi is None else f"{xi:.4f}"
            xis = m["xi_ratio_sem"]
            xis_txt = "-" if xis is None else f"{xis:.4f}"
            cells.append(f"  T={T}: U4={m['u4_mean']:.4f}+-{m['u4_sem']:.4f}"
                         f"  xi/L={xi_txt}+-{xis_txt}")
        lines.append(f"  L={L}:")
        lines.extend(cells)
    lines.append("")
    lines.append("--- Paar-Analyse (D = O(L2) - O(L1)) ---")
    for obs in ("xi_ratio", "u4"):
        lines.append(f"  Observable {obs}:")
        for key, res in report["pair_analysis"][obs].items():
            ts = res["T_splay"]
            ts_txt = "None" if ts is None else f"{ts:.2f}"
            ci = res["splay_ci"]
            lines.append(f"    {key}: T_splay={ts_txt}  "
                         f"Bootstrap-CI=[{ci['ci_low']:.3f},"
                         f"{ci['ci_high']:.3f}]  "
                         f"none_frac={ci['none_frac']:.2f}")
            if res["crossings"]:
                for c in res["crossings"]:
                    lines.append(f"      Crossing T={c['T_cross']:.4f} "
                                 f"(Flanken-Signifikanz "
                                 f"{c['significance']:.1f} sigma)")
            else:
                lines.append("      keine adjacent-Crossings im Fenster")
    lines.append("")
    lines.append("--- PASS-Gates (Pipeline-Integritaet, kein Physik-Gate) ---")
    for k, v in report["pass_gates"].items():
        lines.append(f"  [{'PASS' if v else 'FAIL'}] {k}")
    lines.append(f"  OVERALL: {'PASS' if report['overall_pass'] else 'FAIL'}")
    lines.append(f"  Laufzeit: {report['runtime_s']:.1f}s")
    lines.append("")
    lines.append("--- FINDING / Interpretation (Vertrag: Spec §6) ---")
    largest = f"({report['lattices_L'][-2]},{report['lattices_L'][-1]})"
    ts = report["pair_analysis"]["xi_ratio"][largest]["T_splay"]
    if ts is None:
        lines.append("  NR-PHY043-01: kein persistenter xi_2/L-Splay im "
                     "Fenster — bei diesem")
        lines.append("  Budget nicht diskriminierend; O1 bleibt offen "
                     "(naechster Pfad:")
        lines.append("  Konventions-Nachweis je Referenz in SOURCES.md).")
    else:
        lines.append(f"  xi_2/L-Splay {largest} beginnt bei T={ts:.2f} "
                     f"(CI siehe oben).")
        lines.append(f"  Einordnung: intern per-Site {INTERNAL_PER_SITE} | "
                     f"Referenz {REF_TRIANGULAR}.")
        lines.append("  Splay-Temperaturen liegen konstruktionsbedingt "
                     "OBERHALB von T_BKT")
        lines.append("  (obere Schranke der Merge-Region); ein Splay-Beginn "
                     "deutlich unter")
        lines.append("  1.40 waere ein Signal GEGEN die Referenz-Lage "
                     "gewesen.")
    lines.append("")
    lines.append("Grenzen (ehrlich):")
    lines.append("  - Kein T_BKT-Bestwert; Crossing-/Splay-Lagen driften "
                 "logarithmisch.")
    lines.append("  - n_seeds=4, L<=19: keine 1%-Diskriminierung erwartbar "
                 "(Spec §6).")
    lines.append("  - Universeller (xi_2/L)*-Anker bewusst NICHT verwendet "
                 "(Rhombus-Torus,")
    lines.append("    tau=exp(i pi/3) — Hasenbusch-Wert gilt fuer tau=i).")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    report = run_phy043()
    out = (_SRC.parent / "results" /
           "260808 PHY043 triangular convention-free crossing report.txt")
    out.parent.mkdir(exist_ok=True)
    write_report(report, out)
    print(f"\nReport geschrieben: {out}")
    print("\n--- JSON-Report ---")
    print(json.dumps(_clean(report), indent=2, allow_nan=False))
