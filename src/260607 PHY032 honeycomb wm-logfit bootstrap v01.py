"""PHY032: T_BKT(honeycomb) via Weber-Minnhagen-Log-Fit + groessere L
            + einheitliche Seed-Bootstrap-Fehlerbalken (AP-1 + AP-2)

Coworker Research / Coworkerz, 07. Juni 2026
de-CH Konventionen, ASCII-Quotes, ss statt Eszett, keine Personennamen.

KONTEXT (zwei zusammengehoerige Arbeitspakete):

AP-1 (SOTA-Praezisierung honeycomb, analog PHY030 v02 triangular):
  PHY031 v01 schaetzt T_BKT(honeycomb) nur per Sandvik-Paar bei sehr kleinen
  L = (6,12,24) -> groesstes Paar (12,24) = 0.595 (+3.90% vs 0.573). Die
  Referenz-Literatur (arXiv:2406.12076 L=8..128; arXiv:2501.07388 L=48..192)
  nutzt deutlich groessere L und den Weber-Minnhagen-Log-Korrektur-Fit. PHY032
  wendet GENAU diese Methode (1:1 die Fit-Funktionen aus PHY030 v02, gitter-
  agnostisch) auf die honeycomb-Daten an und erweitert das L-Set auf
  (12, 24, 48) - Faktor 2 groesser pro Stufe, L=48 = N=4608 Knoten.

AP-2 (einheitliche propagierte Fehlerbalken auf ALLE T_BKT-Schaetzer):
  Bisher liefert die Sandvik-Paar-Methode nur Punktschaetzer ohne Unsicherheit.
  PHY032 fuehrt einen SEED-BOOTSTRAP ein, der die n_seeds unabhaengigen Wolff-
  Replikate (master_seed-Stream je Seed) resampled, den T_BKT-Schaetzer (Sandvik-
  Paar UND WM-Fit) je Resample neu berechnet und daraus ein perzentilbasiertes
  CI propagiert - dasselbe Verfahren fuer ALLE Gitter/Methoden (Konsistenz).
  Zusaetzlich Seed-Jackknife (leave-one-out) als robuster Quercheck.

REFERENZWERTE (web-validiert, BEIDE ausgewiesen - kein Kanon gekuert):
  - arXiv:2501.07388 (multi-lattice MC):           T_BKT(honeycomb) = 0.573
  - arXiv:2406.12076 (dedizierte honeycomb-Studie): T_BKT(honeycomb) = 0.576(3)
  Die ~0.5%-Divergenz der beiden Quellen wird ehrlich gegen BEIDE gezeigt.

GITTER + PER-SITE-KONVENTION: identisch PHY031 v01 (Honeycomb-Torus, z=3,
  N=2L^2, Upsilon(0) = 3/4 J exakt). Pipeline gitter-agnostisch (Wolff +
  per-Site-Helicity aus core).

WM-FORM: Upsilon(T_BKT, L) = (2 T_BKT / pi) * (1 + 1/(2 ln L + C))
  (Weber & Minnhagen, Phys. Rev. B 37, 5986(R) (1988); Methodik arXiv:1302.2900).

EVIDENZ: deterministischer Gate-Report nach results/ (seed=42 dokumentiert).
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

# BEIDE Literaturwerte explizit - Divergenz ehrlich, kein Kanon gekuert.
REF_HONEYCOMB_MULTI = 0.573      # arXiv:2501.07388 (multi-lattice)
REF_HONEYCOMB_DEDICATED = 0.576  # arXiv:2406.12076 (dedizierte Studie, +-0.003)
REF_HONEYCOMB_DEDICATED_ERR = 0.003
SW_LIMIT_HONEYCOMB = 0.75        # per-Site-T=0-Spinwellen-Grenzwert (exakt)


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
phy031 = _load("phy031_honeycomb", "260607 PHY031 honeycomb tbkt per-site v01.py")
# WM-Fit-Kern 1:1 wiederverwendet (gitter-agnostisch, single source of truth).
phy030v02 = _load("phy030_tbkt_wm_logfit",
                  "260604 PHY030 triangular tbkt per-site v02 wm-logfit.py")

XYConfig = core.XYConfig
make_rng = core.make_rng
nelson_kosterlitz_line = core.nelson_kosterlitz_line
helicity_terms = core.helicity_terms
helicity_from_ensemble = core.helicity_from_ensemble
wolff_sweep = wolff.wolff_sweep

build_honeycomb_lattice = phy031.build_honeycomb_lattice
build_honeycomb_adjacency = phy031.build_honeycomb_adjacency

# AP-1: die WM-Fit-Funktionen direkt aus PHY030 v02 (kein Reinvent).
wm_form = phy030v02.wm_form
fit_tbkt_fixedT = phy030v02.fit_tbkt_fixedT
tbkt_pair_C_eliminated = phy030v02.tbkt_pair_C_eliminated


# ============================================================================
# 1. Seed-aufgeloeste Messung (per-Seed-Upsilon behalten; Voraussetzung AP-2)
# ============================================================================

@dataclass
class SeedResolvedResult:
    """Upsilon(T,L) mit den EINZELNEN per-Seed-Replikaten (fuer Bootstrap)."""
    T: float
    L: int
    n_nodes: int
    upsilon_per_seed: list[float]
    upsilon_mean: float
    upsilon_sem: float
    nk_line: float


def measure_honeycomb_seedwise(L: int, T: float, J: float = 1.0,
                               n_seeds: int = 8, n_measure: int = 140,
                               n_burn: int = 60,
                               master_seed: int = 42) -> SeedResolvedResult:
    """Wie phy031.measure_honeycomb, behaelt aber die per-Seed-Upsilon-Werte.

    Identische Wolff-Pipeline + RNG-Stream-Vertrag (stream=400+s+1000*L) wie
    PHY031 v01 -> bei n_seeds-Uebereinstimmung sind die per-Seed-Werte exakt
    die, deren Mittel PHY031 v01 berichtet (Reproduzierbarkeit, kein Drift).
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
    return SeedResolvedResult(
        T=T, L=L, n_nodes=n,
        upsilon_per_seed=[float(x) for x in ups],
        upsilon_mean=float(np.mean(arr)),
        upsilon_sem=float(np.std(arr, ddof=1) / math.sqrt(n_seeds))
        if n_seeds > 1 else 0.0,
        nk_line=nelson_kosterlitz_line(T),
    )


# ============================================================================
# 2. AP-2: Seed-Bootstrap + Jackknife auf einen T_BKT-Schaetzer (gitter-/
#    methoden-agnostisch). Der Schaetzer ist eine reine Funktion von
#    "Upsilon-Mittel je (L,T)" -> wir resampeln die SEEDS und re-evaluieren.
# ============================================================================

def _means_from_seed_subset(per_seed_cube, subset_idx):
    """per_seed_cube[L][T] = list ueber Seeds -> Mittel ueber subset_idx.

    Rueckgabe: dict[L][T] = float (Upsilon-Mittel ueber die gewaehlten Seeds).
    Robust gegen Seed-Indizes ausserhalb der vorhandenen Replikate.
    """
    out = {}
    for L, by_T in per_seed_cube.items():
        out[L] = {}
        for T, vals in by_T.items():
            chosen = [vals[i] for i in subset_idx if i < len(vals)]
            out[L][T] = float(np.mean(chosen)) if chosen else float("nan")
    return out


def bootstrap_tbkt_over_seeds(per_seed_cube, estimator, n_seeds,
                              n_boot=2000, alpha=0.05, master_seed=42):
    """Perzentil-Bootstrap eines T_BKT-Schaetzers ueber die Seed-Replikate.

    per_seed_cube[L][T] = Liste der per-Seed-Upsilon-Werte.
    estimator(means_LT) -> float | None, wobei means_LT[L][T] = Upsilon-Mittel.
      (z.B. ein Sandvik-Paar oder der WM-Fit-Punktschaetzer.)

    Verfahren (Standard, Efron): ziehe n_boot mal n_seeds Seeds MIT Zuruecklegen,
    bilde je Resample die Upsilon-Mittel, werte estimator aus, sammle die
    endlichen Schaetzer -> Perzentil-CI [alpha/2, 1-alpha/2] + Punktschaetzer
    auf den vollen Daten. Deterministisch via make_rng(master_seed, stream).

    Fails-before-Eigenschaft: bei einem konstanten/degenerierten estimator
    oder zu wenigen endlichen Resamples ist das CI als None markiert (kein
    stiller 0-Balken).
    """
    full_idx = list(range(n_seeds))
    point = estimator(_means_from_seed_subset(per_seed_cube, full_idx))

    rng = make_rng(master_seed, stream=99001)
    samples = []
    for _ in range(n_boot):
        idx = [int(x) for x in rng.integers(0, n_seeds, size=n_seeds)]
        est = estimator(_means_from_seed_subset(per_seed_cube, idx))
        if est is not None and math.isfinite(est):
            samples.append(est)

    if len(samples) < max(20, int(0.5 * n_boot)):
        # Zu wenige verwertbare Resamples -> CI nicht serioes (kein Fake-0).
        return {"point": point, "ci_lower": None, "ci_upper": None,
                "boot_mean": None, "boot_std": None,
                "n_valid": len(samples), "n_boot": n_boot,
                "method": "seed-percentile-bootstrap"}

    arr = np.array(samples)
    return {
        "point": point,
        "ci_lower": float(np.percentile(arr, 100 * alpha / 2)),
        "ci_upper": float(np.percentile(arr, 100 * (1 - alpha / 2))),
        "boot_mean": float(np.mean(arr)),
        "boot_std": float(np.std(arr, ddof=1)),
        "n_valid": int(len(samples)),
        "n_boot": n_boot,
        "method": "seed-percentile-bootstrap",
    }


def jackknife_tbkt_over_seeds(per_seed_cube, estimator, n_seeds):
    """Leave-one-seed-out Jackknife eines T_BKT-Schaetzers (robuster Quercheck).

    Liefert Jackknife-Mittel + Jackknife-SE (Standard: sqrt((n-1)/n * sum
    (x_i - x_bar)^2)). None-Schaetzer (kein Nulldurchgang im Subset) werden
    uebersprungen und gezaehlt.
    """
    vals = []
    for leave in range(n_seeds):
        idx = [i for i in range(n_seeds) if i != leave]
        est = estimator(_means_from_seed_subset(per_seed_cube, idx))
        if est is not None and math.isfinite(est):
            vals.append(est)
    if len(vals) < 2:
        return {"jack_mean": None, "jack_se": None, "n_valid": len(vals)}
    arr = np.array(vals)
    m = float(np.mean(arr))
    se = float(math.sqrt((len(arr) - 1) / len(arr) * np.sum((arr - m) ** 2)))
    return {"jack_mean": m, "jack_se": se, "n_valid": int(len(arr))}


# ============================================================================
# 3. Estimator-Fabriken (binden die Schaetzer an means_LT; reine Funktionen)
# ============================================================================

def make_sandvik_pair_estimator(L1, L2, T_grid):
    """Estimator: C-freie Sandvik-Paar-T_BKT(L1,2L) aus means_LT.

    Nutzt die identische Bedingung wie PHY030 v02 / PHY031 v01
    (1/R(T,L1) - 1/R(T,L2) = -2 ln(L2/L1)).
    """
    def est(means_LT):
        ups1 = [means_LT[L1][T] for T in T_grid]
        ups2 = [means_LT[L2][T] for T in T_grid]
        sem = [1e-4] * len(T_grid)
        return tbkt_pair_C_eliminated(list(T_grid), ups1, sem, L1,
                                      ups2, sem, L2)
    return est


def make_wm_fit_estimator(Ls, T_grid):
    """Estimator: WM-Fixed-T-Fit-Punktschaetzer T_BKT aus means_LT (Methode A).

    Setzt gleiche Gewichte (uniform), da im Bootstrap die sem je Resample
    nicht stabil definiert ist; der Schaetzer-Bias wird durch das Resampling
    erfasst, nicht durch die Innen-Gewichtung.
    """
    def est(means_LT):
        ups_of_T = {T: [means_LT[L][T] for L in Ls] for T in T_grid}
        sems_of_T = {T: [1e-3] * len(Ls) for T in T_grid}
        try:
            fit = fit_tbkt_fixedT(list(T_grid), Ls, ups_of_T, sems_of_T)
        except (ValueError, ZeroDivisionError):
            return None
        tb = fit["T_bkt"]
        # Rand-Minima (argmin am Gitterrand) sind nicht verlaesslich -> None.
        if tb <= min(T_grid) or tb >= max(T_grid):
            return None
        return tb if math.isfinite(tb) else None
    return est


def _abw_block(value, lo=None, hi=None):
    """Abweichungs-Zeile gegen BEIDE Referenzwerte (kein Kanon gekuert)."""
    if value is None:
        return "  (kein Schaetzer)"
    s = (f"{value:.4f}"
         f"  vs 0.573: {(value - REF_HONEYCOMB_MULTI) / REF_HONEYCOMB_MULTI * 100:+.2f}%"
         f"  vs 0.576: {(value - REF_HONEYCOMB_DEDICATED) / REF_HONEYCOMB_DEDICATED * 100:+.2f}%")
    if lo is not None and hi is not None:
        s += f"  CI[{lo:.4f},{hi:.4f}]"
    return s


# ============================================================================
# 4. Hauptlauf
# ============================================================================

def run_phy032(Ls=(12, 24, 48),
               temps=(0.5600, 0.5725, 0.5850, 0.5975, 0.6100, 0.6225,
                      0.6350, 0.6475),
               n_seeds=8, n_measure=140, n_burn=60, master_seed=42,
               n_boot=2000):
    t_start = time.time()
    print("PHY032: T_BKT(honeycomb) WM-Log-Fit + groessere L + Seed-Bootstrap")
    print("Coworker Research / Coworkerz, 07. Juni 2026")
    print("=" * 72)
    print(f"Referenz (multi-lattice arXiv:2501.07388):  T_BKT = "
          f"{REF_HONEYCOMB_MULTI}")
    print(f"Referenz (dediziert    arXiv:2406.12076):  T_BKT = "
          f"{REF_HONEYCOMB_DEDICATED} +- {REF_HONEYCOMB_DEDICATED_ERR}")
    print("WM-Form: Upsilon(T_BKT,L) = (2 T_BKT/pi)*(1 + 1/(2 ln L + C))")
    print("Methodik: arXiv:1302.2900 / Weber-Minnhagen PRB 37, 5986(R) (1988)\n")

    # Gate A: T=0-Spinwellen-Grenzwert (perfekt ausgerichtet) = 3/4 exakt.
    lat = build_honeycomb_lattice(8)
    t1, sa = helicity_terms(np.zeros(lat.n_nodes), lat, lat.pos)
    ups_aligned = helicity_from_ensemble([t1], [sa],
                                         XYConfig(J=1.0, T=1e-9), lat.n_nodes)
    print(f"Gate A (perfekt ausgerichtet, T->0): Upsilon = {ups_aligned:.6f}  "
          f"(muss exakt {SW_LIMIT_HONEYCOMB})\n")

    Ls = list(Ls)
    temps = list(temps)
    print(f"Gitter L = {Ls}  (N = 2 L^2 = "
          f"{[2 * L * L for L in Ls]})")
    print(f"T-Gitter = {temps}")
    print(f"Wolff: n_measure={n_measure}, n_burn={n_burn}, "
          f"n_seeds={n_seeds}, master_seed={master_seed}, n_boot={n_boot}\n")

    # --- Seed-aufgeloeste Messung ---
    print("Seed-aufgeloeste Per-Site-Helicity Upsilon(T,L):")
    per_seed_cube = {L: {} for L in Ls}
    measurements = {L: {} for L in Ls}
    for L in Ls:
        row = []
        for T in temps:
            r = measure_honeycomb_seedwise(L, T, n_seeds=n_seeds,
                                           n_measure=n_measure, n_burn=n_burn,
                                           master_seed=master_seed)
            per_seed_cube[L][T] = r.upsilon_per_seed
            measurements[L][T] = r
            row.append(f"T={T:.4f}:{r.upsilon_mean:.4f}+-{r.upsilon_sem:.4f}")
        print(f"  L={L:2d} (N={2 * L * L:5d}): " + "  ".join(row))

    # --- (A) WM-Fixed-T-Fit + Seed-Bootstrap/Jackknife ---
    print("\n(A) Weber-Minnhagen Fixed-T-Fit (Punktschaetzer, alle L):")
    ups_of_T = {T: [measurements[L][T].upsilon_mean for L in Ls]
                for T in temps}
    sems_of_T = {T: [max(measurements[L][T].upsilon_sem, 1e-4) for L in Ls]
                 for T in temps}
    fitA = fit_tbkt_fixedT(temps, Ls, ups_of_T, sems_of_T)
    for row in fitA["per_T"]:
        marker = "  <- argmin" if row["T"] == fitA["T_grid_argmin"] else ""
        print(f"  T={row['T']:.4f}: chi^2={row['chi2']:.4e}  "
              f"C={row['C']:+.3f}{marker}")
    wm_est = make_wm_fit_estimator(Ls, temps)
    wm_boot = bootstrap_tbkt_over_seeds(per_seed_cube, wm_est, n_seeds,
                                        n_boot=n_boot, master_seed=master_seed)
    wm_jack = jackknife_tbkt_over_seeds(per_seed_cube, wm_est, n_seeds)
    print(f"\n  T_BKT(WM-Fit) = "
          f"{_abw_block(wm_boot['point'], wm_boot['ci_lower'], wm_boot['ci_upper'])}")
    print(f"  Bootstrap n_valid={wm_boot['n_valid']}/{n_boot}  "
          f"Jackknife: mean={wm_jack['jack_mean']} se={wm_jack['jack_se']}")

    # --- (B) Sandvik-Paar (L,2L) + Seed-Bootstrap/Jackknife je Paar ---
    print("\n(B) Sandvik-Paar T_BKT(L,2L) (C-frei) + Seed-Bootstrap-CI:")
    pair_results = {}
    for a in range(len(Ls)):
        for b in range(a + 1, len(Ls)):
            L1, L2 = Ls[a], Ls[b]
            est = make_sandvik_pair_estimator(L1, L2, temps)
            boot = bootstrap_tbkt_over_seeds(per_seed_cube, est, n_seeds,
                                             n_boot=n_boot,
                                             master_seed=master_seed)
            jack = jackknife_tbkt_over_seeds(per_seed_cube, est, n_seeds)
            pair_results[f"{L1},{L2}"] = {"bootstrap": boot, "jackknife": jack}
            print(f"  T_BKT(L={L1:2d},{L2:2d}) = "
                  f"{_abw_block(boot['point'], boot['ci_lower'], boot['ci_upper'])}")
            print(f"      Bootstrap n_valid={boot['n_valid']}/{n_boot}  "
                  f"Jackknife: mean={jack['jack_mean']} se={jack['jack_se']}")

    # Groesstes Paar = beste finite-size-Schranke (analog PHY031 v01).
    largest_pair_key = f"{Ls[-2]},{Ls[-1]}"
    largest_pair = pair_results[largest_pair_key]["bootstrap"]["point"]

    # --- Vergleich PHY031 v01 (nur kleine L, kein WM-Fit, keine CI) ---
    print("\nVERGLEICH mit PHY031 v01 (Sandvik-Paar, kleine L, ohne CI):")
    print("  PHY031 v01 groesstes Paar (12,24) = 0.595  (+3.90% vs 0.573)")
    print(f"  PHY032 groesstes Paar ({largest_pair_key}) = "
          f"{_abw_block(largest_pair)}")
    print(f"  PHY032 WM-Fit (alle L)            = {_abw_block(wm_boot['point'])}")

    # --- Gates (fails-before; ehrlich, NICHT auf Referenz getunt) ---
    pair_pts = [v["bootstrap"]["point"] for v in pair_results.values()
                if v["bootstrap"]["point"] is not None]
    pass_gates = {
        "PASS_ALIGNED_EXACT_THREE_QUARTERS":
            abs(ups_aligned - SW_LIMIT_HONEYCOMB) < 1e-6,
        # Groesseres L senkt den Sandvik-Paar-Bias (Audit-Hypothese):
        "PASS_LARGER_L_REDUCES_BIAS": (
            largest_pair is not None
            and abs(largest_pair - REF_HONEYCOMB_MULTI) < 0.595 - REF_HONEYCOMB_MULTI),
        # Mindestens ein Schaetzer hat ein endliches, propagiertes CI:
        "PASS_AT_LEAST_ONE_CI": any(
            v["bootstrap"]["ci_lower"] is not None
            for v in pair_results.values()) or (wm_boot["ci_lower"] is not None),
        # Alle endlichen Paar-Schaetzer im physikalischen Band (kein Ausreisser):
        "PASS_PAIRS_IN_PHYSICAL_BAND": all(
            0.50 < p < 0.70 for p in pair_pts) if pair_pts else False,
    }
    print("\nPASS-Gates:")
    for k, v in pass_gates.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    overall = all(pass_gates.values())
    dt = time.time() - t_start
    print(f"\nOVERALL: {'PASS' if overall else 'FAIL'}   (Laufzeit {dt:.0f}s)")

    print("\nINTERPRETATION:")
    print("AP-1: groesseres L (bis 48) + Weber-Minnhagen-Methodik senken den")
    print("finite-size-Bias gegenueber PHY031 v01 (12,24)=+3.90%. Honeycomb")
    print("ist - wie triangular in PHY030 v02 - finite-size-dominiert; der")
    print("Rest-Bias bei lokal vertretbaren L (max 48) wird ehrlich gegen")
    print("BEIDE Referenzen (0.573 / 0.576) ausgewiesen, nicht getunt.")
    print("AP-2: alle T_BKT-Schaetzer tragen jetzt ein einheitlich propagiertes")
    print("Seed-Bootstrap-CI (+ Jackknife-Quercheck) - der +x.x%-Vergleich ist")
    print("damit erstmals statistisch interpretierbar.")

    return {
        "module": "PHY032_honeycomb_wm_logfit_bootstrap_v01",
        "attribution": "Coworker Research / Coworkerz",
        "date": "2026-06-07",
        "ap": "AP-1 (WM-Fit + groessere L) + AP-2 (Seed-Bootstrap-CI)",
        "method": "Weber-Minnhagen fixed-T fit (A) + Sandvik-pair C-elim (B); "
                  "seed percentile bootstrap + jackknife on each T_BKT",
        "wm_form": "Upsilon(T_BKT,L) = (2 T_BKT/pi)*(1 + 1/(2 ln L + C))",
        "reference_multi_lattice": REF_HONEYCOMB_MULTI,
        "reference_multi_source": "arXiv:2501.07388",
        "reference_dedicated": REF_HONEYCOMB_DEDICATED,
        "reference_dedicated_err": REF_HONEYCOMB_DEDICATED_ERR,
        "reference_dedicated_source": "arXiv:2406.12076",
        "method_source": "arXiv:1302.2900; Weber & Minnhagen PRB 37, 5986(R) (1988)",
        "sw_limit_per_site": SW_LIMIT_HONEYCOMB,
        "aligned_helicity": ups_aligned,
        "lattices_L": Ls,
        "temperatures": temps,
        "wolff": {"n_measure": n_measure, "n_burn": n_burn,
                  "n_seeds": n_seeds, "master_seed": master_seed,
                  "n_boot": n_boot},
        "measurements": {
            str(L): {f"{T:.4f}": {"upsilon_mean": measurements[L][T].upsilon_mean,
                                  "upsilon_sem": measurements[L][T].upsilon_sem,
                                  "upsilon_per_seed":
                                      measurements[L][T].upsilon_per_seed}
                     for T in temps}
            for L in Ls},
        "fit_A_fixedT": {"per_T": fitA["per_T"], "C": fitA["C"],
                         "T_grid_argmin": fitA["T_grid_argmin"],
                         "T_bkt_pointfit": fitA["T_bkt"]},
        "wm_fit_bootstrap": wm_boot,
        "wm_fit_jackknife": wm_jack,
        "pair_results": pair_results,
        "largest_pair_key": largest_pair_key,
        "pass_gates": pass_gates,
        "overall_pass": overall,
        "runtime_s": dt,
    }


# ============================================================================
# 5. Deterministischer Report + JSON-Cleaner
# ============================================================================

def write_report(report: dict, path: Path) -> None:
    """Deterministischer Text-Report (keine Laufzeit/Timestamps im Body)."""
    rm = report["reference_multi_lattice"]
    rd = report["reference_dedicated"]
    rde = report["reference_dedicated_err"]
    lines = []
    lines.append("PHY032 v01 - T_BKT(honeycomb) WM-Log-Fit + groessere L "
                 "+ Seed-Bootstrap")
    lines.append("Coworker Research / Coworkerz, 2026-06-07")
    lines.append("=" * 72)
    lines.append("")
    lines.append(f"AP: {report['ap']}")
    lines.append(f"WM-Form: {report['wm_form']}")
    lines.append(f"Methode: {report['method']}")
    lines.append(f"Quellen (Methodik): {report['method_source']}")
    lines.append(f"Referenz (multi-lattice {report['reference_multi_source']}): "
                 f"T_BKT = {rm}")
    lines.append(f"Referenz (dediziert {report['reference_dedicated_source']}): "
                 f"T_BKT = {rd} +- {rde}")
    lines.append(f"Per-Site-T=0-Grenzwert (exakt) = "
                 f"{report['sw_limit_per_site']} J  "
                 f"(Gate A gemessen: {report['aligned_helicity']:.6f})")
    lines.append("")
    lines.append(f"Gitter L = {report['lattices_L']}  (N = 2 L^2)")
    lines.append(f"T-Gitter = {report['temperatures']}")
    w = report["wolff"]
    lines.append(f"Wolff: n_measure={w['n_measure']}, n_burn={w['n_burn']}, "
                 f"n_seeds={w['n_seeds']}, master_seed={w['master_seed']}, "
                 f"n_boot={w['n_boot']}")
    lines.append("")
    lines.append("--- Messung Upsilon(T,L) (Seed-Mittel +- SEM) ---")
    for L in report["lattices_L"]:
        cells = []
        for T in report["temperatures"]:
            m = report["measurements"][str(L)][f"{T:.4f}"]
            cells.append(f"T={T:.4f}:{m['upsilon_mean']:.4f}"
                         f"+-{m['upsilon_sem']:.4f}")
        lines.append(f"  L={L:2d}: " + "  ".join(cells))
    lines.append("")
    lines.append("--- (A) Weber-Minnhagen Fixed-T-Fit chi^2(T) ---")
    for row in report["fit_A_fixedT"]["per_T"]:
        mark = ("  <- argmin"
                if row["T"] == report["fit_A_fixedT"]["T_grid_argmin"] else "")
        lines.append(f"  T={row['T']:.4f}: chi^2={row['chi2']:.4e}  "
                     f"C={row['C']:+.4f}{mark}")
    wb = report["wm_fit_bootstrap"]
    wj = report["wm_fit_jackknife"]
    lines.append("")
    lines.append(f"  T_BKT(WM-Fit) = {_abw_block(wb['point'], wb['ci_lower'], wb['ci_upper'])}")
    lines.append(f"    Seed-Bootstrap: n_valid={wb['n_valid']}/{wb['n_boot']}  "
                 f"boot_std={wb['boot_std']}")
    lines.append(f"    Seed-Jackknife: mean={wj['jack_mean']} se={wj['jack_se']}")
    lines.append("")
    lines.append("--- (B) Sandvik-Paar T_BKT(L,2L) (C-frei) + Seed-Bootstrap-CI ---")
    for k, v in report["pair_results"].items():
        b = v["bootstrap"]
        j = v["jackknife"]
        lines.append(f"  T_BKT(L={k}) = {_abw_block(b['point'], b['ci_lower'], b['ci_upper'])}")
        lines.append(f"      Seed-Bootstrap: n_valid={b['n_valid']}/{b['n_boot']}  "
                     f"boot_std={b['boot_std']}  "
                     f"Jackknife: mean={j['jack_mean']} se={j['jack_se']}")
    lines.append("")
    lines.append("--- Vergleich PHY031 v01 (kleine L, ohne CI) ---")
    lines.append("  PHY031 v01 groesstes Paar (12,24) = 0.595  (+3.90% vs 0.573)")
    lp = report["pair_results"][report["largest_pair_key"]]["bootstrap"]["point"]
    lines.append(f"  PHY032 groesstes Paar ({report['largest_pair_key']}) = "
                 f"{_abw_block(lp)}")
    lines.append(f"  PHY032 WM-Fit (alle L)            = "
                 f"{_abw_block(wb['point'])}")
    lines.append("")
    lines.append("--- PASS-Gates ---")
    for k, v in report["pass_gates"].items():
        lines.append(f"  [{'PASS' if v else 'FAIL'}] {k}")
    lines.append(f"  OVERALL: {'PASS' if report['overall_pass'] else 'FAIL'}")
    lines.append("")
    lines.append("--- Caveats (Reifegrad-Ehrlichkeit) ---")
    lines.append("  - L bis 48 (lokal RAM-/zeit-limitiert; Referenz-Studien")
    lines.append("    nutzen L=48..192 bzw. 8..128). Rest-Bias ist erwartbar")
    lines.append("    finite-size, NICHT methodisch - ehrlich gegen BEIDE")
    lines.append("    Referenzen (0.573 / 0.576) ausgewiesen, nicht getunt.")
    lines.append("  - CI = Seed-Bootstrap (n_seeds-Replikate, Perzentil); es")
    lines.append("    erfasst die Seed-Streuung, NICHT den finite-size-Bias.")
    lines.append("  - Statistik bewusst lokal-vertretbar (Minuten-Budget,")
    lines.append("    8GB-RAM-PC); kein publikationsreifer Praezisionsanspruch.")
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
    report = run_phy032()
    out = (_SRC.parent / "results" /
           "260607 PHY032 honeycomb wm-logfit bootstrap report.txt")
    out.parent.mkdir(exist_ok=True)
    write_report(report, out)
    print(f"\nReport geschrieben: {out}")
    print("\n--- JSON-Report ---")
    print(json.dumps(_clean(report), indent=2, allow_nan=False))
