"""PHY052 - Mixing-/Thermalisierungs-Pilot vor der W4-v03-Produktion (PHY049).

Coworker Research / Coworkerz, 21. November 2026
de-CH Konventionen, ASCII-Quotes, ss statt Eszett, keine Personennamen.

Vertrag: `spec/261121 PHI HEX phy052 w4 v03 mixing pilot preregistration v01.md`
(vor jeder Pilot-Datennahme committet). Der Pilot ist NUR Diagnostik: er
misst, ob der PHY049-Kernel bei den vorregistrierten Sweep-Zahlen
(1000 Thermalisierung + 4000 Messung) aus geordnetem UND zufaelligem Start
dieselbe stationaere Verteilung erreicht. Er liefert keine Physikzahl und
kann von PHY050 nicht adjudiziert werden (eigene Seed-Namensraeume).

Herkunft: Design (3x3-Zellen, 4 Ketten, kalt/heiss, Seeds 79M/89M,
Einflussreihe des Ratio-of-Means) aus dem Cross-Family-Bundle
`PhiHex_W4_v03_Audit_Hardening_v3_2026-10-04.zip` (007acc2, 05_REVIEW_BUNDLE).
Die Diagnostik ist NICHT aus dem Bundle uebernommen, sondern eine
zeilengenaue Portierung der Stan-Referenz `stan-dev/posterior`
(R/convergence.R: z_scale, fold_draws, .split_chains, .rhat, .ess, rhat,
ess_bulk). Dabei korrigiert gegenueber dem Bundle: Normal-Score-Nenner
S + 1/4 (posterior: S - 2c + 1, c = 3/8; das Bundle hatte S - 1/4), biased
Autokovarianz 1/N (Geyer 1992, posterior), Multi-Ketten-ESS statt
Ketten-Einzel-ESS, Geyer "improved estimate" und ESS-Kappung wie posterior.

Literatur: Vehtari, Gelman, Simpson, Carpenter, Buerkner, Bayesian Analysis
16, 667 (2021), doi:10.1214/20-BA1221; Geyer, Stat. Sci. 7, 473 (1992).
"""
from __future__ import annotations

import os

# ERZWUNGEN vor dem numpy-Import (AGENTS.md): reproduzierbare Reduktionen.
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import hashlib  # noqa: E402
import importlib.util  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import platform  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.stats import norm, rankdata  # noqa: E402

_SRC = Path(__file__).resolve().parent
_ROOT = _SRC.parent


def _load(name: str, filename: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _SRC / filename)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


_p49 = _load("phy049_corr_ratio",
             "260927 PHY049 honeycomb correlation ratio v01.py")
_p45 = _p49._p45

# ---- vorregistrierte Konstanten (Spec v01; Test-gebunden) -----------------
PILOT_LS = (48, 96, 192)
PILOT_TS = (0.560, 0.575, 0.595)
PILOT_N_CHAINS = 4                       # 2 geordnet + 2 zufaellig
PILOT_N_THERM = 1000                     # = W4V3_N_THERM
PILOT_N_MEAS = 4000                      # = W4V3_N_MEAS
PILOT_SEED_BASE = 79_000_000
PILOT_INIT_SEED_BASE = 89_000_000
MAX_RHAT = 1.01
MIN_ESS_BULK = 400.0
MAX_COLD_HOT_Z = 3.5                     # |R_kalt - R_heiss| / SE
GATE_OBSERVABLE = "ratio_influence"
SPEC = "spec/261121 PHI HEX phy052 w4 v03 mixing pilot preregistration v01.md"
REPORT_STEM = "261121 PHY052 honeycomb correlation-ratio mixing pilot report"
_RUNTIME_SOURCE_PATHS = (
    "src/260607 PHY031 honeycomb tbkt per-site v01.py",
    "src/260926 PHY045 helicity normalization O1 test v01.py",
    "src/260927 PHY049 honeycomb correlation ratio v01.py",
    "src/261121 PHY052 honeycomb correlation-ratio mixing pilot v01.py",
    SPEC,
)


# ============================================================================
# 1. Diagnostik: Portierung von stan-dev/posterior (R/convergence.R)
#    Konvention hier: x hat Form (n_chains, n_draws).
# ============================================================================

def _chains(x) -> np.ndarray | None:
    """(n_chains, n_draws)-Array oder None (fail-closed, wie posterior
    should_return_NA: NA/Inf oder konstante Eingabe)."""
    a = np.asarray(x, dtype=float)
    if a.ndim != 2 or a.shape[0] < 1 or a.shape[1] < 4:
        return None
    if not np.all(np.isfinite(a)):
        return None
    if float(np.max(a) - np.min(a)) <= np.finfo(float).eps:
        return None
    return a


def split_chains(x: np.ndarray) -> np.ndarray:
    """posterior .split_chains: erste floor(n/2) und letzte Haelfte je
    Kette; bei ungeradem n faellt die mittlere Ziehung weg."""
    n = x.shape[1]
    half = n / 2.0
    first = x[:, :math.floor(half)]
    second = x[:, math.ceil(half + 1) - 1:]
    return np.vstack([first, second])


def z_scale(x: np.ndarray) -> np.ndarray:
    """posterior z_scale: Raenge ueber ALLE Ziehungen (ties average),
    Phi^-1((r - c) / (S - 2c + 1)), c = 3/8 (Nenner S + 1/4)."""
    c = 3.0 / 8.0
    r = rankdata(x.ravel(), method="average")
    s = r.size
    return norm.ppf((r - c) / (s - 2.0 * c + 1.0)).reshape(x.shape)


def fold_draws(x: np.ndarray) -> np.ndarray:
    """posterior fold_draws: |x - median(x)| ueber alle Ziehungen."""
    return np.abs(x - np.median(x))


def _rhat_basic(x: np.ndarray) -> float:
    """posterior .rhat auf (n_chains, n_draws)."""
    n = x.shape[1]
    chain_mean = x.mean(axis=1)
    chain_var = x.var(axis=1, ddof=1)
    var_between = n * float(np.var(chain_mean, ddof=1))
    var_within = float(np.mean(chain_var))
    if not var_within > 0.0:
        return math.nan
    return math.sqrt((var_between / var_within + n - 1.0) / n)


def autocovariance(x: np.ndarray) -> np.ndarray:
    """posterior autocovariance: FFT, "biased" (1/N) nach Geyer (1992)."""
    x = np.asarray(x, dtype=float)
    n = x.size
    varx = float(np.var(x, ddof=1))
    if varx == 0.0:
        return np.zeros(n)
    m = 1 << (n - 1).bit_length()           # >= n (posterior: nextn)
    yc = np.zeros(2 * m)
    yc[:n] = x - x.mean()
    f = np.fft.fft(yc)
    ac = np.fft.ifft(np.abs(f) ** 2).real[:n]
    return ac / ac[0] * varx * (n - 1) / n


def _ess(x: np.ndarray) -> float:
    """posterior .ess: Multi-Ketten-ESS mit Geyers initial-positive und
    initial-monotone Folge, "improved estimate" und Kappung."""
    n_chains, n = x.shape
    if n < 3:
        return math.nan
    acov = np.array([autocovariance(c) for c in x])          # (chains, n)
    acov_means = acov.mean(axis=0)
    mean_var = acov_means[0] * n / (n - 1)
    var_plus = mean_var * (n - 1) / n
    if n_chains > 1:
        var_plus += float(np.var(x.mean(axis=1), ddof=1))
    rho = np.zeros(n)
    t = 0
    rho_even = 1.0
    rho[0] = rho_even
    rho_odd = 1.0 - (mean_var - acov_means[1]) / var_plus
    rho[1] = rho_odd
    while (t < n - 5 and not math.isnan(rho_even + rho_odd)
           and rho_even + rho_odd > 0):
        t += 2
        rho_even = 1.0 - (mean_var - acov_means[t]) / var_plus
        rho_odd = 1.0 - (mean_var - acov_means[t + 1]) / var_plus
        if rho_even + rho_odd >= 0:
            rho[t] = rho_even
            rho[t + 1] = rho_odd
    max_t = t
    if rho_even > 0:
        rho[max_t] = rho_even
    t = 0
    while t <= max_t - 4:
        t += 2
        if rho[t] + rho[t + 1] > rho[t - 2] + rho[t - 1]:
            rho[t] = (rho[t - 2] + rho[t - 1]) / 2.0
            rho[t + 1] = rho[t]
    ess = float(n_chains * n)
    tau_hat = -1.0 + 2.0 * float(np.sum(rho[:max_t])) + float(rho[max_t])
    tau_hat = max(tau_hat, 1.0 / math.log10(ess))
    return ess / tau_hat


def rhat(x) -> float:
    """posterior rhat: max(rank-normalisiertes split-Rhat, gefaltetes);
    NaN (= Gate rot) fuer ungueltige oder konstante Eingabe."""
    a = _chains(x)
    if a is None or a.shape[0] < 2:
        return math.nan
    bulk = _rhat_basic(z_scale(split_chains(a)))
    tail = _rhat_basic(z_scale(split_chains(fold_draws(a))))
    return max(bulk, tail) if all(map(math.isfinite, (bulk, tail))) \
        else math.nan


def ess_bulk(x) -> float:
    """posterior ess_bulk: ESS der rank-normalisierten Split-Ketten."""
    a = _chains(x)
    if a is None:
        return math.nan
    return _ess(z_scale(split_chains(a)))


def ratio_influence(q: np.ndarray, h: np.ndarray) -> tuple[np.ndarray, float]:
    """Einflussreihe des Ratio-of-Means R = <h>/<q> (gepoolt ueber alle
    Ketten): u = (h - R q) / <q>. Var(R_hat) ~ Var(mean u); das ist die
    entscheidungsrelevante Groesse von PHY049, nicht das Mittel der
    Einzel-Quotienten."""
    qbar = float(np.mean(q))
    if not (math.isfinite(qbar) and abs(qbar) > _p49.W4V3_MIN_DEN):
        raise ValueError("ratio denominator mean is zero/non-finite")
    ratio = float(np.mean(h)) / qbar
    return (h - ratio * q) / qbar, ratio


# ============================================================================
# 2. Pilot-Design und Ketten
# ============================================================================

def pilot_seed_for(L: int, t_idx: int, s: int) -> int:
    return PILOT_SEED_BASE + 1000 * L + 100 * t_idx + s


def pilot_init_seed_for(L: int, t_idx: int, s: int) -> int:
    return PILOT_INIT_SEED_BASE + 1000 * L + 100 * t_idx + s


def init_mode_for(s: int) -> str:
    """Kette 0, 2: geordnet (kalt); 1, 3: gleichverteilt zufaellig (heiss)."""
    return "ordered" if s % 2 == 0 else "random"


def t_index(T: float) -> int:
    """Index von T im W4-v03-Gitter (Seed-Koordinate wie PHY049)."""
    grid = [round(x, 4) for x in _p49.W4V3_T_GRID]
    return grid.index(round(float(T), 4))


def initial_theta(n: int, mode: str, init_seed: int) -> np.ndarray:
    if mode == "ordered":
        return np.zeros(n)
    if mode == "random":
        rng = np.random.default_rng(init_seed)
        return rng.uniform(-math.pi, math.pi, size=n)
    raise ValueError(mode)


def _chain_job(args: tuple) -> dict:
    """Eine Kette: (L, T, s, n_therm, n_meas) -> Messreihen g(L/4), g(L/2).
    Der Kernel kommt aus PHY049 und ist VAL-BIT-gebunden."""
    L, T, s, n_therm, n_meas = args
    _p49._validate_L(L)
    lat = _p45.build("honeycomb", L)
    if lat.n != 2 * L * L:
        raise RuntimeError("unexpected honeycomb indexing")
    nbr, deg = _p45._nbr_arrays(lat)
    k = t_index(T)
    seed = pilot_seed_for(L, k, s)
    init_seed = pilot_init_seed_for(L, k, s)
    mode = init_mode_for(s)
    th = initial_theta(lat.n, mode, init_seed)
    run = _p49._production_run_backend()
    t0 = time.perf_counter()
    data = np.asarray(run(th, nbr, deg, 1.0 / T, n_therm, n_meas, seed, L),
                      dtype=float)
    if data.shape != (n_meas, 2) or not np.all(np.isfinite(data)):
        raise RuntimeError("kernel returned malformed pilot data")
    return {"L": L, "T": T, "s": s, "seed": seed, "init_seed": init_seed,
            "init_mode": mode, "q": data[:, 0], "h": data[:, 1],
            "wall_s": time.perf_counter() - t0}


def cell_diagnostics(q: np.ndarray, h: np.ndarray, modes: list) -> dict:
    """Gate-Diagnostik einer (L, T)-Zelle aus (n_chains, n_meas)-Reihen."""
    u, ratio = ratio_influence(q, h)
    r_u, e_u = rhat(u), ess_bulk(u)
    cold = [i for i, m in enumerate(modes) if m == "ordered"]
    hot = [i for i, m in enumerate(modes) if m == "random"]

    def group_ratio(idx):
        if not idx:                      # Design verletzt -> z NaN -> rot
            return math.nan
        return float(np.mean(h[idx])) / float(np.mean(q[idx]))

    def group_se(idx):
        if not idx:
            return math.nan
        ug = (h[idx] - ratio * q[idx]) / float(np.mean(q))
        ess = ess_bulk(ug)
        return (float(np.std(ug, ddof=1)) / math.sqrt(ess)
                if math.isfinite(ess) and ess > 0 else math.nan)

    d_r = group_ratio(cold) - group_ratio(hot)
    se = math.hypot(group_se(cold), group_se(hot))
    z = d_r / se if (math.isfinite(se) and se > 0) else math.nan
    reasons = []
    if not (math.isfinite(r_u) and r_u <= MAX_RHAT):
        reasons.append(f"rhat({GATE_OBSERVABLE}) = {r_u} > {MAX_RHAT}")
    if not (math.isfinite(e_u) and e_u >= MIN_ESS_BULK):
        reasons.append(f"ess_bulk({GATE_OBSERVABLE}) = {e_u} < {MIN_ESS_BULK}")
    # Zweistart-Test (Spec Abschnitt 4): Rhat erkennt einen Versatz erst ab
    # ~0.3 sigma je Ziehung; fuer Bias gegen den Standardfehler des
    # Mittelwerts ist der Vergleich kalt/heiss in Einheiten des MCSE noetig.
    if not (math.isfinite(z) and abs(z) <= MAX_COLD_HOT_Z):
        reasons.append(f"|cold_minus_hot_z| = {z} > {MAX_COLD_HOT_Z}")
    return {
        "ratio_of_means": ratio,
        "gate_observable": GATE_OBSERVABLE,
        "rhat": r_u, "ess_bulk": e_u,
        "ess_fraction": e_u / u.size if math.isfinite(e_u) else math.nan,
        "passed": not reasons, "reasons": reasons,
        "cold_minus_hot_ratio": d_r, "cold_minus_hot_se": se,
        "cold_minus_hot_z": z,
        # Diagnostik, NICHT entscheidend (Spec Abschnitt 4):
        "rhat_g_quarter": rhat(q), "ess_bulk_g_quarter": ess_bulk(q),
        "rhat_g_half": rhat(h), "ess_bulk_g_half": ess_bulk(h),
    }


def run_pilot(Ls=PILOT_LS, Ts=PILOT_TS, n_chains=PILOT_N_CHAINS,
              n_therm=PILOT_N_THERM, n_meas=PILOT_N_MEAS,
              max_workers: int = 4) -> dict:
    jobs = [(L, T, s, n_therm, n_meas) for L in sorted(Ls, reverse=True)
            for T in Ts for s in range(n_chains)]
    backend = ("numba" if _p49._production_run_backend()
               is not _p49._py_run_corr else "python")
    t0 = time.perf_counter()
    res = _p45.map_jobs(_chain_job, jobs, max_workers)
    wall = time.perf_counter() - t0
    cells = []
    for L in Ls:
        for T in Ts:
            ch = sorted((r for r in res if r["L"] == L and r["T"] == T),
                        key=lambda r: r["s"])
            q = np.array([r["q"] for r in ch])
            h = np.array([r["h"] for r in ch])
            modes = [r["init_mode"] for r in ch]
            cells.append({
                "L": L, "T": T, "t_idx": t_index(T),
                "seeds": [r["seed"] for r in ch],
                "init_seeds": [r["init_seed"] for r in ch],
                "init_modes": modes,
                "chain_mean_g_quarter": [float(x) for x in q.mean(axis=1)],
                "chain_mean_g_half": [float(x) for x in h.mean(axis=1)],
                "wall_s_by_chain": [r["wall_s"] for r in ch],
                **cell_diagnostics(q, h, modes)})
    return {"Ls": list(Ls), "Ts": list(Ts), "n_chains": n_chains,
            "n_therm": n_therm, "n_meas": n_meas, "max_workers": max_workers,
            "backend": backend, "wall_s": wall,
            "cpu_s": float(sum(r["wall_s"] for r in res)), "cells": cells}


# ============================================================================
# 3. Entscheid, Freigabe-Gate, Report
# ============================================================================

def _file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def runtime_provenance() -> dict:
    return {
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "numba_version": (str(_p45.numba.__version__)
                          if _p45.HAVE_NUMBA else None),
        "sys_platform": sys.platform, "machine": platform.machine(),
        "thread_env": {v: os.environ.get(v) for v in (
            "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")},
        "source_sha256": {p: _file_sha256(_ROOT / p)
                          for p in _RUNTIME_SOURCE_PATHS},
    }


def is_preregistered_design(pilot: dict) -> bool:
    """Nur der vorregistrierte Plan kann PASS ergeben (keine verkleinerten
    Probelaeufe, keine Zusatz-/Fehlzellen)."""
    return (pilot.get("Ls") == list(PILOT_LS)
            and pilot.get("Ts") == list(PILOT_TS)
            and pilot.get("n_chains") == PILOT_N_CHAINS
            and pilot.get("n_therm") == PILOT_N_THERM
            and pilot.get("n_meas") == PILOT_N_MEAS
            and len(pilot.get("cells", [])) == len(PILOT_LS) * len(PILOT_TS)
            and {(c["L"], c["T"]) for c in pilot["cells"]}
            == {(L, T) for L in PILOT_LS for T in PILOT_TS})


def decision(pilot: dict) -> dict:
    design = is_preregistered_design(pilot)
    failed = [{"L": c["L"], "T": c["T"], "reasons": c["reasons"]}
              for c in pilot.get("cells", []) if c.get("passed") is not True]
    ok = design and not failed
    return {"preregistered_design": design, "failed_cells": failed,
            "label": ("MIXING_PILOT_PASS" if ok
                      else "MIXING_PILOT_FAIL")}


def production_release(pilot_decision: dict, phy049_preflight: dict) -> str:
    """Freigabe-Gate (Audit 2026-10-04, release gate v01): PHY049 darf nur
    starten, wenn der Pilot PASS ist UND der native PHY049-Preflight exakt
    die fuenf Gates literal True traegt. Kein Physik-Urteil."""
    gates = phy049_preflight.get("gates")
    native = (isinstance(gates, dict)
              and set(gates) == set(_p49._p50.W4V3_REQUIRED_PREFLIGHT_GATES)
              and all(v is True for v in gates.values())
              and phy049_preflight.get("production_measurement_eligible")
              is True)
    if pilot_decision.get("label") == "MIXING_PILOT_PASS" and native:
        return "ELIGIBLE_TO_START_PHY049_MEASUREMENT"
    return "DO_NOT_START_PHY049_MEASUREMENT"


def build_report(pilot: dict) -> dict:
    dec = decision(pilot)
    return {"module": "PHY052_honeycomb_correlation_ratio_mixing_pilot_v01",
            "attribution": "Coworker Research / Coworkerz",
            "spec": SPEC,
            "claim_ceiling": ("PREPRODUCTION_MIXING_DIAGNOSTIC_ONLY: no "
                              "T_BKT, not adjudicable by PHY050"),
            "physics_interpretation_enabled": False,
            "thresholds": {"max_rhat": MAX_RHAT,
                           "min_ess_bulk": MIN_ESS_BULK,
                           "max_cold_hot_z": MAX_COLD_HOT_Z,
                           "gate_observable": GATE_OBSERVABLE},
            "seed_namespaces": {"pilot": PILOT_SEED_BASE,
                                "init": PILOT_INIT_SEED_BASE,
                                "production": _p49.W4V3_SEED_BASE},
            "runtime_provenance": runtime_provenance(),
            "pilot": pilot, "decision": dec}


def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(x) for x in o]
    if isinstance(o, (np.floating, float)):
        v = float(o)
        return v if math.isfinite(v) else None
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def write_text_report(rep: dict, path: Path) -> None:
    p, d = rep["pilot"], rep["decision"]

    def f(x, fmt=".4f"):
        return "null" if x is None else format(x, fmt)
    lines = ["PHY052 - Mixing-Pilot W4 v03 (PHY049-Kernel, kalt/heiss)",
             "Coworker Research / Coworkerz, 2026-11-21", "=" * 70,
             f"Spec: {rep['spec']}",
             f"L {p['Ls']}, T {p['Ts']}, {p['n_chains']} Ketten je Zelle "
             f"(geordnet/zufaellig), {p['n_therm']}+{p['n_meas']} Sweeps, "
             f"Backend {p['backend']}",
             f"Wall {p['wall_s']:.0f} s auf {p['max_workers']} Prozessen, "
             f"CPU {p['cpu_s']:.0f} s",
             f"Gate: rhat({rep['thresholds']['gate_observable']}) <= "
             f"{rep['thresholds']['max_rhat']}, ess_bulk >= "
             f"{rep['thresholds']['min_ess_bulk']:.0f}, |z kalt-heiss| <= "
             f"{rep['thresholds']['max_cold_hot_z']}", "",
             "   L      T     R      rhat_u  ess_u   rhat_q  rhat_h  "
             "kalt-heiss z  Gate"]
    for c in p["cells"]:
        lines.append(
            f"  {c['L']:3d}  {c['T']:.3f}  {f(c['ratio_of_means'])}  "
            f"{f(c['rhat'])}  {f(c['ess_bulk'], '6.0f')}  "
            f"{f(c['rhat_g_quarter'])}  {f(c['rhat_g_half'])}  "
            f"{f(c['cold_minus_hot_z'], '+6.2f')}        "
            f"{'PASS' if c['passed'] else 'FAIL'}")
    lines += ["", f"Entscheid: {d['label']} (Design vorregistriert: "
              f"{d['preregistered_design']})"]
    for fc in d["failed_cells"]:
        lines.append(f"  FAIL L={fc['L']} T={fc['T']}: {fc['reasons']}")
    lines += ["", "Grenzen: Diagnostik, keine Physikzahl; ein PASS belegt",
              "Stationaritaet/Durchmischung im Pilot, nicht die Richtigkeit",
              "von T_BKT. Die Pilot-Seeds sind fuer PHY050 nicht verwertbar."]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv=None) -> int:
    rep = _clean(build_report(run_pilot()))
    res = _ROOT / "results"
    (res / f"{REPORT_STEM}.json").write_text(
        json.dumps(rep, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    write_text_report(rep, res / f"{REPORT_STEM}.txt")
    print(f"PHY052: {rep['decision']['label']}; Report results/"
          f"{REPORT_STEM}.{{json,txt}}")
    return 0 if rep["decision"]["label"] == "MIXING_PILOT_PASS" else 3


if __name__ == "__main__":
    raise SystemExit(main())
