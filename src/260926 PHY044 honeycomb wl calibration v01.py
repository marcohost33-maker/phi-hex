"""PHY044 - W4-Kalibrierung honeycomb Wang-Landau: Kosten, Sampler-Praezision,
Budget-Hebel - BLIND (kein T_BKT-Lagewert) + bit-identischer Numba-Kernel.

Coworker Research / Coworkerz, 26. September 2026
de-CH Konventionen, ASCII-Quotes, ss statt Eszett, keine Personennamen.

================================ ZWECK =================================
Issue #45 ("Before any large run"): vor jedem W4-Produktionslauf (L>=64)
muessen die realen Kosten und die WL-Konvergenz GEMESSEN sein - bisher gab
es nur eine Extrapolation (L^3.0..3.4 as-run, ~L^5.4 konvergiert, UNBELEGT).
PHY044 misst:

  K1 Kosten: wl_sweeps (Zeit bis zur 1/t-Phase), Produktion, Spin-Updates
     und Wall-Zeit je Walker fuer L in einer Leiter; Exponenten-Fit.
  K2 Praezision: Walker-Spread-Validitaets-Domaene (PHY042-Semantik, >= 3
     Walker an JEDEM L) beim PHY042-Budget-Rezept; Walker-Streuung der
     C-eliminierten Paar-Crossings (NUR die Streuung, nie die Lage).
  K3 Hebel (L=48): welcher Budget-Knopf erweitert die Domaene - 4x
     Produktion, lnf_final 1e-6, und eine Zerlegung des Spreads in
     g(E)- und Produktions-Anteil (feste g(E), unabhaengige Produktionen).

Vertrag: `spec/260926 PHI HEX w4 honeycomb preregistration v01.md` (§5
Kalibrier-Protokoll, VOR den Laeufen festgelegt).

BLIND-VERTRAG: PHY044 berechnet Paar-Crossings nur intern, um deren
Walker-Streuung zu messen. Lage-Werte (T_BKT) werden weder ausgegeben noch
gespeichert - die Kalibrierung legt Budgets fest, nicht Physik.

================================ KERNEL ================================
`wl_entropic_fast` ist ein bit-identischer Zwilling von PHY041
`wl_entropic_lattice`: identischer RNG-Verbrauch (dieselben numpy-Block-
Ziehungen je Sweep, in derselben Reihenfolge), identische Arithmetik je
Spin-Update (math.cos/log/exp, Float-Floor-Division), identische
Resync-/Phasenlogik; nur die innere Update-Schleife ist Numba-kompiliert.
Die mikrokanonischen Aggregate laufen weiterhin ueber dieselben numpy-
Funktionen (np.dot-Summationsreihenfolge unveraendert).
Gates (tests/test_phy044_wl_calibration.py + Report hier):
  - VAL-BIT-S: Python-Fallback UND Numba == Original auf kleinen Gittern
    (lng, Maske, alle Aggregate, wl_sweeps exakt gleich);
  - VAL-BIT-P: Numba reproduziert ALLE 7 committeten PHY042-Walker
    (L=24/32/48) - Trajektorie exakt (wl_sweeps, Bin-Belegung); die auf
    einer anderen Plattform committeten Y2/Y4-Kurven bis rtol 1e-12
    (ULP-Ebene der numpy-Reduktionen np.dot/np.sum, nicht des Samplers);
  - Stufe speed: Original-Python vs Numba auf DERSELBEN Maschine,
    identischer L=24-Produktionsjob, bit-genau (results_identical).
Ohne numba laeuft derselbe Code als reines Python (numpy-Arrays) - korrekt,
aber langsam; `backend="numba"` verlangt numba explizit.

Provenance: Wang & Landau PRL 86, 2050 (2001); Belardinelli & Pereyra JCP
127, 184105 (2007); Lam, Pitrou, Seibert, "Numba", LLVM-HPC 2015.
"""
from __future__ import annotations

import os

# BLAS-Threads VOR dem numpy-Import pinnen (Befund 2026-09-26, erster
# Leiter-Lauf): jeder Prozess-Pool-Worker startete einen eigenen OpenBLAS-
# Pool (7 Threads/Prozess, Last 13 auf 4 Kernen) - das verfaelscht die
# K1-Wall-Zeiten, und OpenBLAS teilt ddot auf mehrere Threads auf: die
# Reduktionsreihenfolge der Aggregate haengt dann von der Thread-Zahl ab
# (verifiziert 2026-09-26, OpenBLAS 0.3.34: np.dot bitweise verschieden bei
# 4 vs 1 Thread ab n=12288 = honeycomb-Bindungen L=64; bis n=6912 = L=48
# identisch -> committete PHY042-Laeufe L<=48 unberuehrt).
# setdefault: eine bewusste Umgebungs-Vorgabe gewinnt und wird im Report
# protokolliert. Wirkt nur, wenn numpy noch nicht geladen ist (CLI-Lauf).
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import hashlib  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import time  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402

import numpy as np  # noqa: E402

import importlib.util  # noqa: E402
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

try:  # optionaler Beschleuniger (requirements-dev); Fallback = reines Python
    import numba
    HAVE_NUMBA = True
except ImportError:  # pragma: no cover - Fallback-Pfad
    numba = None
    HAVE_NUMBA = False

_SRC = Path(__file__).resolve().parent
_ROOT = _SRC.parent


def _load(name: str, filename: str):
    """Portabler Loader (Muster PHY030 v02 / PHY031-033 / PHY040-042)."""
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _SRC / filename)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


_p41 = _load("phy041_honeycomb_wl",
             "260702 PHY041 honeycomb wang-landau entropic helicity v01.py")
_p42 = _load("phy042_honeycomb_wl_fss", "260706 PHY042 honeycomb wl fss v01.py")

make_rng = _p41.make_rng
honeycomb_arrays = _p41.honeycomb_arrays
WLResult = _p41.WLResult
OBS = _p41.OBS
_full_energy = _p41._full_energy
_aggregates_dir = _p41._aggregates_dir
upsilon_curves = _p41.upsilon_curves
tbkt_pair_from_curves = _p41.tbkt_pair_from_curves
_STREAM_WL = _p41._STREAM_WL

# RNG-Stream-Vertrag PHY044 (kollisionsfrei zu 300+..700+ inkl. +1000L-
# Termen, 900+ PHY043 (+1000L +100000 t_idx, max ~ 1.8e6), 90000+1000w
# PHY042): K2/K3-Walker 3_000_000 + 10_000*w (+L im Kernel); Zusatz-
# Produktionen der Zerlegung 3_500_000 + 10_000*p (+L).
_STREAM_CAL_BASE = 3_000_000
_STREAM_DECOMP_BASE = 3_500_000

# Blind-Vertrag: der Kalibrier-Report darf keine Lage-Werte enthalten.
BLIND_FORBIDDEN_KEYS = ("tbkt", "pair_tbkt", "T_BKT", "t_bkt")

# W4-Vorregistrierung (spec/260926 ... w4 honeycomb preregistration v01.md);
# Test-gebunden an die Spec-Zahlen (Drift-Guard).
W4_BAND = (0.560, 0.580)
W4_MIN_WALKERS = 3
W4_MIN_PAIRS = 3
W4_SIGMA_MAX = 0.010
W4_SIGMA_FSS_FLOOR = 0.003
W4_DOMAIN_THRESHOLD = 0.04
W4_T_REQ = 0.62
W4_LADDER = (32, 48, 64, 96)


def w4_verdict(t_w4: float | None, sigma_sampler: float | None,
               sigma_fss: float | None, n_quotable_pairs: int) -> dict:
    """Entscheidungsregel §7 der W4-Vorregistrierung (Reihenfolge fix).

    sigma_FSS wird auf die Untergrenze W4_SIGMA_FSS_FLOOR angehoben (§6);
    fehlende Groessen fuehren fail-closed zu NEGATIVE_RESULT."""
    if n_quotable_pairs < W4_MIN_PAIRS or t_w4 is None:
        return {"verdict": "NEGATIVE_RESULT", "rule": "S1",
                "sigma_tot": None}
    if sigma_sampler is None or sigma_fss is None:
        return {"verdict": "NEGATIVE_RESULT", "rule": "S2",
                "sigma_tot": None}
    s_fss = max(float(sigma_fss), W4_SIGMA_FSS_FLOOR)
    s_tot = math.hypot(float(sigma_sampler), s_fss)
    if s_tot > W4_SIGMA_MAX:
        return {"verdict": "NEGATIVE_RESULT", "rule": "S2",
                "sigma_tot": s_tot}
    lo, hi = t_w4 - 2.0 * s_tot, t_w4 + 2.0 * s_tot
    consistent = hi >= W4_BAND[0] and lo <= W4_BAND[1]
    return {"verdict": "CONSISTENT" if consistent else "INCONSISTENT",
            "rule": "C" if consistent else "I", "sigma_tot": s_tot,
            "interval": [lo, hi]}


# ============================================================================
# 1. Kernel: innere Update-Schleifen (Numba-kompiliert, sonst reines Python)
# ============================================================================

def _nbr_array(nbr_list) -> tuple[np.ndarray, np.ndarray]:
    """Adjazenz-Tupel -> (nbr[n, zmax] int64, deg[n] int64), Reihenfolge je
    Knoten exakt wie im Tupel (Summationsreihenfolge = Original)."""
    n = len(nbr_list)
    zmax = max(len(t) for t in nbr_list)
    nbr = np.zeros((n, zmax), dtype=np.int64)
    deg = np.zeros(n, dtype=np.int64)
    for i, t in enumerate(nbr_list):
        deg[i] = len(t)
        nbr[i, :len(t)] = t
    return nbr, deg


def _py_anneal_sweep(th, nbr, deg, idx, news, us, E, ba):
    """Ein Anneal-Sweep, Arithmetik 1:1 wie PHY041 (_local_delta)."""
    for k in range(idx.shape[0]):
        i = idx[k]
        old = th[i]
        new = news[k]
        d = 0.0
        for m in range(deg[i]):
            tj = th[nbr[i, m]]
            d -= math.cos(new - tj) - math.cos(old - tj)
        if d < 0 or us[k] < math.exp(-ba * d):
            th[i] = new
            E += d
    return E


def _py_wl_sweep(th, nbr, deg, idx, news, us, E, b, lng, H, lnf, E_lo, binw,
                 nbins):
    """Ein WL-Sweep (Akzeptanz + ln g/H-Update), 1:1 wie PHY041."""
    for k in range(idx.shape[0]):
        i = idx[k]
        old = th[i]
        new = news[k]
        d = 0.0
        for m in range(deg[i]):
            tj = th[nbr[i, m]]
            d -= math.cos(new - tj) - math.cos(old - tj)
        Enew = E + d
        bnew = int((Enew - E_lo) // binw)
        if 0 <= bnew < nbins and math.log(us[k] + 1e-300) <= lng[b] - lng[bnew]:
            th[i] = new
            E = Enew
            b = bnew
        lng[b] += lnf
        H[b] += 1.0
    return E, b


def _py_prod_sweep(th, nbr, deg, idx, news, us, E, b, lng, E_lo, binw, nbins):
    """Ein Produktions-Sweep (g fix), 1:1 wie PHY041."""
    for k in range(idx.shape[0]):
        i = idx[k]
        old = th[i]
        new = news[k]
        d = 0.0
        for m in range(deg[i]):
            tj = th[nbr[i, m]]
            d -= math.cos(new - tj) - math.cos(old - tj)
        Enew = E + d
        bnew = int((Enew - E_lo) // binw)
        if 0 <= bnew < nbins and math.log(us[k] + 1e-300) <= lng[b] - lng[bnew]:
            th[i] = new
            E = Enew
            b = bnew
    return E, b


if HAVE_NUMBA:
    # fastmath AUS (IEEE-Semantik, keine Umordnung), kein parallel (seriell
    # wie das Original); cache=False, damit kein Stale-Cache die Gates
    # umgeht.
    _nb = numba.njit(cache=False, fastmath=False, nogil=True)
    _nb_anneal_sweep = _nb(_py_anneal_sweep)
    _nb_wl_sweep = _nb(_py_wl_sweep)
    _nb_prod_sweep = _nb(_py_prod_sweep)


def _sweep_fns(backend: str):
    """(anneal, wl, prod) fuer backend in {"auto", "numba", "python"}."""
    if backend not in ("auto", "numba", "python"):
        raise ValueError(f"unbekanntes backend {backend!r}")
    if backend == "numba" and not HAVE_NUMBA:
        raise RuntimeError("backend='numba' verlangt, numba nicht installiert")
    if backend == "python" or not HAVE_NUMBA:
        return _py_anneal_sweep, _py_wl_sweep, _py_prod_sweep
    return _nb_anneal_sweep, _nb_wl_sweep, _nb_prod_sweep


# ============================================================================
# 2. Treiber: Phasen getrennt (Anneal+WL | Produktion), Logik 1:1 PHY041
# ============================================================================

class WLState:
    """Zustand nach der WL-Phase (fuer Produktion und Spread-Zerlegung)."""

    def __init__(self, th, E, b, lng_arr, sweeps, sweeps_1t, max_drift,
                 E_lo, binw, nbins, rng):
        self.th = th
        self.E = E
        self.b = b
        self.lng = lng_arr
        self.sweeps = sweeps
        self.sweeps_1t = sweeps_1t
        self.max_drift = max_drift
        self.E_lo = E_lo
        self.binw = binw
        self.nbins = nbins
        self.rng = rng


def wl_phase(nbr_list, ei, ej, n: int, L: int, e_lo_ps: float,
             e_hi_ps: float, lnf_final: float = 1e-5, seed: int = 42,
             stream: int = _STREAM_WL, resync_every: int = 2000,
             backend: str = "auto", nbr_deg=None) -> WLState:
    """Anneal + WL-g(E) (Standard-WL -> echtes 1/t, B&P) - identisch zu
    PHY041 wl_entropic_lattice bis einschliesslich der WL-Schleife."""
    anneal, wl_sweep, _ = _sweep_fns(backend)
    nbr, deg = nbr_deg if nbr_deg is not None else _nbr_array(nbr_list)
    E_lo, E_hi = e_lo_ps * n, e_hi_ps * n
    nbins = int(round(E_hi - E_lo))
    binw = (E_hi - E_lo) / nbins
    rng = make_rng(seed, stream=stream + L)
    two_pi = 2.0 * math.pi

    th = rng.uniform(0, two_pi, n)
    E = _full_energy(th, ei, ej)
    ba = 0.0
    while not (E_lo <= E < E_hi):
        ba += 0.05
        if ba > 400.0:
            raise RuntimeError(
                f"Anneal erreicht Energie-Fenster [{E_lo},{E_hi}) nicht "
                f"(L={L}, E={E:.1f}); Fenster pruefen.")
        idx = rng.integers(0, n, size=n)
        news = rng.uniform(0.0, two_pi, size=n)
        us = rng.random(size=n)
        E = anneal(th, nbr, deg, idx, news, us, E, ba)
    b = min(max(int((E - E_lo) // binw), 0), nbins - 1)

    lng = np.zeros(nbins)
    H = np.zeros(nbins)
    lnf = 1.0
    use_1t = False
    sweeps = 0
    sweeps_1t = None
    max_drift = 0.0
    while lnf > lnf_final:
        idx = rng.integers(0, n, size=n)
        news = rng.uniform(0.0, two_pi, size=n)
        us = rng.random(size=n)
        E, b = wl_sweep(th, nbr, deg, idx, news, us, E, b, lng, H, lnf,
                        E_lo, binw, nbins)
        sweeps += 1
        if sweeps > 5_000_000:
            raise RuntimeError(
                f"WL-Phase konvergiert nicht (L={L}, sweeps={sweeps}, "
                f"lnf={lnf:.2e}); Fenster/Bin-Aufloesung pruefen.")
        if sweeps % resync_every == 0:
            E_exact = _full_energy(th, ei, ej)
            max_drift = max(max_drift, abs(E_exact - E))
            E = E_exact
            b = min(max(int((E - E_lo) // binw), 0), nbins - 1)
        if use_1t:
            lnf = 1.0 / sweeps
        elif H.min() >= 1:
            lnf *= 0.5
            H[:] = 0.0
            if lnf <= 1.0 / sweeps:
                use_1t = True
                sweeps_1t = sweeps
                lnf = 1.0 / sweeps
    lng_arr = lng.copy()
    lng_arr -= lng_arr.max()
    return WLState(th, E, b, lng_arr, sweeps, sweeps_1t, max_drift, E_lo,
                   binw, nbins, rng)


def production_phase(state: WLState, nbr_list, ei, ej, ax, ay, n: int,
                     L: int, prod_sweeps: int, resync_every: int = 2000,
                     backend: str = "auto", rng=None, nbr_deg=None
                     ) -> WLResult:
    """Produktion (g fix) mit mikrokanonischen Aggregaten, identisch PHY041.
    rng=None -> weiter mit dem RNG der WL-Phase (= Original); ein anderer
    rng (Spread-Zerlegung K3) startet eine unabhaengige Produktion auf
    einer KOPIE des Zustands."""
    _, _, prod_sweep = _sweep_fns(backend)
    nbr, deg = nbr_deg if nbr_deg is not None else _nbr_array(nbr_list)
    independent = rng is not None
    rng = state.rng if rng is None else rng
    th = state.th.copy() if independent else state.th
    E, b = state.E, state.b
    E_lo, binw, nbins = state.E_lo, state.binw, state.nbins
    lngl = state.lng
    max_drift = state.max_drift
    two_pi = 2.0 * math.pi

    acc_x = {k: np.zeros(nbins) for k in OBS}
    acc_y = {k: np.zeros(nbins) for k in OBS}
    cnt = np.zeros(nbins)
    for sw in range(prod_sweeps):
        idx = rng.integers(0, n, size=n)
        news = rng.uniform(0.0, two_pi, size=n)
        us = rng.random(size=n)
        E, b = prod_sweep(th, nbr, deg, idx, news, us, E, b, lngl, E_lo,
                          binw, nbins)
        if (sw + 1) % resync_every == 0:
            E_exact = _full_energy(th, ei, ej)
            max_drift = max(max_drift, abs(E_exact - E))
            E = E_exact
            b = min(max(int((E - E_lo) // binw), 0), nbins - 1)
        phi = th[ei] - th[ej]
        cph = np.cos(phi)
        sph = np.sin(phi)
        ox = _aggregates_dir(cph, sph, ax)
        oy = _aggregates_dir(cph, sph, ay)
        for k in OBS:
            acc_x[k][b] += ox[k]
            acc_y[k][b] += oy[k]
        cnt[b] += 1
    if not independent:
        state.E, state.b, state.max_drift = E, b, max_drift
    centers = E_lo + (np.arange(nbins) + 0.5) * binw
    mask = cnt > 0
    denom = np.maximum(cnt, 1)
    micro_x = {k: np.where(mask, acc_x[k] / denom, 0.0) for k in OBS}
    micro_y = {k: np.where(mask, acc_y[k] / denom, 0.0) for k in OBS}
    return WLResult(L=L, n=n, centers=centers, lng=state.lng.copy(),
                    mask=mask, micro_x=micro_x, micro_y=micro_y,
                    wl_sweeps=state.sweeps, prod_sweeps=prod_sweeps,
                    one_over_t_engaged=state.sweeps_1t is not None,
                    sweeps_at_1t=state.sweeps_1t or 0)


def wl_entropic_fast(nbr_list, ei, ej, ax, ay, n: int, L: int,
                     e_lo_ps: float, e_hi_ps: float,
                     lnf_final: float = 1e-5, prod_sweeps: int = 30000,
                     seed: int = 42, stream: int = _STREAM_WL,
                     resync_every: int = 2000, backend: str = "auto"
                     ) -> WLResult:
    """Drop-in fuer PHY041 wl_entropic_lattice (gleiche Signatur bis auf
    verbose/backend; bit-identisches Ergebnis - Gates VAL-BIT-S/-P)."""
    nd = _nbr_array(nbr_list)
    st = wl_phase(nbr_list, ei, ej, n, L, e_lo_ps, e_hi_ps, lnf_final, seed,
                  stream, resync_every, backend, nbr_deg=nd)
    return production_phase(st, nbr_list, ei, ej, ax, ay, n, L, prod_sweeps,
                            resync_every, backend, nbr_deg=nd)


def results_identical(a, b) -> bool:
    """Bit-Gleichheit zweier WLResult (lng, Maske, Zentren, alle Aggregate,
    Sweep-Zaehler)."""
    if (a.wl_sweeps != b.wl_sweeps or a.prod_sweeps != b.prod_sweeps
            or a.n != b.n or a.L != b.L):
        return False
    if not (np.array_equal(a.lng, b.lng) and np.array_equal(a.mask, b.mask)
            and np.array_equal(a.centers, b.centers)):
        return False
    return all(np.array_equal(a.micro_x[k], b.micro_x[k])
               and np.array_equal(a.micro_y[k], b.micro_y[k]) for k in OBS)


# ============================================================================
# 3. Kalibrier-Jobs (top-level, picklebar; RNG rein (seed, stream)-bestimmt)
# ============================================================================

def _prod_rule(nbins: int, L: int, factor: int = 1) -> int:
    """PHY042-Rezept max(30000, 60*nbins) (L=24: 30000), optional x factor."""
    return factor * _p42._prod_sweeps_for(nbins, L)


def _cal_stream(walker: int) -> int:
    return _STREAM_CAL_BASE + 10_000 * walker


def _job(args: tuple) -> dict:
    """Ein Kalibrier-Job. kind:
    "phy042"  - PHY042-Walker mit Original-Stream + committetem Fenster
                (VAL-BIT-P);
    "cal"     - K1/K2/K3-Walker (lnf_final, prod-Faktor), eigener Stream;
    "decomp"  - K3-Zerlegung: EINE g(E) (Walker 0), n_prod unabhaengige
                Produktionen."""
    kind, L, w, e_lo, e_hi, opts = args
    _, nbr_list, ei, ej, ax, ay, n = honeycomb_arrays(L)
    nbins = int(round((e_hi - e_lo) * n))
    nd = _nbr_array(nbr_list)
    seed = opts.get("seed", 42)
    lnf_final = opts.get("lnf_final", 1e-5)
    backend = opts.get("backend", "auto")
    if kind == "phy042":
        stream = _p42._stream_for(w)
    else:
        stream = _cal_stream(w)
    prod = _prod_rule(nbins, L, opts.get("prod_factor", 1))
    t0 = time.perf_counter()
    st = wl_phase(nbr_list, ei, ej, n, L, e_lo, e_hi, lnf_final, seed,
                  stream, backend=backend, nbr_deg=nd)
    t_wl = time.perf_counter() - t0
    out = {"kind": kind, "L": L, "walker": w, "n": n, "nbins": nbins,
           "stream": stream, "lnf_final": lnf_final, "prod_sweeps": prod,
           "prod_factor": opts.get("prod_factor", 1),
           "wl_sweeps": st.sweeps, "sweeps_1t": st.sweeps_1t,
           "t_wl_s": t_wl, "backend": backend if backend != "auto" else
           ("numba" if HAVE_NUMBA else "python")}
    if kind == "decomp":
        curves, unc = [], []
        t1 = time.perf_counter()
        for p in range(opts["n_prod"]):
            rng = make_rng(seed, stream=_STREAM_DECOMP_BASE + 10_000 * p + L)
            res = production_phase(st, nbr_list, ei, ej, ax, ay, n, L, prod,
                                   backend=backend, rng=rng, nbr_deg=nd)
            curves.append(upsilon_curves(res, _p42._T_GRID))
            unc.append(_uncovered_profile(res))
        out["t_prod_s"] = time.perf_counter() - t1
        out["curves"] = [_curve_dict(c) for c in curves]
        out["uncovered_mass"] = unc
        return out
    t1 = time.perf_counter()
    res = production_phase(st, nbr_list, ei, ej, ax, ay, n, L, prod,
                           backend=backend, nbr_deg=nd)
    out["t_prod_s"] = time.perf_counter() - t1
    out["covered"] = [int(res.mask.sum()), int(len(res.mask))]
    out["leak_max"] = max(_p41.canonical_edge_leak(res, float(T))
                          for T in _p42._T_GRID)
    out["uncovered_mass"] = _uncovered_profile(res)
    out["curve"] = _curve_dict(upsilon_curves(res, _p42._T_GRID))
    return out


def _uncovered_profile(res) -> list:
    """Kanonische Masse auf produktions-unbesetzten Bins je T des Gitters
    (PHY042-Vertrag `_uncovered_mass`; Gate bindet nur in-Domaene)."""
    return [_p42._uncovered_mass(res, float(T)) for T in _p42._T_GRID]


def _curve_dict(c: dict) -> dict:
    return {"T": c["T"].tolist(), "y2": c["y2"].tolist(),
            "y4_scaled": c["y4_scaled"].tolist(), "E": c["E"].tolist()}


def _run_jobs(jobs: list, max_workers: int) -> list:
    jobs = sorted(jobs, key=lambda j: (-j[1], j[2]))   # grosse L zuerst
    if max_workers <= 1:
        return [_job(j) for j in jobs]
    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        return list(ex.map(_job, jobs))


# ============================================================================
# 4. Auswertung (blind)
# ============================================================================

def spread_domain(y2_curves: list, t_grid) -> dict:
    """PHY042-Semantik (_walker_domain), JSON-tauglich."""
    d = _p42._walker_domain(y2_curves, np.asarray(t_grid))
    return _p42._domain_status(d)


def pair_crossing_spread(curves_a: list, La: int, curves_b: list, Lb: int,
                         t_grid) -> dict:
    """BLIND: Walker-Streuung (max-min) der C-eliminierten Paar-Crossings
    ueber alle Walker-Kombinationen - OHNE Lage-Wert. Rueckgabe: n_combos,
    n_crossing, spread (None bei < 2 Crossings)."""
    t = np.asarray(t_grid)
    vals = []
    for ca in curves_a:
        for cb in curves_b:
            tb = tbkt_pair_from_curves(t, np.asarray(ca), La,
                                       np.asarray(cb), Lb)
            if tb is not None:
                vals.append(tb)
    return {"n_combos": len(curves_a) * len(curves_b),
            "n_crossing": len(vals),
            "spread": (max(vals) - min(vals)) if len(vals) >= 2 else None}


def fit_power_law(Ls, ys) -> dict:
    """Kleinste Quadrate in log-log: y = A * L^p. Rueckgabe p, A und die
    lokalen Exponenten benachbarter L (Kruemmungs-Diagnostik)."""
    x = np.log(np.asarray(Ls, dtype=float))
    y = np.log(np.asarray(ys, dtype=float))
    p, lnA = np.polyfit(x, y, 1)
    local = [float((y[i + 1] - y[i]) / (x[i + 1] - x[i]))
             for i in range(len(x) - 1)]
    return {"p": float(p), "A": float(math.exp(lnA)), "local_p": local}


def assert_blind(obj, path: str = "") -> None:
    """Fails-closed: kein Lage-Schluessel irgendwo im Report."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            if any(f in str(k) for f in BLIND_FORBIDDEN_KEYS):
                raise AssertionError(f"Blind-Vertrag verletzt: {path}/{k}")
            assert_blind(v, f"{path}/{k}")
    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            assert_blind(v, f"{path}[{i}]")


def _clean(o):
    if o is None:
        return None
    if isinstance(o, (np.bool_, bool)):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, (np.floating, float)):
        v = float(o)
        return v if math.isfinite(v) else None
    if isinstance(o, np.ndarray):
        return [_clean(x) for x in o.tolist()]
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(x) for x in o]
    return o


# ============================================================================
# 5. Stufen (je Stufe ein JSON-Teilergebnis; `report` fuegt zusammen)
# ============================================================================

PHY042_REPORT = _ROOT / _p42.PHY042_REPORT_V01
LADDER_LS = (24, 32, 48, 64)
LADDER_WALKERS = 3            # W4-Vertrag: >= 3 Walker an JEDEM L
LEVER_L = 48


def _windows_for(Ls, master_seed: int = 42) -> dict:
    """Fenster je L exakt nach PHY042-Verfahren (Wolff-Anker T=0.50/0.70)."""
    out = {}
    for L in Ls:
        m_lo, s_lo = _p41.wolff_anchor(L, 0.50, master_seed=master_seed)
        m_hi, s_hi = _p41.wolff_anchor(L, 0.70, master_seed=master_seed)
        out[L] = list(_p41.window_from_anchors(m_lo, s_lo, m_hi, s_hi))
    return out


CROSS_PLATFORM_RTOL = 1e-12   # ULP-Ebene der numpy-Reduktionen


def _max_rel_diff(a, b) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return float(np.max(np.abs(a - b) / np.maximum(np.abs(b), 1e-300)))


def stage_valbit(max_workers: int = 4) -> dict:
    """VAL-BIT-P: alle 7 committeten PHY042-Walker (L=24/32/48) mit dem
    Numba-Kernel - wl_sweeps und gespeicherte Y2/Y4-Kurven EXAKT."""
    raw = PHY042_REPORT.read_bytes().replace(b"\r\n", b"\n")
    rep = json.loads(raw)
    jobs = []
    for L in rep["Ls"]:
        e_lo, e_hi = rep["windows"][str(L)]["window_ps"]
        for w in range(int(rep["n_walkers"][str(L)])):
            jobs.append(("phy042", int(L), w, e_lo, e_hi,
                         {"backend": "numba"}))
    t0 = time.perf_counter()
    outs = _run_jobs(jobs, max_workers)
    wall = time.perf_counter() - t0
    rows = []
    for o in outs:
        key = f"{o['L']}_{o['walker']}"
        c = rep["curves"][key]
        dy2 = _max_rel_diff(o["curve"]["y2"], c["y2"])
        dy4 = _max_rel_diff(o["curve"]["y4_scaled"], c["y4_scaled"])
        rows.append({"key": key, "wl_sweeps": o["wl_sweeps"],
                     "sweeps_1t": o["sweeps_1t"],
                     "wl_sweeps_committed": rep["wl_sweeps"][key],
                     "wl_sweeps_equal": o["wl_sweeps"] == rep["wl_sweeps"][key],
                     "covered_equal": o["covered"] == rep["covered"][key],
                     "y2_bit_equal": o["curve"]["y2"] == c["y2"],
                     "y4_bit_equal": o["curve"]["y4_scaled"] == c["y4_scaled"],
                     "y2_max_rel_diff": dy2, "y4_max_rel_diff": dy4,
                     "t_wl_s": o["t_wl_s"], "t_prod_s": o["t_prod_s"],
                     "n": o["n"], "prod_sweeps": o["prod_sweeps"]})
    # Gate-Semantik (Befund 2026-09-26): gegen die auf ANDERER Plattform
    # committete Evidenz sind die Trajektorien identisch (wl_sweeps +
    # Belegung exakt - jede Abweichung im Sampler waere chaotisch sichtbar),
    # die Kurven aber nur bis auf die numpy-Reduktionsreihenfolge (np.dot/
    # np.sum, ULP-Ebene). Bit-Gleichheit Original vs Numba wird auf DERSELBEN
    # Maschine separat gepinnt (Stufe speed + Tests).
    ok = all(r["wl_sweeps_equal"] and r["covered_equal"]
             and r["y2_max_rel_diff"] <= CROSS_PLATFORM_RTOL
             and r["y4_max_rel_diff"] <= CROSS_PLATFORM_RTOL for r in rows)
    return {"stage": "valbit", "source_sha256":
            hashlib.sha256(raw).hexdigest(), "jobs": sorted(
                rows, key=lambda r: r["key"]), "wall_s": wall,
            "max_workers": max_workers, "rtol": CROSS_PLATFORM_RTOL,
            "PASS_VAL_BIT_PHY042": ok}


def stage_speed() -> dict:
    """Ein identischer Produktions-Job (PHY042 L=24 Walker 0) einmal mit dem
    Original-Python-Kernel, einmal mit Numba - Speedup auf echter Groesse,
    plus Bit-Gleichheit dieses Jobs."""
    rep = json.loads(PHY042_REPORT.read_text(encoding="utf-8"))
    L = 24
    e_lo, e_hi = rep["windows"][str(L)]["window_ps"]
    _, nbr_list, ei, ej, ax, ay, n = honeycomb_arrays(L)
    nbins = int(round((e_hi - e_lo) * n))
    prod = _p42._prod_sweeps_for(nbins, L)
    t0 = time.perf_counter()
    ref = _p41.wl_entropic_lattice(nbr_list, ei, ej, ax, ay, n, L, e_lo, e_hi,
                                   lnf_final=1e-5, prod_sweeps=prod,
                                   seed=42, stream=_STREAM_WL, verbose=False)
    t_py = time.perf_counter() - t0
    wl_entropic_fast(nbr_list, ei, ej, ax, ay, n, L, e_lo, e_hi,
                     lnf_final=0.5, prod_sweeps=1, backend="numba")  # JIT
    t0 = time.perf_counter()
    fast = wl_entropic_fast(nbr_list, ei, ej, ax, ay, n, L, e_lo, e_hi,
                            lnf_final=1e-5, prod_sweeps=prod, seed=42,
                            stream=_STREAM_WL, backend="numba")
    t_nb = time.perf_counter() - t0
    updates = (ref.wl_sweeps + prod) * n
    return {"stage": "speed", "L": L, "n": n, "wl_sweeps": ref.wl_sweeps,
            "prod_sweeps": prod, "spin_updates": updates,
            "t_python_s": t_py, "t_numba_s": t_nb,
            "speedup": t_py / t_nb,
            "python_updates_per_s": updates / t_py,
            "numba_updates_per_s": updates / t_nb,
            "PASS_BIT_IDENTICAL_L24": results_identical(ref, fast)}


def stage_ladder(Ls=LADDER_LS, n_walkers: int = LADDER_WALKERS,
                 max_workers: int = 4) -> dict:
    """K1 Kosten + K2 Praezision beim PHY042-Budget-Rezept."""
    windows = _windows_for(Ls)
    jobs = [("cal", L, w, windows[L][0], windows[L][1],
             {"backend": "numba"})
            for L in Ls for w in range(n_walkers)]
    t0 = time.perf_counter()
    outs = _run_jobs(jobs, max_workers)
    return {"stage": "ladder", "windows": windows, "jobs": outs,
            "wall_s": time.perf_counter() - t0, "max_workers": max_workers}


def stage_levers(L: int = LEVER_L, n_walkers: int = 3, n_prod: int = 3,
                 max_workers: int = 4) -> dict:
    """K3 Hebel bei L=48: 4x Produktion, lnf_final 1e-6, Spread-Zerlegung."""
    e_lo, e_hi = _windows_for((L,))[L]
    jobs = [("cal", L, w, e_lo, e_hi, {"backend": "numba", "prod_factor": 4})
            for w in range(n_walkers)]
    jobs += [("cal", L, w, e_lo, e_hi, {"backend": "numba",
                                        "lnf_final": 1e-6})
             for w in range(n_walkers)]
    jobs.append(("decomp", L, 0, e_lo, e_hi, {"backend": "numba",
                                              "n_prod": n_prod}))
    t0 = time.perf_counter()
    outs = _run_jobs(jobs, max_workers)
    return {"stage": "levers", "L": L, "window": [e_lo, e_hi], "jobs": outs,
            "wall_s": time.perf_counter() - t0, "max_workers": max_workers}


REPORT_STEM = "260926 PHY044 honeycomb wl calibration report"


def _group(jobs: list, **match) -> list:
    return [j for j in jobs if all(j.get(k) == v for k, v in match.items())]


def _walker_block(jobs: list, t_grid) -> dict:
    """Domaene + Spread-Profil (blind) einer Walker-Gruppe gleichen L, dazu
    die unbesetzte kanonische Masse IN der Domaene (PHY042-Gate), die
    Bin-Belegung und ob die 1/t-Phase je Walker gegriffen hat."""
    jobs = sorted(jobs, key=lambda j: j["walker"])
    y2 = [j["curve"]["y2"] for j in jobs]
    d = _p42._walker_domain(y2, np.asarray(t_grid))
    dom = _p42._domain_status(d)
    spread = (np.max(y2, axis=0) - np.min(y2, axis=0)) if len(y2) > 1 else None
    unc = [max((u for u, m in zip(j["uncovered_mass"], d["mask"]) if m),
               default=0.0) for j in jobs]
    return {"n_walkers": len(y2), "domain": dom,
            "meets_T_req": (dom["tmax"] is not None
                            and dom["tmax"] >= W4_T_REQ - 1e-12),
            "spread_profile": None if spread is None else spread.tolist(),
            "uncovered_mass_in_domain_max": float(max(unc)),
            "covered_fraction": [j["covered"][0] / j["covered"][1]
                                 for j in jobs],
            "one_over_t_engaged": [j["sweeps_1t"] is not None for j in jobs],
            "sweeps_1t": [j["sweeps_1t"] for j in jobs],
            "wl_sweeps": [j["wl_sweeps"] for j in jobs]}


def analyse(stages: dict) -> dict:
    """Blinde Auswertung K1-K3 + Gates aus den Stufen-JSONs."""
    t_grid = [float(t) for t in _p42._T_GRID]
    out: dict = {"module": "PHY044_honeycomb_wl_calibration",
                 "attribution": "Coworker Research / Coworkerz",
                 "date": "2026-09-26",
                 "spec": "spec/260926 PHI HEX w4 honeycomb preregistration "
                         "v01.md (§5)",
                 "blind_contract": "no pair-crossing location is computed "
                                   "into, printed to or stored in this "
                                   "report; only spreads",
                 "T_grid": t_grid, "T_req": W4_T_REQ,
                 "domain_threshold": W4_DOMAIN_THRESHOLD}
    gates = {}
    vb = stages.get("valbit")
    if vb:
        out["val_bit_phy042"] = {k: vb[k] for k in (
            "source_sha256", "jobs", "wall_s", "max_workers", "rtol")}
        gates["PASS_VAL_BIT_PHY042_ALL_7_WALKERS"] = bool(
            vb["PASS_VAL_BIT_PHY042"] and len(vb["jobs"]) == 7)
    sp = stages.get("speed")
    if sp:
        out["speed"] = {k: sp[k] for k in sp if k not in (
            "stage", "PASS_BIT_IDENTICAL_L24")}
        gates["PASS_BIT_IDENTICAL_L24_SPEED_JOB"] = bool(
            sp["PASS_BIT_IDENTICAL_L24"])
    lad = stages.get("ladder")
    if lad:
        jobs = lad["jobs"]
        Ls = sorted({j["L"] for j in jobs})
        k1 = {}
        for L in Ls:
            g = _group(jobs, L=L)
            upd = [(j["wl_sweeps"] + j["prod_sweeps"]) * j["n"] for j in g]
            wall = [j["t_wl_s"] + j["t_prod_s"] for j in g]
            k1[str(L)] = {
                "n": g[0]["n"], "nbins": g[0]["nbins"],
                "prod_sweeps": g[0]["prod_sweeps"],
                "wl_sweeps": [j["wl_sweeps"] for j in g],
                "sweeps_1t": [j["sweeps_1t"] for j in g],
                "spin_updates_mean": float(np.mean(upd)),
                "wall_numba_s_mean": float(np.mean(wall)),
                "wall_numba_s_max": float(np.max(wall))}
        fit_upd = fit_power_law(Ls, [k1[str(L)]["spin_updates_mean"]
                                     for L in Ls])
        fit_wall = fit_power_law(Ls, [k1[str(L)]["wall_numba_s_mean"]
                                      for L in Ls])
        fit_wl = fit_power_law(Ls, [float(np.mean(k1[str(L)]["wl_sweeps"]))
                                    for L in Ls])
        for L in Ls:
            g = _group(jobs, L=L)
            k1[str(L)]["wall_wl_s_mean"] = float(np.mean(
                [j["t_wl_s"] for j in g]))
            k1[str(L)]["wall_prod_s_mean"] = float(np.mean(
                [j["t_prod_s"] for j in g]))
        fit_wall_wl = fit_power_law(Ls, [k1[str(L)]["wall_wl_s_mean"]
                                         for L in Ls])
        fit_wall_prod = fit_power_law(Ls, [k1[str(L)]["wall_prod_s_mean"]
                                           for L in Ls])
        proj = {}
        L_last = Ls[-1]
        for Lp in (96, 128):
            # Projektion vom groessten gemessenen L mit dem LOKALEN Exponenten
            # des letzten Intervalls (konservativer als der Gesamt-Fit, wenn
            # die Kurve steiler wird) und mit dem Gesamt-Fit
            p_loc = fit_wall["local_p"][-1]
            w_last = k1[str(L_last)]["wall_numba_s_mean"]
            w_loc = w_last * (Lp / L_last) ** p_loc
            w_fit = fit_wall["A"] * Lp ** fit_wall["p"]
            row = {"wall_numba_h_local_exp": w_loc / 3600.0,
                   "wall_numba_h_global_fit": w_fit / 3600.0,
                   "local_exponent_used": p_loc}
            if sp:
                upd_loc = (k1[str(L_last)]["spin_updates_mean"]
                           * (Lp / L_last) ** fit_upd["local_p"][-1])
                row["wall_python_h_model"] = (
                    upd_loc / sp["python_updates_per_s"] / 3600.0)
            proj[str(Lp)] = row
        out["K1_cost"] = {"per_L": k1, "fit_spin_updates": fit_upd,
                          "fit_wall_numba": fit_wall,
                          "fit_wall_wl_numba": fit_wall_wl,
                          "fit_wall_prod_numba": fit_wall_prod,
                          "fit_wl_sweeps": fit_wl, "projection": proj,
                          "wall_s": lad["wall_s"],
                          "max_workers": lad["max_workers"],
                          "windows": lad["windows"]}
        k2 = {str(L): _walker_block(_group(jobs, L=L), t_grid) for L in Ls}
        pairs = {}
        for i, La in enumerate(Ls):
            for Lb in Ls[i + 1:]:
                ca = [j["curve"]["y2"] for j in _group(jobs, L=La)]
                cb = [j["curve"]["y2"] for j in _group(jobs, L=Lb)]
                pairs[f"{La}_{Lb}"] = pair_crossing_spread(ca, La, cb, Lb,
                                                           t_grid)
        out["K2_precision"] = {"per_L": k2, "pair_crossing_spread": pairs}
        gates["PASS_WALKER_CONTRACT_GE3_EVERY_L"] = all(
            k2[str(L)]["n_walkers"] >= W4_MIN_WALKERS for L in Ls)
        gates["PASS_NO_EDGE_LEAK_LADDER"] = all(
            j["leak_max"] < 1e-3 for j in jobs)
        # PHY042-Vertrag: unbesetzte kanonische Masse in-Domaene < 1e-3.
        # (Hier stand zunaechst "volle Bin-Belegung" - strenger als jeder
        # Vertrag und fuer Rand-Bins ohne kanonisches Gewicht falsch
        # gestellt; die Belegung wird als Befund berichtet, nicht gegated.)
        gates["PASS_NO_UNCOVERED_MASS_IN_DOMAIN_LADDER"] = all(
            k2[str(L)]["uncovered_mass_in_domain_max"] < 1e-3 for L in Ls)
        # 1/t-Eintritt ist eine Eigenschaft des BUDGETS (das misst die
        # Kalibrierung), kein Pipeline-Defekt -> Befund, kein Gate.
        out["findings"] = {"one_over_t_engaged_per_L": {
            str(L): k2[str(L)]["one_over_t_engaged"] for L in Ls}}
    lev = stages.get("levers")
    if lev and lad:
        L = lev["L"]
        jobs = lev["jobs"]
        base = _walker_block(_group(lad["jobs"], L=L), t_grid)
        p4 = _walker_block(_group(jobs, kind="cal", prod_factor=4), t_grid)
        l6 = _walker_block([j for j in jobs if j["kind"] == "cal"
                            and j["lnf_final"] == 1e-6], t_grid)
        dec = _group(jobs, kind="decomp")[0]
        dec_y2 = [c["y2"] for c in dec["curves"]]
        dec_block = {"n_productions": len(dec_y2),
                     "wl_sweeps": dec["wl_sweeps"],
                     "sweeps_1t": dec["sweeps_1t"],
                     "domain": spread_domain(dec_y2, t_grid),
                     "spread_profile": (np.max(dec_y2, axis=0)
                                        - np.min(dec_y2, axis=0)).tolist()}
        variants = {"baseline_1x_lnf1e-5": base, "prod_x4_lnf1e-5": p4,
                    "prod_x1_lnf1e-6": l6}
        cost = {}
        for name, grp in (("baseline_1x_lnf1e-5", _group(lad["jobs"], L=L)),
                          ("prod_x4_lnf1e-5", _group(jobs, kind="cal",
                                                     prod_factor=4)),
                          ("prod_x1_lnf1e-6", [j for j in jobs
                                               if j["kind"] == "cal"
                                               and j["lnf_final"] == 1e-6])):
            cost[name] = {
                "spin_updates_mean": float(np.mean(
                    [(j["wl_sweeps"] + j["prod_sweeps"]) * j["n"]
                     for j in grp])),
                "wall_numba_s_mean": float(np.mean(
                    [j["t_wl_s"] + j["t_prod_s"] for j in grp]))}
        meeting = [n for n in variants if variants[n]["meets_T_req"]]
        cheapest = (min(meeting, key=lambda n: cost[n]["spin_updates_mean"])
                    if meeting else None)
        # Anteil des Produktions-Rauschens am Walker-Spread (feste g(E) vs
        # unabhaengige g(E)+Produktion), je T bis T_req
        k_req = int(np.argmin(np.abs(np.asarray(t_grid) - W4_T_REQ)))
        prod_share = (float(np.max(dec_block["spread_profile"][:k_req + 1]))
                      / float(np.max(base["spread_profile"][:k_req + 1])))
        dec_block["max_spread_upto_T_req_ratio_vs_baseline"] = prod_share
        recipe = None
        if cheapest is not None and "K1_cost" in out:
            fac = 4 if cheapest.startswith("prod_x4") else 1
            k1c = out["K1_cost"]
            L_last = max(int(x) for x in k1c["per_L"])
            r_last = k1c["per_L"][str(L_last)]
            p_wl = k1c["fit_wall_wl_numba"]["local_p"][-1]
            p_pr = k1c["fit_wall_prod_numba"]["local_p"][-1]
            proj = {}
            for Lp in (48, 64, 96, 128):
                wl = r_last["wall_wl_s_mean"] * (Lp / L_last) ** p_wl
                pr = fac * r_last["wall_prod_s_mean"] * (Lp / L_last) ** p_pr
                proj[str(Lp)] = {"wall_numba_h_per_walker": (wl + pr) / 3600,
                                 "measured": Lp == L}
            proj["48"]["wall_numba_h_per_walker_measured"] = (
                cost[cheapest]["wall_numba_s_mean"] / 3600)
            recipe = {"variant": cheapest, "prod_factor": fac,
                      "lnf_final": 1e-5,
                      "local_exponents": {"wl": p_wl, "prod": p_pr},
                      "projection": proj,
                      "note": "L=64/96/128 extrapoliert (K1-Exponenten des "
                              "letzten Intervalls); T_req muss im W4-Lauf "
                              "selbst je L gemessen werden (Spec §4.4)"}
        out["K3_levers"] = {"L": L, "variants": variants, "cost": cost,
                            "decomposition_fixed_gE": dec_block,
                            "budget_rule": "smallest calibrated budget with "
                                           "T_max(L) >= T_req (Spec §5)",
                            "cheapest_meeting_T_req": cheapest,
                            "w4_recipe": recipe,
                            "w4_go_L_ge_48": cheapest is not None}
        gates["PASS_NO_EDGE_LEAK_LEVERS"] = all(
            j["leak_max"] < 1e-3 for j in jobs if j["kind"] == "cal")
        gates["PASS_NO_UNCOVERED_MASS_IN_DOMAIN_LEVERS"] = all(
            v["uncovered_mass_in_domain_max"] < 1e-3
            for v in variants.values())
    out["pass_gates"] = gates
    out["overall_pass"] = bool(gates) and all(gates.values())
    assert_blind(out)
    return out


def _fmt(x, f=".4f"):
    return "null" if x is None else format(x, f)


def write_text_report(rep: dict, path: Path) -> None:
    L_ = []
    L_.append("PHY044 - W4-Kalibrierung honeycomb Wang-Landau (BLIND)")
    L_.append("Coworker Research / Coworkerz, 2026-09-26")
    L_.append("=" * 70)
    L_.append(f"Spec: {rep['spec']}")
    L_.append(f"Blind-Vertrag: {rep['blind_contract']}")
    L_.append(f"T_req = {rep['T_req']}, Domaenen-Schwelle Walker-Spread "
              f"< {rep['domain_threshold']}")
    if "val_bit_phy042" in rep:
        vb = rep["val_bit_phy042"]
        L_.append("")
        L_.append("--- VAL-BIT-P: Numba reproduziert committeten PHY042 ---")
        L_.append(f"Quelle sha256 {vb['source_sha256']}")
        L_.append(f"Trajektorie exakt (wl_sweeps, Belegung); Kurven bis "
                  f"rtol {vb['rtol']:g} (numpy-Reduktionen, Plattform)")
        for j in vb["jobs"]:
            L_.append(f"  {j['key']:>5}: wl_sweeps {j['wl_sweeps']} "
                      f"(committed {j['wl_sweeps_committed']}) "
                      f"Belegung {'==' if j['covered_equal'] else '!='} "
                      f"dY2 {j['y2_max_rel_diff']:.1e} "
                      f"dY4 {j['y4_max_rel_diff']:.1e} "
                      f"t={j['t_wl_s'] + j['t_prod_s']:.0f}s")
        L_.append(f"  Wall gesamt {vb['wall_s']:.0f}s auf "
                  f"{vb['max_workers']} Prozessen (PHY042-Original 2026-07-06:"
                  " ~3498 s auf 4)")
    if "speed" in rep:
        s = rep["speed"]
        L_.append("")
        L_.append("--- Speedup auf identischem Produktions-Job (L=24) ---")
        L_.append(f"  Python {s['t_python_s']:.1f}s | Numba "
                  f"{s['t_numba_s']:.1f}s | Speedup {s['speedup']:.1f}x | "
                  f"{s['python_updates_per_s']:.3g} vs "
                  f"{s['numba_updates_per_s']:.3g} Spin-Updates/s")
    if "K1_cost" in rep:
        k = rep["K1_cost"]
        L_.append("")
        L_.append("--- K1 Kosten je Walker (PHY042-Rezept, Numba) ---")
        L_.append("     L      n  nbins   wl_sweeps (3 Walker)        "
                  "prod   Updates    Wall[s]")
        for L, r in k["per_L"].items():
            L_.append(f"  {L:>4} {r['n']:>6} {r['nbins']:>6}   "
                      f"{str(r['wl_sweeps']):<26} {r['prod_sweeps']:>7} "
                      f"{r['spin_updates_mean']:.3e} "
                      f"{r['wall_numba_s_mean']:>8.0f}")
        f = k["fit_spin_updates"]
        L_.append(f"  Fit Spin-Updates ~ L^{f['p']:.2f} (lokal: "
                  f"{', '.join(f'{p:.2f}' for p in f['local_p'])})")
        f = k["fit_wl_sweeps"]
        L_.append(f"  Fit wl_sweeps    ~ L^{f['p']:.2f} (lokal: "
                  f"{', '.join(f'{p:.2f}' for p in f['local_p'])})")
        f = k["fit_wall_numba"]
        L_.append(f"  Fit Wall (Numba) ~ L^{f['p']:.2f} (lokal: "
                  f"{', '.join(f'{p:.2f}' for p in f['local_p'])})")
        for Lp, r in k["projection"].items():
            py = r.get("wall_python_h_model")
            L_.append(f"  Projektion L={Lp}: Numba {r['wall_numba_h_local_exp']:.1f} h"
                      f" (lokaler Exp. {r['local_exponent_used']:.2f}) / "
                      f"{r['wall_numba_h_global_fit']:.1f} h (Gesamt-Fit)"
                      + ("" if py is None else f"; Python-Modell {py:.0f} h"))
    if "K2_precision" in rep:
        k = rep["K2_precision"]
        L_.append("")
        L_.append("--- K2 Praezision (Walker-Spread-Domaene, blind) ---")
        for L, b in k["per_L"].items():
            d = b["domain"]
            L_.append(f"  L={L:>3}: T_max={_fmt(d['tmax'])} "
                      f"(max Spread {_fmt(d['max_spread'])}) "
                      f"T_req erfuellt: {b['meets_T_req']}")
            L_.append(f"         1/t gegriffen {b['one_over_t_engaged']} "
                      f"(ab Sweep {b['sweeps_1t']}); Belegung "
                      f"{', '.join(f'{c:.2f}' for c in b['covered_fraction'])}"
                      f"; unbesetzte Masse in-Domaene "
                      f"{b['uncovered_mass_in_domain_max']:.1e}")
        L_.append("  Walker-Streuung der Paar-Crossings (max-min, KEINE Lage):")
        for key, r in k["pair_crossing_spread"].items():
            L_.append(f"    {key}: {r['n_crossing']}/{r['n_combos']} Kombos "
                      f"mit Crossing, Streuung {_fmt(r['spread'])}")
    if "K3_levers" in rep:
        k = rep["K3_levers"]
        L_.append("")
        L_.append(f"--- K3 Hebel bei L={k['L']} ---")
        for name, b in k["variants"].items():
            c = k["cost"][name]
            L_.append(f"  {name:<22} T_max={_fmt(b['domain']['tmax'])} "
                      f"max Spread {_fmt(b['domain']['max_spread'])} "
                      f"Updates {c['spin_updates_mean']:.3e} "
                      f"Wall {c['wall_numba_s_mean']:.0f}s "
                      f"T_req: {b['meets_T_req']}")
            L_.append(f"  {'':<22} 1/t {b['one_over_t_engaged']} "
                      f"wl_sweeps {b['wl_sweeps']}")
        d = k["decomposition_fixed_gE"]
        L_.append(f"  Zerlegung (1 g(E), {d['n_productions']} Produktionen): "
                  f"T_max={_fmt(d['domain']['tmax'])} "
                  f"max Spread {_fmt(d['domain']['max_spread'])}")
        L_.append(f"  -> Anteil Produktions-Rauschen am Spread (bis T_req): "
                  f"{d['max_spread_upto_T_req_ratio_vs_baseline']:.2f}")
        L_.append(f"  Budget-Regel: {k['budget_rule']}")
        L_.append(f"  -> guenstigste Variante mit T_req: "
                  f"{k['cheapest_meeting_T_req']}; W4-GO fuer L>=48: "
                  f"{k['w4_go_L_ge_48']}")
        r = k.get("w4_recipe")
        if r:
            L_.append(f"  W4-Rezept: prod x{r['prod_factor']}, lnf_final "
                      f"{r['lnf_final']:g}; Numba-Wall je Walker (Projektion):")
            for Lp, v in r["projection"].items():
                L_.append(f"    L={Lp}: {v['wall_numba_h_per_walker']:.2f} h")
    L_.append("")
    L_.append("--- PASS-Gates (Integritaet, keine Physik) ---")
    for g, v in rep["pass_gates"].items():
        L_.append(f"  [{'PASS' if v else 'FAIL'}] {g}")
    L_.append(f"  OVERALL: {'PASS' if rep['overall_pass'] else 'FAIL'}")
    L_.append("")
    L_.append("Grenzen (ehrlich): Kosten-Projektionen sind Extrapolationen "
              "aus L<=64;")
    L_.append("die Domaenen-Aussage gilt fuer das kalibrierte Budget und diese "
              "Streams;")
    L_.append("kein T_BKT-Wert, keine Band-Aussage (Blind-Vertrag).")
    path.write_text("\n".join(L_) + "\n", encoding="utf-8")


def build_report(stage_dir: Path) -> None:
    stages = {}
    for name in STAGE_NAMES:
        f = Path(stage_dir) / f"{name}.json"
        if f.exists():
            stages[name] = json.loads(f.read_text(encoding="utf-8"))
    rep = _clean(analyse(stages))
    rep["stage_meta"] = {n: {k: stages[n].get(k) for k in (
        "numba_version", "numpy_version", "stage_wall_s",
        "blas_threads_env")} for n in stages}
    out = _ROOT / "results"
    (out / f"{REPORT_STEM}.json").write_text(
        json.dumps(rep, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    write_text_report(rep, out / f"{REPORT_STEM}.txt")
    print(f"Report: results/{REPORT_STEM}.{{json,txt}} "
          f"(OVERALL {'PASS' if rep['overall_pass'] else 'FAIL'})")


STAGES = {"valbit": stage_valbit, "speed": stage_speed,
          "ladder": stage_ladder, "levers": stage_levers}
STAGE_NAMES = tuple(STAGES)


def _main(argv: list) -> None:
    """python "<this file>" <stage> --out DIR [--workers N]
    Stufen: valbit | speed | ladder | levers (je DIR/<stage>.json);
    report DIR  -> results/260926 PHY044 ... report.{json,txt}."""
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=[*STAGES, "report"])
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--workers", type=int,
                    default=max(1, min(4, os.cpu_count() or 1)))
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    if a.stage == "report":
        build_report(a.out)
        return
    fn = STAGES[a.stage]
    kw = {} if a.stage == "speed" else {"max_workers": a.workers}
    t0 = time.perf_counter()
    res = fn(**kw)
    res["numba_version"] = numba.__version__ if HAVE_NUMBA else None
    res["numpy_version"] = np.__version__
    res["blas_threads_env"] = {v: os.environ.get(v) for v in (
        "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")}
    res["stage_wall_s"] = time.perf_counter() - t0
    path = a.out / f"{a.stage}.json"
    path.write_text(json.dumps(_clean(res), indent=1, allow_nan=False),
                    encoding="utf-8")
    print(f"[{a.stage}] geschrieben: {path} ({res['stage_wall_s']:.0f}s)")


if __name__ == "__main__":
    _main(sys.argv[1:])
