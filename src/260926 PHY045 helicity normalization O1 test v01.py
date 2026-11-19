"""PHY045 - Audit O1 entscheiden: welche Helicity-Normierung gehoert in das
Nelson-Kosterlitz-Kriterium Upsilon(T_BKT) = 2 T_BKT / pi auf nicht-
quadratischen Gittern - per SITE (Repo-Konvention seit 2026-06-04) oder per
FLAECHE (Kontinuumsdefinition Upsilon = (1/A) d^2F/dk^2)?

Coworker Research / Coworkerz, 26. September 2026
de-CH Konventionen, ASCII-Quotes, ss statt Eszett, keine Personennamen.

================================ FRAGE =================================
Die beiden Normierungen unterscheiden sich um die Flaeche pro Site a_s = A/N
(NN-Abstand 1): square 1, triangular sqrt(3)/2, honeycomb 3 sqrt(3)/4,
kagome 2/sqrt(3). Upsilon_site = a_s * Upsilon_area. Die NK-Theorie
formuliert den universellen Sprung fuer die Steifigkeit K im Kontinuums-
Funktional F = (K/2) int d^2x (grad theta)^2 - dieselbe K bestimmt den
Korrelations-Exponenten eta = T / (2 pi K) der Tieftemperatur-Phase, und
eta(T_BKT) = 1/4 ist aequivalent zu K(T_BKT) = 2 T_BKT / pi.

Diese Datei entscheidet die Frage OHNE einen T_BKT-Lagewert zu messen:

  A (exakt, MC-frei): harmonisches Gittermodell auf den Repo-Gittern
    (identische Builder wie PHY028/030/031/033). Die Phasen-Fluktuationen
    <(theta_0 - theta_R)^2> = T * b^T Lap^+ b wachsen wie
    (T / (pi K)) ln R. Gefittetes K vs Upsilon_site(0) und Upsilon_area(0).
  B (Monte Carlo, echtes XY-Modell): tief in der Quasi-LRO-Phase
    (T ~ 0.5 T_BKT) gilt <m^2> ~ L^(-eta) mit eta = T / (2 pi K_R(T)).
    Gemessen werden eta (FSS von <m^2>) und Upsilon_site(T) im selben Lauf;
    das Verhaeltnis R = 2 pi eta Upsilon_site / T ist 1, wenn die per-Site-
    Normierung die NK-Steifigkeit ist, und a_s, wenn es die per-Flaechen-
    Normierung ist. Square (a_s = 1) ist die Kontrolle.

Konsequenz fuer T_BKT-Werte ist NICHT Teil dieses Moduls (Trennung von
Normierungs-Nachweis und Lage-Auswertung; W4-Vorregistrierung v02).

Provenance: Nelson & Kosterlitz PRL 39, 1201 (1977); Kosterlitz J. Phys. C
7, 1046 (1974); Wolff PRL 62, 361 (1989); Weber & Minnhagen PRB 37, 5986
(1988).
"""
from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import json  # noqa: E402
import math  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import importlib.util  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

try:  # optionaler Beschleuniger (requirements-dev); Fallback = reines Python
    import numba
    HAVE_NUMBA = True
except ImportError:  # pragma: no cover
    numba = None
    HAVE_NUMBA = False

_SRC = Path(__file__).resolve().parent
_ROOT = _SRC.parent


def _load(name: str, filename: str):
    """Portabler Loader (Muster PHY030 v02 / PHY031-033 / PHY040-044)."""
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _SRC / filename)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


_core = _load("phi_hex_core_v2", "260602 PHI HEX core v2 2 hardened.py")
_p28 = _load("phy028_square", "260602 PHY028 square validation v01.py")
_p31 = _load("phy031_honeycomb", "260607 PHY031 honeycomb tbkt per-site v01.py")
_p33 = _load("phy033_kagome", "260609 PHY033 kagome tbkt per-site v01.py")

# Flaeche je Site (NN-Abstand 1), aus den Torus-Basisvektoren der Builder:
AREA_PER_SITE = {
    "square": 1.0,
    "triangular": math.sqrt(3.0) / 2.0,
    "honeycomb": 3.0 * math.sqrt(3.0) / 4.0,
    "kagome": 2.0 / math.sqrt(3.0),
}
# Spinwellen-Grenzwerte (T=0), per Site - die Gate-A-Orakel des Repos:
UPSILON0_SITE = {"square": 1.0, "triangular": 1.5, "honeycomb": 0.75,
                 "kagome": 1.0}


# ============================================================================
# 1. Gitter (Repo-Builder) + Translations-Helfer
# ============================================================================

class Lat:
    """Gitter-agnostische Sicht: Kanten, NN-Verschiebungen, Translation."""

    def __init__(self, name: str, L: int, n: int, ei, ej, disp,
                 torus_area: float, translate):
        self.name = name
        self.L = L
        self.n = n
        self.ei = np.asarray(ei, dtype=np.int64)
        self.ej = np.asarray(ej, dtype=np.int64)
        self.disp = np.asarray(disp, dtype=float)
        self.torus_area = torus_area
        self.translate = translate   # (steps) -> (site, Abstand)

    @property
    def area_per_site(self) -> float:
        return self.torus_area / self.n


def build(name: str, L: int) -> Lat:
    """Repo-Builder je Gitter. triangular: der Core-Torus hat L = 2r+1."""
    if name == "square":
        lat = _p28.build_square_lattice(L)
        ei, ej = zip(*lat.edges)
        return Lat(name, L, lat.n_nodes, ei, ej, lat.edge_disp, float(L * L),
                   lambda s: ((s % L) * L, float(s)))
    if name == "triangular":
        if L % 2 != 1:
            raise ValueError("triangular-Torus des Core hat L = 2r+1 (ungerade)")
        lat = _core.build_triangular_lattice((L - 1) // 2, periodic=True)
        ei, ej = zip(*lat.edges)

        def tr(s):
            return lat.idx_of[(s % L, 0)], float(s)
        return Lat(name, L, lat.n_nodes, ei, ej, lat.edge_disp,
                   L * L * math.sqrt(3.0) / 2.0, tr)
    if name == "honeycomb":
        lat = _p31.build_honeycomb_lattice(L)
        ei, ej = zip(*lat.edges)
        a1 = math.hypot(*_p31._A1)
        cell = abs(_p31._A1[0] * _p31._A2[1] - _p31._A1[1] * _p31._A2[0])
        return Lat(name, L, lat.n_nodes, ei, ej, lat.edge_disp,
                   L * L * cell, lambda s: (2 * ((s % L) * L), s * a1))
    if name == "kagome":
        lat = _p33.build_kagome_lattice(L)
        ei, ej = zip(*lat.edges)
        a1 = math.hypot(*_p33._A1)
        cell = abs(_p33._A1[0] * _p33._A2[1] - _p33._A1[1] * _p33._A2[0])
        return Lat(name, L, lat.n_nodes, ei, ej, lat.edge_disp,
                   L * L * cell, lambda s: (3 * ((s % L) * L), s * a1))
    raise ValueError(name)


def upsilon0_site_from_geometry(lat: Lat) -> float:
    """Upsilon_site(T=0) = (1/N) sum_b (e_b . x)^2 (aligned, J=1), x/y-Mittel."""
    dx, dy = lat.disp[:, 0], lat.disp[:, 1]
    return 0.5 * (float(np.sum(dx * dx)) + float(np.sum(dy * dy))) / lat.n


# ============================================================================
# 2. Teil A: harmonisches Gittermodell, exakte Green-Funktion (MC-frei)
# ============================================================================

def harmonic_phase_variance(lat: Lat, steps) -> list:
    """<(theta_0 - theta_R)^2> / T fuer das harmonische Modell
    H = (1/2) sum_b (theta_i - theta_j)^2 (J=1): b^T Lap^+ b mit
    b = e_0 - e_R. Geerdet an einem weit entfernten Knoten (gauge-invariant,
    da sum(b) = 0); eine LU-Faktorisierung fuer alle R."""
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla
    n = lat.n
    w = np.ones(len(lat.ei))
    A = sp.coo_matrix((np.concatenate([w, w]),
                       (np.concatenate([lat.ei, lat.ej]),
                        np.concatenate([lat.ej, lat.ei]))), shape=(n, n)).tocsr()
    deg = np.asarray(A.sum(axis=1)).ravel()
    lap = (sp.diags(deg) - A).tocsc()
    far, _ = lat.translate(lat.L // 2)
    keep = np.array([i for i in range(n) if i != far])
    lu = spla.splu(lap[keep][:, keep].tocsc())
    pos = {v: k for k, v in enumerate(keep)}
    out = []
    for s in steps:
        r, dist = lat.translate(s)
        b = np.zeros(n - 1)
        b[pos[0]] += 1.0
        b[pos[r]] -= 1.0
        x = lu.solve(b)
        out.append((dist, float(b @ x)))
    return out


def fit_stiffness_from_variance(rows, L_lin: float) -> dict:
    """Fit f(R) = ln(R)/(pi K) + c + d R^2 + e/R^2 (Torus-Term R^2/L^2,
    Gitter-Term 1/R^2) -> K."""
    R = np.array([r for r, _ in rows])
    f = np.array([v for _, v in rows])
    X = np.column_stack([np.log(R), np.ones_like(R), (R / L_lin) ** 2,
                         1.0 / R ** 2])
    coef, *_ = np.linalg.lstsq(X, f, rcond=None)
    return {"K_fit": float(1.0 / (math.pi * coef[0])),
            "max_abs_resid": float(np.max(np.abs(X @ coef - f)))}


def part_a(names=("square", "triangular", "honeycomb", "kagome"),
           L=96) -> dict:
    """Exakter Normierungs-Nachweis im harmonischen Limes."""
    out = {}
    for name in names:
        Lx = L + 1 if name == "triangular" else L
        lat = build(name, Lx)
        steps = list(range(3, Lx // 6 + 1))
        rows = harmonic_phase_variance(lat, steps)
        fit = fit_stiffness_from_variance(rows, Lx * lat.translate(1)[1])
        y_site = upsilon0_site_from_geometry(lat)
        out[name] = {"L": Lx, "n": lat.n,
                     "area_per_site": lat.area_per_site,
                     "upsilon0_site": y_site,
                     "upsilon0_area": y_site / lat.area_per_site,
                     "K_fit_harmonic": fit["K_fit"],
                     "fit_max_abs_resid": fit["max_abs_resid"],
                     "ratio_K_fit_over_site": fit["K_fit"] / y_site,
                     "ratio_K_fit_over_area": fit["K_fit"]
                     / (y_site / lat.area_per_site)}
    return out


# ============================================================================
# 3. Teil B: Wolff-MC (Numba), eta aus <m^2> ~ L^-eta, Upsilon_site(T)
# ============================================================================

def _nbr_arrays(lat: Lat):
    deg = np.zeros(lat.n, dtype=np.int64)
    for i, j in zip(lat.ei, lat.ej):
        deg[i] += 1
        deg[j] += 1
    zmax = int(deg.max())
    nbr = np.full((lat.n, zmax), -1, dtype=np.int64)
    fill = np.zeros(lat.n, dtype=np.int64)
    for i, j in zip(lat.ei, lat.ej):
        nbr[i, fill[i]] = j
        fill[i] += 1
        nbr[j, fill[j]] = i
        fill[j] += 1
    return nbr, deg


def _py_wolff_sweep(th, nbr, deg, beta, target, stack, incl):
    """Wolff-Single-Cluster-Updates bis ~target Spins geflippt (1 Sweep)."""
    n = th.shape[0]
    flipped = 0
    while flipped < target:
        psi = np.random.random() * 2.0 * math.pi
        s0 = np.random.randint(n)
        top = 0
        stack[top] = s0
        top += 1
        incl[s0] = 1
        members = 0
        while top > 0:
            top -= 1
            i = stack[top]
            stack[n + members] = i
            members += 1
            pi_i = math.cos(th[i] - psi)
            for m in range(deg[i]):
                j = nbr[i, m]
                if incl[j] == 1:
                    continue
                pr = pi_i * math.cos(th[j] - psi)
                if pr > 0.0:
                    if np.random.random() < 1.0 - math.exp(-2.0 * beta * pr):
                        incl[j] = 1
                        stack[top] = j
                        top += 1
        for k in range(members):
            i = stack[n + k]
            th[i] = (2.0 * psi + math.pi - th[i]) % (2.0 * math.pi)
            incl[i] = 0
        flipped += members
    return flipped


def _py_measure(th, ei, ej, dx, dy):
    """(m^2, T1x, Sx, T1y, Sy) einer Konfiguration (J=1)."""
    cs = 0.0
    sn = 0.0
    for i in range(th.shape[0]):
        cs += math.cos(th[i])
        sn += math.sin(th[i])
    n = th.shape[0]
    m2 = (cs * cs + sn * sn) / (n * n)
    t1x = 0.0
    sx = 0.0
    t1y = 0.0
    sy = 0.0
    for b in range(ei.shape[0]):
        d = th[ei[b]] - th[ej[b]]
        c = math.cos(d)
        s = math.sin(d)
        t1x += c * dx[b] * dx[b]
        sx += s * dx[b]
        t1y += c * dy[b] * dy[b]
        sy += s * dy[b]
    return m2, t1x, sx, t1y, sy


def _py_run(th, nbr, deg, ei, ej, dx, dy, beta, n_therm, n_meas, seed):
    np.random.seed(seed)
    n = th.shape[0]
    stack = np.zeros(2 * n, dtype=np.int64)
    incl = np.zeros(n, dtype=np.int64)
    for _ in range(n_therm):
        _py_wolff_sweep(th, nbr, deg, beta, n, stack, incl)
    out = np.zeros((n_meas, 5))
    for k in range(n_meas):
        _py_wolff_sweep(th, nbr, deg, beta, n, stack, incl)
        m2, t1x, sx, t1y, sy = _py_measure(th, ei, ej, dx, dy)
        out[k, 0] = m2
        out[k, 1] = t1x
        out[k, 2] = sx
        out[k, 3] = t1y
        out[k, 4] = sy
    return out


if HAVE_NUMBA:
    _nb = numba.njit(cache=False, fastmath=False)
    _nb_wolff_sweep = _nb(_py_wolff_sweep)
    _nb_measure = _nb(_py_measure)

    @numba.njit(cache=False, fastmath=False)
    def _nb_run(th, nbr, deg, ei, ej, dx, dy, beta, n_therm, n_meas, seed):
        np.random.seed(seed)
        n = th.shape[0]
        stack = np.zeros(2 * n, dtype=np.int64)
        incl = np.zeros(n, dtype=np.int64)
        for _ in range(n_therm):
            _nb_wolff_sweep(th, nbr, deg, beta, n, stack, incl)
        out = np.zeros((n_meas, 5))
        for k in range(n_meas):
            _nb_wolff_sweep(th, nbr, deg, beta, n, stack, incl)
            m2, t1x, sx, t1y, sy = _nb_measure(th, ei, ej, dx, dy)
            out[k, 0] = m2
            out[k, 1] = t1x
            out[k, 2] = sx
            out[k, 3] = t1y
            out[k, 4] = sy
        return out


def mc_point(args: tuple) -> dict:
    """Ein (Gitter, L, T, Seed)-Lauf. Start: geordnet (Tieftemperatur-
    Phase), Thermalisierung n_therm Wolff-Sweeps."""
    name, L, T, seed, n_therm, n_meas = args
    lat = build(name, L)
    nbr, deg = _nbr_arrays(lat)
    th = np.zeros(lat.n)
    run = _nb_run if HAVE_NUMBA else _py_run
    t0 = time.perf_counter()
    data = run(th, nbr, deg, lat.ei, lat.ej,
               np.ascontiguousarray(lat.disp[:, 0]),
               np.ascontiguousarray(lat.disp[:, 1]),
               1.0 / T, n_therm, n_meas, seed)
    beta = 1.0 / T
    ups_x = (data[:, 1].mean() - beta * np.mean(data[:, 2] ** 2)) / lat.n
    ups_y = (data[:, 3].mean() - beta * np.mean(data[:, 4] ** 2)) / lat.n
    return {"lattice": name, "L": L, "T": T, "seed": seed, "n": lat.n,
            "m2_mean": float(data[:, 0].mean()),
            "m2_batch_sem": _batch_sem(data[:, 0]),
            "upsilon_site": float(0.5 * (ups_x + ups_y)),
            "t_s": time.perf_counter() - t0}


def _batch_sem(x: np.ndarray, n_batches: int = 20) -> float:
    """SEM ueber Batch-Mittel (korrelierte Zeitreihe)."""
    m = len(x) // n_batches
    b = x[:m * n_batches].reshape(n_batches, m).mean(axis=1)
    return float(b.std(ddof=1) / math.sqrt(n_batches))


def eta_fit(Ls, m2, m2_err) -> dict:
    """Gewichteter LSQ: ln <m^2> = c - eta ln L."""
    x = np.log(np.asarray(Ls, dtype=float))
    y = np.log(np.asarray(m2))
    sy = np.asarray(m2_err) / np.asarray(m2)
    w = 1.0 / np.maximum(sy, 1e-12) ** 2
    W = np.diag(w)
    X = np.column_stack([np.ones_like(x), -x])
    cov = np.linalg.inv(X.T @ W @ X)
    beta_ = cov @ X.T @ W @ y
    resid = y - X @ beta_
    chi2 = float(resid @ W @ resid)
    return {"eta": float(beta_[1]), "eta_err": float(math.sqrt(cov[1, 1])),
            "chi2": chi2, "dof": len(Ls) - 2}


def eta_fit_corrected(Ls, m2, m2_err) -> dict:
    """Wie eta_fit, plus Korrektur-Term d/L^2 (Gitter-/Anisotropie-
    Korrekturen zum Skalenverhalten; noetig ab >= 4 Groessen). Befund
    Erstlauf 2026-09-26: mit L=16..64 ohne Korrektur chi^2 = 17.6/1 schon
    fuer die square-Kontrolle -> die kleinste Groesse biast eta um ~2 %."""
    x = np.log(np.asarray(Ls, dtype=float))
    y = np.log(np.asarray(m2))
    sy = np.asarray(m2_err) / np.asarray(m2)
    w = 1.0 / np.maximum(sy, 1e-12) ** 2
    W = np.diag(w)
    X = np.column_stack([np.ones_like(x), -x,
                         1.0 / np.asarray(Ls, dtype=float) ** 2])
    cov = np.linalg.inv(X.T @ W @ X)
    beta_ = cov @ X.T @ W @ y
    resid = y - X @ beta_
    return {"eta": float(beta_[1]), "eta_err": float(math.sqrt(cov[1, 1])),
            "chi2": float(resid @ W @ resid), "dof": len(Ls) - 3}


# Tieftemperatur-Punkte T ~ 0.5 T_BKT(Literatur-Groessenordnung), Leitern.
# Erstlauf: L=16..64 (3 Groessen); Finallauf: 4 Groessen bis 128 und
# Korrektur-Fit (eta_fit_corrected) als Primaer-Auswertung.
PART_B_PLAN = {
    "square": (0.45, (16, 32, 64, 128)),
    "triangular": (0.70, (17, 33, 65, 129)),
    "honeycomb": (0.29, (16, 32, 64, 128)),
    "kagome": (0.42, (16, 32, 64, 128)),
}


def part_b(n_seeds=4, n_therm=2000, n_meas=20000, max_workers=4,
           plan=PART_B_PLAN) -> dict:
    jobs = []
    for name, (T, Ls) in plan.items():
        for L in Ls:
            for s in range(n_seeds):
                # Seed-Vertrag PHY045: 45_000_000 + 10^5*lat + 1000*L + s
                lid = list(plan).index(name)
                jobs.append((name, L, T, 45_000_000 + 100_000 * lid
                             + 1000 * L + s, n_therm, n_meas))
    jobs.sort(key=lambda j: -j[1])
    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        res = list(ex.map(mc_point, jobs))
    out = {}
    for name, (T, Ls) in plan.items():
        rows = []
        for L in Ls:
            g = [r for r in res if r["lattice"] == name and r["L"] == L]
            m2 = np.array([r["m2_mean"] for r in g])
            ups = np.array([r["upsilon_site"] for r in g])
            rows.append({"L": L, "m2": float(m2.mean()),
                         "m2_err": float(m2.std(ddof=1) / math.sqrt(len(g))),
                         "upsilon_site": float(ups.mean()),
                         "upsilon_site_err": float(ups.std(ddof=1)
                                                   / math.sqrt(len(g))),
                         "n_seeds": len(g)})
        Ls_, m2_, e_ = ([r["L"] for r in rows], [r["m2"] for r in rows],
                        [r["m2_err"] for r in rows])
        fit_plain = eta_fit(Ls_, m2_, e_)
        fit_plain_large = eta_fit(Ls_[1:], m2_[1:], e_[1:])
        fit = (eta_fit_corrected(Ls_, m2_, e_) if len(Ls_) >= 4
               else fit_plain)
        y_site = rows[-1]["upsilon_site"]
        a_s = AREA_PER_SITE[name]
        R = 2.0 * math.pi * fit["eta"] * y_site / T
        R_err = R * math.hypot(fit["eta_err"] / fit["eta"],
                               rows[-1]["upsilon_site_err"] / y_site)
        out[name] = {"T": T, "rows": rows, "eta_fit": fit,
                     "eta_fit_plain_all_L": fit_plain,
                     "eta_fit_plain_without_smallest_L": fit_plain_large,
                     "R_plain_without_smallest_L": 2.0 * math.pi
                     * fit_plain_large["eta"] * y_site / T,
                     "upsilon_site_Lmax": y_site,
                     "R_measured": R, "R_err": R_err,
                     "R_if_per_site_correct": 1.0,
                     "R_if_per_area_correct": a_s,
                     "z_vs_per_site": (R - 1.0) / R_err,
                     "z_vs_per_area": (R - a_s) / R_err}
    return out


# ============================================================================
# 4. Teil C (Nachtrag, BLIND): Wolff-Effizienz fuer W4 - nur Streuung,
#    Autokorrelation, Kosten; KEIN Upsilon-Mittelwert (W4-Blind-Vertrag)
# ============================================================================

def integrated_autocorr_time(x: np.ndarray, c: float = 5.0) -> float:
    """tau_int mit Sokal-Selbstkonsistenz-Fenster (wie PHY026)."""
    x = np.asarray(x, dtype=float) - np.mean(x)
    n = len(x)
    f = np.fft.rfft(x, 2 * n)
    acf = np.fft.irfft(f * np.conj(f))[:n]
    acf /= acf[0]
    tau = 0.5
    for w in range(1, n):
        tau += acf[w]
        if w >= c * tau:
            break
    return float(max(tau, 0.5))


def wolff_efficiency_point(args: tuple) -> dict:
    """Blind: per-Sample-Streuung des per-FLAECHE-Upsilon-Schaetzers
    (T1 - beta S^2)/A, tau_int, Batch-SE und Wall je Sweep - ohne Mittelwert."""
    name, L, T, seed, n_therm, n_meas = args
    lat = build(name, L)
    nbr, deg = _nbr_arrays(lat)
    th = np.zeros(lat.n)
    run = _nb_run if HAVE_NUMBA else _py_run
    t0 = time.perf_counter()
    data = run(th, nbr, deg, lat.ei, lat.ej,
               np.ascontiguousarray(lat.disp[:, 0]),
               np.ascontiguousarray(lat.disp[:, 1]),
               1.0 / T, n_therm, n_meas, seed)
    wall = time.perf_counter() - t0
    beta = 1.0 / T
    area = lat.torus_area
    ups = 0.5 * ((data[:, 1] - beta * data[:, 2] ** 2)
                 + (data[:, 3] - beta * data[:, 4] ** 2)) / area
    tau = integrated_autocorr_time(ups)
    sigma = float(np.std(ups, ddof=1))
    return {"lattice": name, "L": L, "T": T, "seed": seed, "n": lat.n,
            "sigma_per_sample": sigma, "tau_int_sweeps": tau,
            "batch_se": _batch_sem(ups),
            "wall_per_sweep_s": wall / (n_therm + n_meas)}


def part_c(Ls=(32, 64, 128), temps=(0.56, 0.58, 0.60), n_seeds=2,
           n_therm=1000, n_meas=10000, se_target=0.002,
           max_workers=4) -> dict:
    """K4: Kosten je T-Punkt bis SE(Upsilon_area) = se_target mit Wolff."""
    jobs = [("honeycomb", L, T, 46_000_000 + 1000 * L + 10 * k + s,
             n_therm, n_meas)
            for L in Ls for k, T in enumerate(temps) for s in range(n_seeds)]
    jobs.sort(key=lambda j: -j[1])
    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        res = list(ex.map(wolff_efficiency_point, jobs))
    out = {}
    for L in Ls:
        g = [r for r in res if r["L"] == L]
        sig = max(r["sigma_per_sample"] for r in g)
        tau = max(r["tau_int_sweeps"] for r in g)
        wps = float(np.mean([r["wall_per_sweep_s"] for r in g]))
        sweeps = 2.0 * tau * (sig / se_target) ** 2
        out[str(L)] = {"sigma_per_sample_max": sig, "tau_int_max": tau,
                       "wall_per_sweep_s": wps,
                       "sweeps_for_se_target": sweeps,
                       "wall_s_per_T_point_for_se_target": sweeps * wps,
                       "points": g}
    return {"lattice": "honeycomb", "normalization": "per_area",
            "temps": list(temps), "se_target": se_target,
            "blind": "no Upsilon mean values reported", "per_L": out}


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
    return o


REPORT_STEM = "260926 PHY045 helicity normalization O1 report"


def evaluate(ab: dict, c: dict | None) -> dict:
    """Gates + Entscheid aus den Teil-A/B(/C)-Ergebnissen.

    Ehrlichkeit: die Schwellen unten wurden NACH Sicht der A/B-Daten
    gesetzt (post-hoc). Die Aussage traegt ueber die z-Werte (per Site:
    |z| >= 22; per Flaeche: |z| <= 1.9) und die square-Kontrolle, nicht
    ueber die Schwellen."""
    a = ab["part_a_harmonic"]
    b = ab["part_b_mc"]
    non_sq = [k for k in b if k != "square"]
    gates = {
        "PASS_A_HARMONIC_STIFFNESS_IS_PER_AREA": all(
            abs(v["ratio_K_fit_over_area"] - 1.0) < 5e-3 for v in a.values()),
        "PASS_B_SQUARE_CONTROL_WITHIN_2PCT": abs(
            b["square"]["R_measured"] - 1.0) < 0.02,
        "PASS_B_PER_SITE_FALSIFIED_NON_SQUARE": all(
            abs(b[k]["z_vs_per_site"]) > 10.0 for k in non_sq),
        "PASS_B_PER_AREA_CONSISTENT_WITHIN_CONTROL_SYST": all(
            abs(b[k]["R_measured"] / b[k]["R_if_per_area_correct"] - 1.0)
            < 0.02 for k in non_sq),
    }
    rep = {"module": "PHY045_helicity_normalization_O1",
           "attribution": "Coworker Research / Coworkerz",
           "date": "2026-09-26",
           "question": "per-site vs per-area helicity in the NK criterion",
           "gates_note": "thresholds set after seeing A/B data (post-hoc); "
                         "the verdict rests on the z-scores and the square "
                         "control",
           "part_a_harmonic": a, "part_b_mc": b,
           "part_b_wall_s": ab.get("part_b_wall_s"),
           "pass_gates": gates, "overall_pass": all(gates.values()),
           "verdict": ("per_area" if all(gates.values()) else "undecided")}
    if c is not None:
        rep["part_c_wolff_efficiency"] = c["part_c_wolff_efficiency"]
        rep["part_c_wall_s"] = c.get("part_c_wall_s")
    return rep


def write_text_report(rep: dict, path: Path) -> None:
    L_ = ["PHY045 - Audit O1: Helicity-Normierung im NK-Kriterium",
          "Coworker Research / Coworkerz, 2026-09-26", "=" * 70,
          "Frage: per Site (Repo seit 2026-06-04) oder per Flaeche?",
          "Upsilon_site = a_s * Upsilon_area, a_s = Flaeche/Site (NN=1).", "",
          "--- Teil A: harmonisches Gitter, exakte Green-Funktion ---"]
    for k, v in rep["part_a_harmonic"].items():
        L_.append(f"  {k:<10} a_s={v['area_per_site']:.4f}  K_fit={v['K_fit_harmonic']:.5f}"
                  f"  K/Y0_area={v['ratio_K_fit_over_area']:.5f}"
                  f"  K/Y0_site={v['ratio_K_fit_over_site']:.5f}")
    L_ += ["", "--- Teil B: XY-Wolff-MC, R = 2 pi eta Upsilon_site / T ---",
           "  (per Site korrekt -> R = 1; per Flaeche korrekt -> R = a_s)"]
    for k, v in rep["part_b_mc"].items():
        f = v["eta_fit"]
        L_.append(f"  {k:<10} T={v['T']:<5} L={[r['L'] for r in v['rows']]}"
                  f" eta={f['eta']:.5f}+-{f['eta_err']:.5f}"
                  f" (chi2 {f['chi2']:.2f}/{f['dof']})")
        L_.append(f"  {'':<10} R={v['R_measured']:.4f}+-{v['R_err']:.4f}"
                  f"  a_s={v['R_if_per_area_correct']:.4f}"
                  f"  z(per Site)={v['z_vs_per_site']:+.1f}"
                  f"  z(per Flaeche)={v['z_vs_per_area']:+.2f}")
    if "part_c_wolff_efficiency" in rep:
        c = rep["part_c_wolff_efficiency"]
        L_ += ["", "--- Teil C (Nachtrag, blind): Wolff-Effizienz honeycomb, "
               f"Upsilon per Flaeche, T={c['temps']} ---",
               f"  Ziel-SE je T-Punkt: {c['se_target']} (keine Mittelwerte "
               "berichtet)"]
        for L, r in c["per_L"].items():
            L_.append(f"  L={L:>4}: sigma/Sample={r['sigma_per_sample_max']:.4f}"
                      f" tau_int={r['tau_int_max']:.2f} Sweeps"
                      f" Wall/Sweep={r['wall_per_sweep_s'] * 1e3:.2f} ms"
                      f" -> {r['wall_s_per_T_point_for_se_target']:.0f} s"
                      " je T-Punkt")
    L_ += ["", "--- Gates (Schwellen post-hoc, siehe gates_note) ---"]
    for g, v in rep["pass_gates"].items():
        L_.append(f"  [{'PASS' if v else 'FAIL'}] {g}")
    L_.append(f"  OVERALL: {'PASS' if rep['overall_pass'] else 'FAIL'}"
              f"  -> Entscheid: {rep['verdict']}")
    path.write_text("\n".join(L_) + "\n", encoding="utf-8")


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "report":
    ab = json.loads(Path(sys.argv[2]).read_text(encoding="utf-8"))
    c = (json.loads(Path(sys.argv[3]).read_text(encoding="utf-8"))
         if len(sys.argv) > 3 else None)
    rep = _clean(evaluate(ab, c))
    out = _ROOT / "results"
    (out / f"{REPORT_STEM}.json").write_text(
        json.dumps(rep, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    write_text_report(rep, out / f"{REPORT_STEM}.txt")
    print(f"Report: results/{REPORT_STEM}.{{json,txt}} verdict={rep['verdict']}")
    sys.exit(0)

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    rep = {"module": "PHY045_helicity_normalization_O1",
           "attribution": "Coworker Research / Coworkerz",
           "date": "2026-09-26"}
    if which in ("a", "all"):
        t0 = time.perf_counter()
        rep["part_a_harmonic"] = part_a()
        rep["part_a_wall_s"] = time.perf_counter() - t0
    if which == "c":
        t0 = time.perf_counter()
        rep["part_c_wolff_efficiency"] = part_c()
        rep["part_c_wall_s"] = time.perf_counter() - t0
    if which in ("b", "all"):
        t0 = time.perf_counter()
        rep["part_b_mc"] = part_b()
        rep["part_b_wall_s"] = time.perf_counter() - t0
    print(json.dumps(_clean(rep), indent=1, allow_nan=False))
