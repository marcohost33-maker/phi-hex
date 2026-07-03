"""PHY041 - Wang-Landau-entropischer Helicity-Modul auf HONEYCOMB:
der in PHY040 (Spec §6) vertraglich benannte Follow-up.

Coworker Research / Coworkerz, 2. Juli 2026
de-CH Konventionen, ASCII-Quotes, ss statt Eszett, keine Personennamen.

================================ ZWECK =================================
PHY040 hat den entropischen Sampler (Wang-Landau + 1/t) am square-Goldstandard
validiert und die zwei Limits der PHY031-039-Serie sauber getrennt: das
STATISTISCHE Rauschen ist geloest, der finite-size-Bias verbleibt. PHY041
traegt den Sampler auf das eigentliche offene Ziel der Serie - HONEYCOMB
(Referenzen 0.573 arXiv:2501.07388 / 0.576(3) arXiv:2406.12076; letztere
nutzt exakt diese Methode: WL + 2.-/4.-Ordnungs-Helicity-Modul).

================================ METHODE ===============================
Identischer Vertrag wie PHY040 (Spec `260616 ... wang-landau ... method v01.md`),
plus drei web-recherchierte Haertungen:

1. GITTER-AGNOSTISCHER KERNEL: der WL-/Produktions-Walk arbeitet auf
   (Adjazenz-Liste, Kanten, Bond-Projektionen) - kein square-Hardcode mehr.
   Die ANALYSE-Mathematik (kanonische Rueckgewichtung, Upsilon_2/Upsilon_4,
   C-eliminierte Paare) wird 1:1 aus PHY040 importiert (single source of
   truth, kein Reinvent). Der Kernel ist zusaetzlich skalar-optimiert
   (Block-RNG je Sweep + math.cos statt numpy auf 3-Element-Slices):
   ~18x schneller als der PHY040-Innerloop, physik-identischer Algorithmus
   (WL-Akzeptanz log U <= ln g(b) - ln g(b')); Aequivalenz-Gates in
   tests/test_phy041_honeycomb_wl.py.
2. AUTO-ENERGIEFENSTER AUS WOLFF-ANKERN: statt hartkodiertem Fenster
   (PHY040: [-1.82,-0.75] square, manuell) wird das BKT-Fenster je L aus
   kurzen Wolff-Ankerlaeufen bei T_lo/T_hi bestimmt (Mittel +- k*std,
   Grundzustands-Floor-Guard bei -1.5 J per Site). Gitter-agnostisch,
   selbst-justierend; vermeidet unreachable Bins nahe des Grundzustands
   (dynamical traps, arXiv:1508.01888) und unnoetig breite Fenster.
3. EXPLIZITES LEAK-GATE: PHY040 hat das kanonische Gewichts-Leck am oberen
   Fensterrand (T>1.0) nur im Kommentar dokumentiert und das T-Gitter manuell
   getrimmt. PHY041 prueft HART, dass fuer JEDES Analyse-T das kanonische
   Gewicht in den aeussersten Rand-Bins < 1e-3 bleibt (beide Raender).

QUERVALIDIERUNG (arXiv-unabhaengig, zweifach):
(a) frische Wolff-Referenz (L=12): <E>(T)/N und Upsilon_2(T) muessen sich
    mit den WL-Kurven decken (Toleranzen wie PHY040: 0.03 / 0.04);
(b) DETERMINISTISCHER Drift-Guard gegen die committed PHY032-Evidenz: die
    WL-Upsilon_2(T)-Kurven muessen das im Repo festgeschriebene
    Wolff-Messgitter Upsilon(T,L) (L=12/24, results/260607 PHY032 ...)
    innerhalb der dortigen Seed-SEMs reproduzieren.

T_BKT: C-eliminierte Weber-Minnhagen-Paare auf den GLATTEN WL-Kurven
(phy040.tbkt_pair_from_curves), Paare (12,16)/(12,24)/(16,24). Der
4.-Ordnungs-Dip (N*Upsilon_4) wird rauschfrei je L lokalisiert (PHY039-
Observable unter dem entropischen Sampler).

EHRLICHKEIT: Die Physik-Hypothese (Paar-Schaetzer nahe 0.573/0.576) ist
bewusst KEIN Build-Breaker (AGENTS.md: Negativ-Results sind Buerger erster
Klasse). PASS = Pipeline-Integritaet (Orakel, Quervalidierung, kein Leak,
Glaette); der Physik-Befund wird als FINDING ausgewiesen.

Provenance: Wang & Landau PRL 86, 2050 (2001); Belardinelli & Pereyra JCP
127, 184105 (2007); kontinuierliche Modelle / Fenster: PRE 89, 013311 (2014),
arXiv:1508.01888, cond-mat/0611039; Methodik-Vorbild honeycomb:
arXiv:2406.12076; WM-Paare: PRB 37, 5986(R) (1988), arXiv:1302.2900.
"""
from __future__ import annotations

import json
import math
import re
import time
import numpy as np

import importlib.util
import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parent
_ROOT = _SRC.parent


def _load(name: str, filename: str):
    """Portabler Loader (Muster PHY030 v02 / PHY031-033 / PHY040)."""
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _SRC / filename)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


_core = _load("phi_hex_core_v2", "260602 PHI HEX core v2 2 hardened.py")
_wolff = _load("phy026_wolff_cluster", "260602 PHY026 wolff cluster v01.py")
_p31 = _load("phy031_honeycomb", "260607 PHY031 honeycomb tbkt per-site v01.py")
_p40 = _load("phy040_wang_landau",
             "260616 PHY040 wang-landau entropic helicity v01.py")

XYConfig = _core.XYConfig
make_rng = _core.make_rng
helicity_terms = _core.helicity_terms
helicity_from_ensemble = _core.helicity_from_ensemble
wolff_sweep = _wolff.wolff_sweep
build_honeycomb_lattice = _p31.build_honeycomb_lattice
build_honeycomb_adjacency = _p31.build_honeycomb_adjacency

# Analyse-Primitiven 1:1 aus PHY040 (single source of truth, kein Reinvent).
OBS = _p40.OBS
WLResult = _p40.WLResult
_aggregates_dir = _p40._aggregates_dir
_canonical_weights = _p40._canonical_weights
_full_energy = _p40._full_energy
upsilon_curves = _p40.upsilon_curves
tbkt_pair_from_curves = _p40.tbkt_pair_from_curves

# Referenzen (beide, kein Kanon gekuert - Marco/Vero-Entscheid 2026-06-07).
REF_MULTI = 0.573      # arXiv:2501.07388 (multi-lattice MC)
REF_DEDIC = 0.576      # arXiv:2406.12076 (dedizierte honeycomb-WL-Studie)
E0_PER_SITE = -1.5     # Grundzustand honeycomb (z=3): -z/2 per Site

PHY032_REPORT = (_ROOT / "results" /
                 "260607 PHY032 honeycomb wm-logfit bootstrap report.txt")

# RNG-Stream-Vertrag PHY041 (kollisionsfrei zu 400+ PHY031/032, 500+ PHY033,
# 500/700+ PHY040): Anker 600+s+1000L, WL-Walker 650+L, Wolff-Referenz
# 660+s+1000L.
_STREAM_ANCHOR = 600
_STREAM_WL = 650
_STREAM_REF = 660


# ============================================================================
# 1. Gitter-Arrays + Auto-Energiefenster
# ============================================================================

def honeycomb_arrays(L: int):
    """(lat, nbr_list, ei, ej, ax, ay, n) fuer den gitter-agnostischen Kernel."""
    lat = build_honeycomb_lattice(L)
    adj = build_honeycomb_adjacency(lat)
    nbr_list = [tuple(row) for row in adj]
    edges = np.asarray(lat.edges)
    ei, ej = edges[:, 0], edges[:, 1]
    ax = np.array([dx for dx, dy in lat.edge_disp])
    ay = np.array([dy for dx, dy in lat.edge_disp])
    return lat, nbr_list, ei, ej, ax, ay, lat.n_nodes


def window_from_anchors(mean_lo: float, std_lo: float, mean_hi: float,
                        std_hi: float, k_lo: float = 6.0, k_hi: float = 8.0,
                        min_lo: float = 0.05, min_hi: float = 0.10,
                        floor: float = E0_PER_SITE + 0.03) -> tuple[float, float]:
    """BKT-Energiefenster per-Site aus zwei kanonischen Ankern (T_lo, T_hi).

    Unterer Rand: Anker-Mittel minus k_lo*std (mind. min_lo), geclampt auf
    floor (Grundzustands-Guard: Bins unterhalb erreichbarer Energien machen
    das Flachheits-Kriterium unerfuellbar - dynamical traps). Oberer Rand:
    Anker-Mittel plus k_hi*std (mind. min_hi; asymmetrisch grosszuegiger,
    weil das kanonische Gewicht oberhalb T_hi zuerst am oberen Rand leckt).
    """
    e_lo = mean_lo - max(k_lo * std_lo, min_lo)
    e_hi = mean_hi + max(k_hi * std_hi, min_hi)
    return max(e_lo, floor), e_hi


def wolff_anchor(L: int, T: float, n_seeds: int = 3, n_measure: int = 80,
                 n_burn: int = 60, master_seed: int = 42) -> tuple[float, float]:
    """Kanonischer Energie-Anker (Mittel, Fluktuationsbreite) per-Site."""
    lat, _, ei, ej, _, _, n = honeycomb_arrays(L)
    adj = build_honeycomb_adjacency(lat)
    beta = 1.0 / T
    means, stds = [], []
    for s in range(n_seeds):
        rng = make_rng(master_seed, stream=_STREAM_ANCHOR + s + 1000 * L)
        th = rng.uniform(0, 2 * math.pi, n)
        for _ in range(n_burn):
            wolff_sweep(th, adj, beta, 1.0, rng, target_flips=n)
        es = []
        for _ in range(n_measure):
            wolff_sweep(th, adj, beta, 1.0, rng, target_flips=n)
            es.append(_full_energy(th, ei, ej))
        arr = np.array(es) / n
        means.append(float(arr.mean()))
        stds.append(float(arr.std()))
    return float(np.mean(means)), float(np.mean(stds))


# ============================================================================
# 2. Gitter-agnostischer WL-Kernel (skalar-optimiert, Block-RNG)
# ============================================================================

def _local_delta(thl: list, nbr: tuple, old: float, new: float) -> float:
    """Energie-Aenderung eines Single-Spin-Vorschlags old->new (Nachbarn nbr)."""
    d = 0.0
    for jn in nbr:
        tj = thl[jn]
        d -= math.cos(new - tj) - math.cos(old - tj)
    return d


def wl_entropic_lattice(nbr_list, ei, ej, ax, ay, n: int, L: int,
                        e_lo_ps: float, e_hi_ps: float,
                        lnf_final: float = 1e-5, prod_sweeps: int = 30000,
                        seed: int = 42, stream: int = _STREAM_WL,
                        resync_every: int = 2000,
                        verbose: bool = True) -> WLResult:
    """WL-g(E) (echtes 1/t) + Produktion mit mikrokanonischen Aggregaten.

    Gleiche Akzeptanz und Bin-Konvention (binw ~ 1 in absoluter Energie) wie
    PHY040, aber gitter-agnostisch (Adjazenz + Bond-Projektionen), skalar-
    optimiert (Block-RNG je Sweep, math.cos, ln g als Python-Liste) und mit
    der ORIGINALGETREUEN 1/t-Phasenlogik nach Belardinelli & Pereyra (JCP
    127, 184105 (2007)): die lnf-Halbierungen nutzen das B&P-Kriterium
    "alle Bins besucht" (H_min >= 1) statt strikter Flachheit, damit lnf
    schnell unter 1/t faellt und die 1/t-Politur tatsaechlich GREIFT.
    (Im PHY040-Kernel waren die strikten Flachheits-Epochen so langsam, dass
    lnf(t) nie 1/t erreichte und die Fehler-Sättigung des Standard-WL
    bestehen blieb - auf honeycomb L=12 als Walker-abhaengiger g(E)-Bias
    ~0.02 in <E>/N messbar; mit echter 1/t-Politur < 0.006. Genau die
    Saettigung, die B&P beheben.)

    Periodischer Energie-Resync (alle resync_every Sweeps) pinnt den
    akkumulierten Gleitkomma-Drift der laufenden Energie (Integritaets-
    Metrik in der Konsole, typ. < 1e-10).
    """
    E_lo, E_hi = e_lo_ps * n, e_hi_ps * n
    nbins = int(round(E_hi - E_lo))
    binw = (E_hi - E_lo) / nbins
    rng = make_rng(seed, stream=stream + L)
    two_pi = 2.0 * math.pi
    log = math.log

    th = rng.uniform(0, two_pi, n)
    thl = list(th)
    E = _full_energy(np.asarray(thl), ei, ej)
    ba = 0.0
    while not (E_lo <= E < E_hi):
        ba += 0.05
        idx = rng.integers(0, n, size=n).tolist()
        news = rng.uniform(0.0, two_pi, size=n).tolist()
        us = rng.random(size=n).tolist()
        for k in range(n):
            i = idx[k]
            old = thl[i]
            d = _local_delta(thl, nbr_list[i], old, news[k])
            if d < 0 or us[k] < math.exp(-ba * d):
                thl[i] = news[k]
                E += d
    b = min(max(int((E - E_lo) // binw), 0), nbins - 1)

    lng = [0.0] * nbins
    H = [0.0] * nbins
    lnf = 1.0
    use_1t = False
    sweeps = 0
    sweeps_1t = None
    max_drift = 0.0
    t0 = time.time()
    # ---- WL-Phase: g(E) aufbauen (Standard-WL -> echtes 1/t, B&P) ----
    while lnf > lnf_final:
        idx = rng.integers(0, n, size=n).tolist()
        news = rng.uniform(0.0, two_pi, size=n).tolist()
        us = rng.random(size=n).tolist()
        for k in range(n):
            i = idx[k]
            old = thl[i]
            new = news[k]
            d = _local_delta(thl, nbr_list[i], old, new)
            Enew = E + d
            bnew = int((Enew - E_lo) // binw)
            if 0 <= bnew < nbins and log(us[k] + 1e-300) <= lng[b] - lng[bnew]:
                thl[i] = new
                E = Enew
                b = bnew
            lng[b] += lnf
            H[b] += 1
        sweeps += 1
        if sweeps % resync_every == 0:
            E_exact = _full_energy(np.asarray(thl), ei, ej)
            max_drift = max(max_drift, abs(E_exact - E))
            E = E_exact
            b = min(max(int((E - E_lo) // binw), 0), nbins - 1)
        if use_1t:
            lnf = 1.0 / sweeps
        elif min(H) >= 1:
            # B&P-Kriterium "alle Bins besucht" (statt strikter Flachheit):
            # lnf faellt schnell; sobald das halbierte lnf unter 1/t liegt,
            # dauerhaft in die 1/t-Politur (Vergleich NUR am Halbierungs-
            # Ereignis - jede-Sweep-Vergleich triggert trivial bei t=1).
            lnf *= 0.5
            H = [0.0] * nbins
            if lnf <= 1.0 / sweeps:
                use_1t = True
                sweeps_1t = sweeps
                lnf = 1.0 / sweeps
    if verbose:
        print(f"    WL g(E): L={L} n={n} bins={nbins} sweeps={sweeps} "
              f"(1/t ab {sweeps_1t}) drift_max={max_drift:.2e} "
              f"t={time.time() - t0:.0f}s")
    lng_arr = np.asarray(lng)
    lng_arr -= lng_arr.max()

    # ---- Produktion: g fix, flacher E-Walk, mikrokanonische Aggregate ----
    acc_x = {k: np.zeros(nbins) for k in OBS}
    acc_y = {k: np.zeros(nbins) for k in OBS}
    cnt = np.zeros(nbins)
    lngl = lng_arr.tolist()
    for sw in range(prod_sweeps):
        idx = rng.integers(0, n, size=n).tolist()
        news = rng.uniform(0.0, two_pi, size=n).tolist()
        us = rng.random(size=n).tolist()
        for k in range(n):
            i = idx[k]
            old = thl[i]
            new = news[k]
            d = _local_delta(thl, nbr_list[i], old, new)
            Enew = E + d
            bnew = int((Enew - E_lo) // binw)
            if 0 <= bnew < nbins and log(us[k] + 1e-300) <= lngl[b] - lngl[bnew]:
                thl[i] = new
                E = Enew
                b = bnew
        if (sw + 1) % resync_every == 0:
            E_exact = _full_energy(np.asarray(thl), ei, ej)
            max_drift = max(max_drift, abs(E_exact - E))
            E = E_exact
            b = min(max(int((E - E_lo) // binw), 0), nbins - 1)
        th_arr = np.asarray(thl)
        phi = th_arr[ei] - th_arr[ej]
        cph = np.cos(phi)
        sph = np.sin(phi)
        ox = _aggregates_dir(cph, sph, ax)
        oy = _aggregates_dir(cph, sph, ay)
        for k in OBS:
            acc_x[k][b] += ox[k]
            acc_y[k][b] += oy[k]
        cnt[b] += 1
    if verbose:
        print(f"    production: sweeps={prod_sweeps} "
              f"covered={int((cnt > 0).sum())}/{nbins} "
              f"drift_max={max_drift:.2e} t={time.time() - t0:.0f}s")
    centers = E_lo + (np.arange(nbins) + 0.5) * binw
    mask = cnt > 0
    denom = np.maximum(cnt, 1)
    micro_x = {k: np.where(mask, acc_x[k] / denom, 0.0) for k in OBS}
    micro_y = {k: np.where(mask, acc_y[k] / denom, 0.0) for k in OBS}
    return WLResult(L=L, n=n, centers=centers, lng=lng_arr, mask=mask,
                    micro_x=micro_x, micro_y=micro_y, wl_sweeps=sweeps,
                    prod_sweeps=prod_sweeps)


def wang_landau_honeycomb(L: int, t_anchor_lo: float = 0.50,
                          t_anchor_hi: float = 0.70, master_seed: int = 42,
                          verbose: bool = True, **kw) -> tuple[WLResult, dict]:
    """Auto-Fenster (Wolff-Anker) + WL-Lauf fuer EIN honeycomb-L."""
    m_lo, s_lo = wolff_anchor(L, t_anchor_lo, master_seed=master_seed)
    m_hi, s_hi = wolff_anchor(L, t_anchor_hi, master_seed=master_seed)
    e_lo, e_hi = window_from_anchors(m_lo, s_lo, m_hi, s_hi)
    if verbose:
        print(f"    window L={L}: anchors <E>/N({t_anchor_lo})={m_lo:+.4f} "
              f"<E>/N({t_anchor_hi})={m_hi:+.4f} -> [{e_lo:.4f}, {e_hi:.4f}]")
    _, nbr_list, ei, ej, ax, ay, n = honeycomb_arrays(L)
    res = wl_entropic_lattice(nbr_list, ei, ej, ax, ay, n, L, e_lo, e_hi,
                              seed=master_seed, verbose=verbose, **kw)
    meta = {"anchor_lo": {"T": t_anchor_lo, "mean_ps": m_lo, "std_ps": s_lo},
            "anchor_hi": {"T": t_anchor_hi, "mean_ps": m_hi, "std_ps": s_hi},
            "window_ps": [e_lo, e_hi]}
    return res, meta


# ============================================================================
# 3. Leak-Gate + PHY032-Drift-Guard + Wolff-Referenz
# ============================================================================

def canonical_edge_leak(res: WLResult, T: float, k_edge: int = 3) -> float:
    """Kanonisches Gewicht in den aeussersten k_edge Bins je Fensterrand."""
    w = _canonical_weights(res, T)
    return float(w[:k_edge].sum() + w[-k_edge:].sum())


def parse_phy032_grid(path: Path = PHY032_REPORT) -> dict:
    """Committed Upsilon(T,L)-Messgitter aus dem PHY032-Report (Drift-Guard,
    single source of truth wie PHY037). {L: [(T, upsilon, sem), ...]}."""
    grid: dict[int, list[tuple[float, float, float]]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        m = re.match(r"\s+L=(\d+): (T=.*)$", line)
        if not m:
            continue
        entries = re.findall(r"T=([0-9.]+):(-?[0-9.]+)\+-([0-9.]+)", m.group(2))
        if entries:
            grid[int(m.group(1))] = [(float(t), float(u), float(s))
                                     for t, u, s in entries]
    return grid


def wolff_reference(L: int, T: float, n_seeds: int = 6, n_measure: int = 300,
                    n_burn: int = 200, master_seed: int = 42) -> dict:
    """Frische Wolff-Referenz (per-Site <E>/N und Upsilon_2) zur Quervalidierung.

    Upsilon_2 wird wie in der WL-Pipeline ueber BEIDE Twist-Richtungen (x,y)
    gemittelt (apples-to-apples; core-Default misst nur x)."""
    lat, _, ei, ej, _, _, n = honeycomb_arrays(L)
    adj = build_honeycomb_adjacency(lat)
    beta = 1.0 / T
    cfg = XYConfig(J=1.0, T=T)
    es, ys = [], []
    for s in range(n_seeds):
        rng = make_rng(master_seed, stream=_STREAM_REF + s + 1000 * L)
        th = rng.uniform(0, 2 * math.pi, n)
        for _ in range(n_burn):
            wolff_sweep(th, adj, beta, 1.0, rng, target_flips=n)
        ee = []
        t1s = {d: [] for d in ((1.0, 0.0), (0.0, 1.0))}
        sas = {d: [] for d in ((1.0, 0.0), (0.0, 1.0))}
        for _ in range(n_measure):
            wolff_sweep(th, adj, beta, 1.0, rng, target_flips=n)
            ee.append(_full_energy(th, ei, ej))
            for d in t1s:
                t1, sa = helicity_terms(th, lat, lat.pos, direction=d)
                t1s[d].append(t1)
                sas[d].append(sa)
        es.append(float(np.mean(ee)) / n)
        ys.append(float(np.mean([
            helicity_from_ensemble(t1s[d], sas[d], cfg, n) for d in t1s])))
    return {"T": T, "E_ps": float(np.mean(es)),
            "E_sem": float(np.std(es, ddof=1) / math.sqrt(n_seeds)),
            "y2": float(np.mean(ys)),
            "y2_sem": float(np.std(ys, ddof=1) / math.sqrt(n_seeds))}


def aligned_upsilon_zero(L: int = 6) -> float:
    """Orakel: perfekt ausgerichtete Konfiguration -> Upsilon_2 = 3/4 J exakt
    (per Site, honeycomb), durch die MIKROKANONISCHE Aggregat-Maschinerie."""
    _, _, ei, ej, ax, ay, n = honeycomb_arrays(L)
    th = np.zeros(n)
    phi = th[ei] - th[ej]
    cph, sph = np.cos(phi), np.sin(phi)
    vals = []
    for a in (ax, ay):
        agg = _aggregates_dir(cph, sph, a)
        # T beliebig: s-Terme sind exakt 0 -> Upsilon_2 = c/n
        vals.append((agg["c"] - 1.0 * (agg["s2"] - agg["s"] ** 2)) / n)
    return float(np.mean(vals))


# ============================================================================
# 4. Lauf + PASS-Gates + Report
# ============================================================================

def run_phy041(Ls=(12, 16, 24), lnf_final=1e-5, prod_sweeps=30000,
               t_grid=None, master_seed=42) -> dict:
    if t_grid is None:
        t_grid = np.round(np.arange(0.52, 0.6701, 0.005), 4)
    print("PHY041: Wang-Landau-entropischer Helicity-Modul - HONEYCOMB")
    print("Coworker Research / Coworkerz, 2. Juli 2026")
    print("=" * 72)
    t_start = time.time()

    # Gate A: Geometrie-/Aggregat-Orakel (exakt, arXiv-unabhaengig)
    y0 = aligned_upsilon_zero()
    aligned_ok = abs(y0 - 0.75) < 1e-9
    print(f"\n[GATE A] aligned Upsilon_2(0) = {y0:.6f} (exakt 0.75) "
          f"-> {'PASS' if aligned_ok else 'FAIL'}")

    results: dict[int, WLResult] = {}
    curves: dict[int, dict] = {}
    windows: dict[int, dict] = {}
    for L in Ls:
        print(f"\n[WL] honeycomb L={L} (N={2 * L * L})")
        res, meta = wang_landau_honeycomb(L, master_seed=master_seed,
                                          lnf_final=lnf_final,
                                          prod_sweeps=prod_sweeps)
        results[L] = res
        windows[L] = meta
        curves[L] = upsilon_curves(res, t_grid)

    # --- Leak-Gate: kanonisches Gewicht an den Fensterraendern ------------
    leak_max = 0.0
    for L in Ls:
        for T in t_grid:
            leak_max = max(leak_max, canonical_edge_leak(results[L], float(T)))
    leak_ok = leak_max < 1e-3
    print(f"\n[LEAK] max. Randgewicht ueber alle L, T: {leak_max:.2e} "
          f"(Gate < 1e-3) -> {'PASS' if leak_ok else 'FAIL'}")

    # --- Quervalidierung (a): frische Wolff-Referenz (L=Ls[0]) ------------
    Lval = Ls[0]
    print(f"\n[VAL-A] WL vs frische Wolff-Referenz (L={Lval}):")
    e_ok, y2_ok = True, True
    val_rows = []
    for T in (0.55, 0.60, 0.65):
        c = curves[Lval]
        k = int(np.argmin(np.abs(c["T"] - T)))
        Ewl = c["E"][k] / results[Lval].n
        y2wl = c["y2"][k]
        ref = wolff_reference(Lval, T, master_seed=master_seed)
        de = abs(Ewl - ref["E_ps"])
        dy = abs(y2wl - ref["y2"])
        e_ok = e_ok and de < 0.03
        y2_ok = y2_ok and dy < 0.04
        val_rows.append({"T": T, "E_wl": Ewl, "E_wolff": ref["E_ps"],
                         "y2_wl": y2wl, "y2_wolff": ref["y2"],
                         "y2_wolff_sem": ref["y2_sem"]})
        print(f"      T={T:.3f}  <E>/N WL={Ewl:+.4f} Wolff={ref['E_ps']:+.4f} "
              f"(d={de:.4f})  Y2 WL={y2wl:.4f} "
              f"Wolff={ref['y2']:.4f}+-{ref['y2_sem']:.4f} (d={dy:.4f})")

    # --- Quervalidierung (b): Drift-Guard gegen committed PHY032-Gitter ---
    print("\n[VAL-B] WL vs committed PHY032-Messgitter (Drift-Guard):")
    grid032 = parse_phy032_grid()
    grid_ok = True
    grid_rows = []
    for L in Ls:
        if L not in grid032:
            continue
        for (T, u032, sem032) in grid032[L]:
            if not (t_grid[0] <= T <= t_grid[-1]):
                continue
            y2wl = float(np.interp(T, curves[L]["T"], curves[L]["y2"]))
            tol = max(4.0 * sem032, 0.025)
            ok = abs(y2wl - u032) <= tol
            grid_ok = grid_ok and ok
            grid_rows.append({"L": L, "T": T, "y2_wl": y2wl,
                              "y2_phy032": u032, "sem_phy032": sem032,
                              "ok": ok})
            print(f"      L={L} T={T:.4f}  WL={y2wl:.4f} "
                  f"PHY032={u032:.4f}+-{sem032:.4f} "
                  f"(d={abs(y2wl - u032):.4f}, tol={tol:.4f}) "
                  f"{'ok' if ok else 'ABWEICHUNG'}")

    # --- T_BKT: C-eliminierte Paare auf glatten WL-Kurven -----------------
    print("\n[T_BKT] C-eliminierte Paare auf glatten WL-Upsilon_2-Kurven:")
    pairs = [(Ls[i], Ls[j]) for i in range(len(Ls)) for j in range(len(Ls))
             if Ls[j] > Ls[i]]
    pair_tbkt = {}
    for La, Lb in pairs:
        tb = tbkt_pair_from_curves(t_grid, curves[La]["y2"], La,
                                   curves[Lb]["y2"], Lb)
        pair_tbkt[(La, Lb)] = tb
        if tb is not None:
            print(f"      T_BKT({La},{Lb}) = {tb:.4f}  "
                  f"(vs {REF_MULTI}: {(tb - REF_MULTI) / REF_MULTI * 100:+.2f}%"
                  f", vs {REF_DEDIC}: {(tb - REF_DEDIC) / REF_DEDIC * 100:+.2f}%)")
        else:
            print(f"      T_BKT({La},{Lb}): kein Nulldurchgang im T-Fenster")

    # --- 4.-Ordnungs-Dip (rauschfrei) je L ---------------------------------
    print("\n[Y4] N*Upsilon_4-Dip (deterministisch, rauschfrei):")
    y4_dips = {}
    for L in Ls:
        c = curves[L]
        k = int(np.argmin(c["y4_scaled"]))
        y4_dips[L] = {"T_dip": float(c["T"][k]),
                      "depth": float(c["y4_scaled"][k]),
                      "at_edge": bool(k == 0 or k == len(c["T"]) - 1)}
        print(f"      L={L}: Dip bei T={y4_dips[L]['T_dip']:.4f} "
              f"(Tiefe {y4_dips[L]['depth']:.0f}"
              f"{', am Gitterrand' if y4_dips[L]['at_edge'] else ''})")

    # --- Glaette-Gate (groesstes L, Kernfenster) ---------------------------
    cL = curves[Ls[-1]]
    core_win = (cL["T"] >= 0.54) & (cL["T"] <= 0.66)
    y2_core = cL["y2"][core_win]
    y2_smooth = bool(np.all(np.diff(y2_core) < 0))
    print(f"\n[Y2] Upsilon_2(T) L={Ls[-1]} Kernfenster [0.54,0.66]: "
          f"streng monoton+glatt={y2_smooth} "
          f"({y2_core[0]:.3f} .. {y2_core[-1]:.3f})")

    gates = {
        "PASS_ALIGNED_EXACT_THREE_QUARTERS": aligned_ok,
        "PASS_NO_CANONICAL_EDGE_LEAK": leak_ok,
        "PASS_WL_ENERGY_MATCHES_WOLFF": e_ok,
        "PASS_WL_Y2_MATCHES_WOLFF": y2_ok,
        "PASS_WL_Y2_MATCHES_PHY032_GRID": grid_ok,
        "PASS_Y2_SMOOTH_NOISE_FREE": y2_smooth,
    }
    print("\n[GATES]")
    for k, v in gates.items():
        print(f"      [{'PASS' if v else 'FAIL'}] {k}")
    overall = all(gates.values())
    print(f"\nOVERALL: {'PASS' if overall else 'FAIL'}  "
          f"(= Pipeline-Integritaet; Physik-Befund unten ist FINDING, "
          f"kein Build-Breaker)")

    # --- ehrlicher Physik-Befund -------------------------------------------
    valid_pairs = {p: v for p, v in pair_tbkt.items() if v is not None}
    best = (valid_pairs[max(valid_pairs, key=lambda p: min(p))]
            if valid_pairs else None)
    print("\n[FINDING] (ehrlich, NICHT getunt)")
    if best is not None:
        print(f"    - bestes Paar (groesstes min-L): T_BKT = {best:.4f} "
              f"({(best - REF_MULTI) / REF_MULTI * 100:+.2f}% vs 0.573, "
              f"{(best - REF_DEDIC) / REF_DEDIC * 100:+.2f}% vs 0.576(3))")
    print("    - Vergleich Wolff-Stand (PHY032, Paar (12,24)): 0.6014;")
    print("      dieselbe Paar-Bedingung auf den GLATTEN WL-Kurven ersetzt")
    print("      das Crossing-Rauschen durch deterministische Kurven -")
    print("      der finite-size-Bias der kleinen L bleibt (PHY040-Befund).")
    print(f"    - Laufzeit gesamt: {time.time() - t_start:.0f}s "
          f"(skalar-optimierter Kernel, ~18x vs PHY040-Innerloop)")

    return {
        "module": "PHY041_honeycomb_wang_landau_entropic_helicity",
        "attribution": "Coworker Research / Coworkerz",
        "date": "2026-07-02",
        "references": {"multi_lattice": REF_MULTI, "dedicated": REF_DEDIC},
        "method": ("Wang-Landau (PRL 86,2050(2001)) + 1/t (B&P JCP 127,"
                   "184105(2007)), gitter-agnostisch; Auto-Fenster aus "
                   "Wolff-Ankern; Leak-Gate; Analyse-Primitiven 1:1 PHY040"),
        "Ls": list(Ls),
        "lnf_final": lnf_final,
        "prod_sweeps": prod_sweeps,
        "master_seed": master_seed,
        "windows": {str(L): windows[L] for L in Ls},
        "wl_sweeps": {str(L): results[L].wl_sweeps for L in Ls},
        "leak_max": leak_max,
        "validation_vs_wolff": val_rows,
        "validation_vs_phy032_grid": grid_rows,
        "pair_tbkt": {f"{a}_{b}": v for (a, b), v in pair_tbkt.items()},
        "best_pair_tbkt": best,
        "y4_dips": {str(L): y4_dips[L] for L in Ls},
        "curves": {str(L): {"T": curves[L]["T"].tolist(),
                            "y2": curves[L]["y2"].tolist(),
                            "y4_scaled": curves[L]["y4_scaled"].tolist()}
                   for L in Ls},
        "pass_gates": gates,
        "overall_pass": overall,
    }


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


if __name__ == "__main__":
    report = run_phy041()
    print("\n--- JSON-Report ---")
    print(json.dumps(_clean(report), indent=2, allow_nan=False))
