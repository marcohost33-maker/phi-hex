"""Schnelle Korrektheits-Gates fuer PHY044 (W4-Kalibrierung, Numba-Kernel).

(a) VAL-BIT-S: der PHY044-Treiber ist ein BIT-IDENTISCHER Zwilling des
    PHY041-Kernels - im reinen Python-Fallback (laeuft immer) und
    Numba-kompiliert (laeuft, wenn numba installiert ist; requirements-dev).
(b) Treiber-Vertraege: Backend-Wahl fail-closed, Nachbar-Reihenfolge,
    unabhaengige Produktion mutiert den WL-Zustand nicht, Stream-Vertrag.
(c) Blind-Vertrag und Auswertungs-Primitiven (Potenzgesetz, Paar-Streuung).
(d) W4-Vorregistrierung: Konstanten zeichengenau an die Spec gebunden,
    `w4_verdict` gegen Orakel jeder Regel.
Die teuren Laeufe (VAL-BIT-P, K1-K3) sind Gate-Log in results/.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

import conftest

phy044 = conftest._load(
    "phy044_wl_calibration", "260926 PHY044 honeycomb wl calibration v01.py")
p41 = phy044._p41

ROOT = Path(__file__).resolve().parents[1]
SPEC = ROOT / "spec" / "260926 PHI HEX w4 honeycomb preregistration v01.md"


def _small_job(L: int):
    """Kleines honeycomb-Gitter mit festem, erreichbarem Fenster (MC-frei
    bestimmt: zwischen Grundzustand -1.5 und ungeordnet ~0)."""
    _, nbr, ei, ej, ax, ay, n = p41.honeycomb_arrays(L)
    return nbr, ei, ej, ax, ay, n, -1.30, -0.95


# ---------------------------------------------------------------------------
# (a) VAL-BIT-S
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("L,lnf,prod", [(3, 1e-2, 150), (4, 1e-2, 200)])
def test_python_fallback_is_bit_identical_to_phy041(L, lnf, prod):
    nbr, ei, ej, ax, ay, n, lo, hi = _small_job(L)
    ref = p41.wl_entropic_lattice(nbr, ei, ej, ax, ay, n, L, lo, hi,
                                  lnf_final=lnf, prod_sweeps=prod,
                                  verbose=False)
    fast = phy044.wl_entropic_fast(nbr, ei, ej, ax, ay, n, L, lo, hi,
                                   lnf_final=lnf, prod_sweeps=prod,
                                   backend="python")
    assert phy044.results_identical(ref, fast)
    assert ref.wl_sweeps == fast.wl_sweeps > 0


@pytest.mark.parametrize("L,lnf,prod,stream",
                         [(4, 1e-2, 200, 650), (6, 1e-2, 300, 91000)])
def test_numba_kernel_is_bit_identical_to_phy041(L, lnf, prod, stream):
    pytest.importorskip("numba")
    nbr, ei, ej, ax, ay, n, lo, hi = _small_job(L)
    ref = p41.wl_entropic_lattice(nbr, ei, ej, ax, ay, n, L, lo, hi,
                                  lnf_final=lnf, prod_sweeps=prod,
                                  stream=stream, verbose=False)
    fast = phy044.wl_entropic_fast(nbr, ei, ej, ax, ay, n, L, lo, hi,
                                   lnf_final=lnf, prod_sweeps=prod,
                                   stream=stream, backend="numba")
    assert phy044.results_identical(ref, fast)


def test_results_identical_detects_single_bit_flip():
    """Negativ-Kontrolle: der Vergleich ist wirklich bit-genau."""
    nbr, ei, ej, ax, ay, n, lo, hi = _small_job(3)
    a = phy044.wl_entropic_fast(nbr, ei, ej, ax, ay, n, 3, lo, hi,
                                lnf_final=1e-2, prod_sweeps=100,
                                backend="python")
    b = phy044.wl_entropic_fast(nbr, ei, ej, ax, ay, n, 3, lo, hi,
                                lnf_final=1e-2, prod_sweeps=100,
                                backend="python")
    assert phy044.results_identical(a, b)
    k = int(np.argmax(b.mask))
    b.micro_x["c"][k] = np.nextafter(b.micro_x["c"][k], np.inf)
    assert not phy044.results_identical(a, b)


# ---------------------------------------------------------------------------
# (b) Treiber-Vertraege
# ---------------------------------------------------------------------------

def test_backend_selection_fails_closed(monkeypatch):
    with pytest.raises(ValueError):
        phy044._sweep_fns("gpu")
    monkeypatch.setattr(phy044, "HAVE_NUMBA", False)
    with pytest.raises(RuntimeError):
        phy044._sweep_fns("numba")
    assert phy044._sweep_fns("auto")[1] is phy044._py_wl_sweep


def test_nbr_array_preserves_tuple_order_and_degree():
    nbr, deg = phy044._nbr_array([(2, 1), (0,), (1, 0, 3), (2,)])
    assert deg.tolist() == [2, 1, 3, 1]
    assert nbr[0, :2].tolist() == [2, 1]
    assert nbr[2].tolist() == [1, 0, 3]


def test_independent_production_leaves_wl_state_untouched():
    """K3-Zerlegung: eine Produktion mit eigenem RNG arbeitet auf einer
    KOPIE; zweimal derselbe Stream -> bit-gleich; anderer Stream -> anders."""
    L = 3
    nbr, ei, ej, ax, ay, n, lo, hi = _small_job(L)
    st = phy044.wl_phase(nbr, ei, ej, n, L, lo, hi, lnf_final=1e-2,
                         backend="python")
    th0, E0, b0 = st.th.copy(), st.E, st.b
    runs = []
    for s in (7, 7, 8):
        rng = phy044.make_rng(42, stream=s)
        runs.append(phy044.production_phase(st, nbr, ei, ej, ax, ay, n, L,
                                            120, backend="python", rng=rng))
        assert np.array_equal(st.th, th0) and (st.E, st.b) == (E0, b0)
    assert phy044.results_identical(runs[0], runs[1])
    assert not phy044.results_identical(runs[0], runs[2])


def test_stream_contract_collision_free():
    """PHY044-Streams (3e6 + 1e4*w + L, Zerlegung 3.5e6 + 1e4*p + L)
    kollidieren fuer w, p < 10 und L <= 128 mit keinem bisherigen Vertrag
    (max. PHY043: 900 + s + 1000*L + 100000*t_idx < 2e6)."""
    phy043_max = 900 + 7 + 1000 * 25 + 100000 * 17
    assert phy043_max < phy044._STREAM_CAL_BASE
    cal = {phy044._cal_stream(w) + L for w in range(10) for L in range(129)}
    dec = {phy044._STREAM_DECOMP_BASE + 10_000 * p + L
           for p in range(10) for L in range(129)}
    assert not cal & dec
    assert min(cal) > 90000 + 1000 * 9 + 128        # PHY042-Zusatz-Walker


# ---------------------------------------------------------------------------
# (c) Blind-Vertrag + Auswertung
# ---------------------------------------------------------------------------

def test_assert_blind_rejects_location_keys_at_any_depth():
    phy044.assert_blind({"a": [{"spread": 0.01}], "b": {"n": 3}})
    for bad in ({"tbkt": 0.59}, {"x": [{"pair_tbkt": 1}]},
                {"y": {"z": {"T_BKT_est": 0.5}}}):
        with pytest.raises(AssertionError):
            phy044.assert_blind(bad)


def test_pair_crossing_spread_reports_spread_only():
    t = phy044._p42._T_GRID
    base = 1.3 - 2.0 * (t - 0.52)
    a = [base, base + 0.002]
    b = [base - 0.01 * (t - 0.52) * 10, base - 0.012 * (t - 0.52) * 10]
    out = phy044.pair_crossing_spread(a, 24, b, 32, t)
    assert set(out) == {"n_combos", "n_crossing", "spread"}
    assert out["n_combos"] == 4
    phy044.assert_blind(out)


def test_fit_power_law_oracle():
    Ls = [24, 32, 48, 64]
    ys = [3.0 * L ** 3.5 for L in Ls]
    fit = phy044.fit_power_law(Ls, ys)
    assert fit["p"] == pytest.approx(3.5, abs=1e-12)
    assert fit["A"] == pytest.approx(3.0, rel=1e-9)
    assert all(p == pytest.approx(3.5, abs=1e-12) for p in fit["local_p"])


# ---------------------------------------------------------------------------
# (d) W4-Vorregistrierung
# ---------------------------------------------------------------------------

def test_w4_constants_match_preregistration_spec():
    spec = SPEC.read_text(encoding="utf-8")
    lo, hi = phy044.W4_BAND
    assert f"**B = [{lo:.3f}, {hi:.3f}]**" in spec
    assert f"**>= {phy044.W4_MIN_WALKERS} unabhaengige g(E)-Walker" in spec
    assert f"weniger als **{phy044.W4_MIN_PAIRS}** quotierbare" in spec
    assert f"sigma_tot > **{phy044.W4_SIGMA_MAX:.3f}**" in spec
    assert f"sigma_FSS >= {phy044.W4_SIGMA_FSS_FLOOR:.3f}" in spec
    assert f"< {phy044.W4_DOMAIN_THRESHOLD:.2f} bleibt" in spec
    assert f"**T_req = {phy044.W4_T_REQ:.2f}**" in spec
    assert ("{" + ", ".join(str(L) for L in phy044.W4_LADDER) + "}") in spec
    # Band-Grenzen folgen aus dem Referenzband (Addendum §4)
    band = phy044._p42.REF_BAND
    assert lo == band["nn_mc"][0]
    assert hi == band["binder_beta"][0]
    assert phy044.W4_DOMAIN_THRESHOLD == phy044._p42.DOMAIN_THRESHOLD


def test_w4_verdict_rules_in_order():
    v = phy044.w4_verdict
    assert v(0.575, 0.001, 0.002, 2)["rule"] == "S1"          # zu wenig Paare
    assert v(None, 0.001, 0.002, 5)["rule"] == "S1"
    assert v(0.575, None, 0.002, 3)["rule"] == "S2"
    assert v(0.575, 0.009, 0.006, 3)["rule"] == "S2"          # hypot > 0.010
    ok = v(0.586, 0.002, 0.001, 3)                            # Floor 0.003
    assert ok["sigma_tot"] == pytest.approx(math.hypot(0.002, 0.003))
    assert ok["verdict"] == "CONSISTENT"                      # lo = 0.5788
    far = v(0.590, 0.002, 0.001, 3)                           # lo = 0.5828
    assert far["verdict"] == "INCONSISTENT" and far["rule"] == "I"
    assert v(0.570, 0.0, 0.0, 4)["verdict"] == "CONSISTENT"
    assert v(0.555, 0.0, 0.0, 4)["verdict"] == "CONSISTENT"   # hi = 0.561
    assert v(0.550, 0.0, 0.0, 4)["interval"][1] == pytest.approx(0.556)
    assert v(0.550, 0.0, 0.0, 4)["verdict"] == "INCONSISTENT"


def test_committed_calibration_report_is_blind_if_present():
    """Der committete PHY044-Report (Gate-Log) enthaelt keinen Lagewert."""
    rep = ROOT / "results" / "260926 PHY044 honeycomb wl calibration report.json"
    if not rep.exists():
        pytest.skip("Kalibrier-Report noch nicht committet")
    data = json.loads(rep.read_text(encoding="utf-8"))
    phy044.assert_blind(data)
    assert data["overall_pass"] is True
    assert data["pass_gates"]["PASS_VAL_BIT_PHY042_ALL_7_WALKERS"] is True
    assert data["pass_gates"]["PASS_BIT_IDENTICAL_L24_SPEED_JOB"] is True
    # Vorregistrierte Konstanten sind die, mit denen gerechnet wurde
    assert data["T_req"] == phy044.W4_T_REQ
    assert data["domain_threshold"] == phy044.W4_DOMAIN_THRESHOLD
    assert data["val_bit_phy042"]["source_sha256"] == (
        "19a9ce3c799401dbb55e519c09b390bee9dfbe23792f99b2e7c650c3eeefa3cf")
    # >= 3 Walker an jedem L der Leiter; BLAS-Threads gepinnt protokolliert
    assert all(b["n_walkers"] >= phy044.W4_MIN_WALKERS
               for b in data["K2_precision"]["per_L"].values())
    assert all(m["blas_threads_env"]["OPENBLAS_NUM_THREADS"] == "1"
               for m in data["stage_meta"].values())


def _fake_job(kind, L, w, **kw):
    """Synthetischer Stufen-Job (Form wie _job) fuer den Report-Pfad."""
    t = phy044._p42._T_GRID
    y2 = 1.3 - 2.0 * (t - 0.52) - 0.01 * math.log(L)
    y2 = y2 + np.where(t > 0.60, 0.02 * w, 0.0)   # Spread 0.04 ab T>0.60
    job = {"kind": kind, "L": L, "walker": w, "n": 2 * L * L,
           "nbins": 10 * L, "stream": 0, "lnf_final": kw.get("lnf", 1e-5),
           "prod_sweeps": 1000 * kw.get("pf", 1), "prod_factor": kw.get("pf", 1),
           "wl_sweeps": 100000 + 10 * L ** 2, "sweeps_1t": 5000,
           "t_wl_s": 1e-3 * L ** 3, "t_prod_s": 1e-4 * L ** 3,
           "backend": "numba", "covered": [10 * L, 10 * L], "leak_max": 0.0,
           "uncovered_mass": [0.0] * len(t),
           "curve": {"T": t.tolist(), "y2": y2.tolist(),
                     "y4_scaled": (0 * t).tolist(), "E": (0 * t).tolist()}}
    if kind == "decomp":
        job["curves"] = [job["curve"]] * 3
        job["uncovered_mass"] = [job["uncovered_mass"]] * 3
    return job


def test_analyse_and_text_report_are_blind_and_complete(tmp_path):
    ladder = [_fake_job("cal", L, w) for L in (24, 32, 48, 64)
              for w in range(3)]
    levers = ([_fake_job("cal", 48, w, pf=4) for w in range(3)]
              + [_fake_job("cal", 48, w, lnf=1e-6) for w in range(3)]
              + [_fake_job("decomp", 48, 0)])
    stages = {"ladder": {"jobs": ladder, "wall_s": 1.0, "max_workers": 4,
                         "windows": {}},
              "levers": {"L": 48, "jobs": levers},
              "speed": {"t_python_s": 10.0, "t_numba_s": 1.0,
                        "speedup": 10.0, "python_updates_per_s": 3e5,
                        "numba_updates_per_s": 3e6, "L": 24, "n": 1152,
                        "wl_sweeps": 1, "prod_sweeps": 1, "spin_updates": 1,
                        "PASS_BIT_IDENTICAL_L24": True}}
    rep = phy044.analyse(stages)
    phy044.assert_blind(rep)
    assert rep["K2_precision"]["per_L"]["48"]["domain"]["tmax"] == \
        pytest.approx(0.60)
    assert rep["K2_precision"]["per_L"]["48"]["meets_T_req"] is False
    assert rep["K3_levers"]["w4_go_L_ge_48"] is False         # Stop-Regel S0
    assert rep["K1_cost"]["fit_wall_numba"]["p"] == pytest.approx(3.0)
    assert rep["pass_gates"]["PASS_WALKER_CONTRACT_GE3_EVERY_L"] is True
    assert rep["pass_gates"]["PASS_NO_UNCOVERED_MASS_IN_DOMAIN_LADDER"] is True
    assert rep["findings"]["one_over_t_engaged_per_L"]["64"] == [True] * 3
    out = tmp_path / "r.txt"
    phy044.write_text_report(phy044._clean(rep), out)
    text = out.read_text(encoding="utf-8")
    assert "OVERALL: PASS" in text and "W4-GO fuer L>=48: False" in text


def test_uncovered_mass_gate_binds_only_inside_domain():
    """PHY042-Vertrag: unbesetzte kanonische Masse zaehlt nur in-Domaene;
    ausserhalb (hier T > 0.60, wo der Walker-Spread 0.04 erreicht) ist
    sie Befund, kein Gate."""
    t = phy044._p42._T_GRID
    jobs = [_fake_job("cal", 48, w) for w in range(3)]
    for j in jobs:
        j["uncovered_mass"] = [0.5 if T > 0.60 + 1e-9 else 1e-6 for T in t]
    b = phy044._walker_block(jobs, t.tolist())
    assert b["domain"]["tmax"] == pytest.approx(0.60)
    assert b["uncovered_mass_in_domain_max"] == pytest.approx(1e-6)
    jobs[1]["uncovered_mass"][0] = 0.01              # in-Domaene verletzt
    assert phy044._walker_block(jobs, t.tolist())[
        "uncovered_mass_in_domain_max"] == pytest.approx(0.01)
