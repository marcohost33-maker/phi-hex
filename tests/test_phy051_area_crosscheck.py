"""Gates fuer PHY051 (triangular + kagome, Upsilon pro Flaeche, O1-Diskriminator).

(a) Vorregistrierte Konstanten zeichengenau an die Spec v01 gebunden.
(b) Seed-Vertrag kollisionsfrei und disjunkt zu PHY045/046; Paare der Leiter.
(c) T=0-Geometrie-Orakel exakt (MC-frei).
(d) Schaetzer-Orakel auf synthetischen WM-Kurven + Wahrheitstafel fails-closed.
(e) Mini-Produktion end-to-end (klein, schnell); committete Reports gegen Plan.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

import conftest

phy051 = conftest._load(
    "phy051_area_crosscheck",
    "260928 PHY051 triangular kagome area-helicity crosscheck v01.py")

ROOT = Path(__file__).resolve().parents[1]
SPEC = (ROOT / "spec"
        / "260928 PHI HEX phy051 triangular kagome per-area preregistration v01.md")


def test_constants_match_spec_v01():
    spec = SPEC.read_text(encoding="utf-8")
    for name, plan in phy051.PHY051_PLAN.items():
        assert ("L in {" + ", ".join(str(L) for L in plan["ladder"]) + "}"
                in spec), name
        g = plan["t_grid"]
        assert f"{g[0]:.3f} .. {g[-1]:.3f}, Schritt 0.005" in spec, name
        assert np.allclose(np.diff(g), 0.005)
        assert (f"[{plan['band'][0]:.3f}, {plan['band'][1]:.3f}]" in spec), name
        assert f"L_min = {plan['L_min_primary']}" in spec, name
        for a, b in plan["pairs"]:
            assert f"({a},{b})" in spec, (name, a, b)
    assert len(phy051.PHY051_PLAN["triangular"]["t_grid"]) == 27
    assert len(phy051.PHY051_PLAN["kagome"]["t_grid"]) == 17
    assert f"{phy051.PHY051_N_SEEDS} unabhaengige Seeds" in spec
    assert (f"je {phy051.PHY051_N_THERM} Thermalisierungs-" in spec
            and f"{phy051.PHY051_N_MEAS} Mess-Sweeps" in spec)
    assert "`51_000_000 + 1_000_000*lat_id + 1000*L + 10*t_idx + s`" in spec
    assert phy051.PHY051_SEED_BASE == 51_000_000
    assert "Untergrenze **0.5 %" in spec
    assert phy051.PHY051_SIGMA_FSS_FLOOR_REL == 0.005
    assert phy051.PHY051_SIGMA_MAX_REL == 0.010
    assert "3.0 h auf 4 Prozessen" in spec
    assert phy051.PHY051_WALL_BUDGET_H == 3.0
    assert phy051.PHY051_NORMALIZATION == "per_area"
    assert "1.4654(17)" in spec and "1.418(2)" in spec and "0.825" in spec


def test_seed_contract_unique_and_disjoint_from_predecessors():
    seeds = set()
    for name, plan in phy051.PHY051_PLAN.items():
        for L in plan["ladder"]:
            for k in range(len(plan["t_grid"])):
                for s in range(phy051.PHY051_N_SEEDS):
                    seeds.add(phy051.seed_for(name, L, k, s))
    n_expected = sum(len(p["ladder"]) * len(p["t_grid"])
                     * phy051.PHY051_N_SEEDS
                     for p in phy051.PHY051_PLAN.values())
    assert len(seeds) == n_expected
    assert min(seeds) >= 51_000_000 and max(seeds) < 2 ** 32
    # PHY045 (45M), PHY045-C (46M), PHY046 (47M): alle < 48M
    assert min(seeds) >= 48_000_000


def test_pairs_belong_to_ladder_and_are_near_doubling():
    for name, plan in phy051.PHY051_PLAN.items():
        lad = set(plan["ladder"])
        for a, b in plan["pairs"]:
            assert a in lad and b in lad, name
            assert 1.95 <= b / a <= 2.0, (name, a, b)
        assert plan["L_min_primary"] in lad
        n_prim = sum(1 for a, _ in plan["pairs"] if a >= plan["L_min_primary"])
        assert n_prim >= phy051.PHY051_MIN_PRIMARY_PAIRS
    tri = phy051.PHY051_PLAN["triangular"]["ladder"]
    assert all(L % 2 == 1 for L in tri)   # Core-Torus L = 2r+1


def test_t0_geometry_oracle_exact():
    o = phy051.geometry_oracle("triangular", (33,))
    assert o["33"]["upsilon0_site"] == pytest.approx(1.5, abs=1e-12)
    assert o["33"]["upsilon0_area"] == pytest.approx(math.sqrt(3.0), abs=1e-12)
    assert o["33"]["area_per_site"] == pytest.approx(math.sqrt(3.0) / 2, abs=1e-12)
    k = phy051.geometry_oracle("kagome", (32,))
    assert k["32"]["upsilon0_site"] == pytest.approx(1.0, abs=1e-12)
    assert k["32"]["upsilon0_area"] == pytest.approx(math.sqrt(3.0) / 2, abs=1e-12)
    assert k["32"]["n"] == 3 * 32 * 32


def _wm_curves(T0, C, t_grid, ladder, slope):
    """Synthetische per-Flaeche-Kurven, die bei T0 exakt der WM-Form folgen."""
    return {L: np.array([(2 * T0 / math.pi) * (1 + 1 / (2 * math.log(L) + C))
                         - slope * (T - T0) for T in t_grid]) for L in ladder}


def _synthetic_prod(lattice, T0, noise=2e-4, seed=5, ladder=None,
                    t_grid=None):
    """Synthetik: per-Flaeche-Kurven exakt WM-foermig bei T0, physikalisch
    plausible Steigung (Abfall um 2T0/pi ueber ~0.4 T0, wie square/honeycomb
    in PHY046/048); per Site = a_s * per Flaeche. Breites T-Gitter, damit
    beide Kanaele sicher im Fenster crossen (die Auswertung ist gitter-
    agnostisch; der Plan-T-Gitter-Vertrag wird separat getestet)."""
    plan = phy051.PHY051_PLAN[lattice]
    ladder = list(plan["ladder"] if ladder is None else ladder)
    t = list(plan["t_grid"] if t_grid is None else t_grid)
    a_s = phy051.AREA_PER_SITE[lattice]
    slope = 2.0 * T0 / math.pi / (0.4 * T0)
    rng = np.random.default_rng(seed)
    base = _wm_curves(T0, 2.0, t, ladder, slope)
    ns = 4
    area = {str(L): (base[L][:, None]
                     + noise * rng.standard_normal((len(t), ns))).tolist()
            for L in ladder}
    site = {str(L): (a_s * np.asarray(area[str(L)])).tolist() for L in ladder}
    return {"lattice": lattice, "t_grid": t, "ladder": ladder, "n_seeds": ns,
            "n_therm": 1, "n_meas": 1, "ups_area": area, "ups_site": site,
            "unmeasured_L": [], "wall_s": 0.0, "cpu_s": 0.0,
            "max_workers": 1, "have_numba": False}


def test_analyse_recovers_wm_temperature_and_o1_truth_table_triangular():
    """Exakte WM-Kurven pro Flaeche bei T0 in B_tri: per-Flaeche-Kanal muss
    T0 liefern (CONSISTENT); der per-Site-Kanal (x a_s = 0.866) crosst die
    NK-Gerade tiefer und liegt ausserhalb B -> O1_PER_AREA_CORROBORATED."""
    prod = _synthetic_prod("triangular", T0=1.465,
                           t_grid=phy051._grid(1.30, 1.60))
    area = phy051.analyse(prod, "ups_area")
    assert area["estimates"]["n_pairs_with_crossing"] == 5
    for T in area["estimates"]["pairs"].values():
        assert T == pytest.approx(1.465, abs=2e-3)
    assert area["T_P"] == pytest.approx(1.465, abs=2e-3)
    assert area["verdict"]["verdict"] == "CONSISTENT"
    site = phy051.analyse(prod, "ups_site")
    assert site["T_P"] is not None and site["T_P"] < 1.45
    assert site["verdict"]["verdict"] == "INCONSISTENT"
    rep = phy051._clean(phy051.build_report(prod))
    assert rep["o1_label"] == "O1_PER_AREA_CORROBORATED"
    json.dumps(rep, allow_nan=False)


def test_o1_challenged_when_per_site_curves_hit_the_band():
    """Umgekehrtes Szenario: die per-SITE-Kurven folgen der WM-Form bei T0 in B
    (d.h. die Flaechen-Kurven sind um 1/a_s zu hoch) -> per Site CONSISTENT,
    per Flaeche INCONSISTENT -> O1_CHALLENGED."""
    prod = _synthetic_prod("triangular", T0=1.465,
                           t_grid=phy051._grid(1.30, 1.60))
    a_s = phy051.AREA_PER_SITE["triangular"]
    site_curves = prod["ups_area"]                       # WM-Form bei T0
    prod["ups_site"] = site_curves
    prod["ups_area"] = {L: (np.asarray(v) / a_s).tolist()
                        for L, v in site_curves.items()}
    rep = phy051._clean(phy051.build_report(prod))
    assert rep["convention_crosscheck_per_site"]["verdict"]["verdict"] == "CONSISTENT"
    assert rep["primary_per_area"]["verdict"]["verdict"] == "INCONSISTENT"
    assert rep["o1_label"] == "O1_CHALLENGED"


def test_stop_rules_fail_closed():
    band = (1.45, 1.48)
    assert phy051.verdict(1.46, 0.001, 0.001, 5, 3, band,
                          ladder_complete=False)["rule"] == "S0"
    assert phy051.verdict(1.46, 0.001, 0.001, 2, 2, band)["rule"] == "S1"
    assert phy051.verdict(None, 0.001, 0.001, 5, 3, band)["rule"] == "S1"
    assert phy051.verdict(1.46, 0.001, 0.001, 3, 1, band)["rule"] == "S1b"
    assert phy051.verdict(1.46, None, 0.001, 5, 3, band)["rule"] == "S2"
    assert phy051.verdict(1.46, 0.02, 0.001, 5, 3, band)["rule"] == "S2"
    ok = phy051.verdict(1.46, 0.001, 0.001, 5, 3, band)
    assert ok["verdict"] == "CONSISTENT" and ok["rule"] == "C"
    assert ok["sigma_fss_floored"] == pytest.approx(0.005 * 1.46)
    bad = phy051.verdict(1.40, 0.001, 0.001, 5, 3, band)
    assert bad["verdict"] == "INCONSISTENT" and bad["rule"] == "I"
    # Wahrheitstafel
    assert phy051.o1_label("CONSISTENT", "INCONSISTENT") == "O1_PER_AREA_CORROBORATED"
    assert phy051.o1_label("INCONSISTENT", "CONSISTENT") == "O1_CHALLENGED"
    assert phy051.o1_label("CONSISTENT", "CONSISTENT") == "NON_DISCRIMINATING"
    assert phy051.o1_label("INCONSISTENT", "INCONSISTENT") == "TENSION_BOTH_CHANNELS"
    assert phy051.o1_label("NEGATIVE_RESULT", "CONSISTENT") == "NEGATIVE_RESULT"
    assert phy051.o1_label("CONSISTENT", "NEGATIVE_RESULT") == "NEGATIVE_RESULT"


def test_analyse_with_incomplete_ladder_is_negative_s0():
    prod = _synthetic_prod("kagome", T0=0.825, ladder=[32, 48, 64, 96],
                           t_grid=phy051._grid(0.70, 0.95))
    out = phy051.analyse(prod, "ups_area")
    assert out["ladder_complete"] is False
    assert out["verdict"]["rule"] == "S0"


def test_weighted_pair_mean_contract():
    pc = {(33, 65): 1.40, (65, 129): 1.46, (129, 257): 1.47}
    sig = {(33, 65): 0.001, (65, 129): 0.001, (129, 257): 0.003}
    r = phy051.weighted_pair_mean(pc, sig, 65)
    assert r["n_pairs"] == 2 and r["weighted"] is True
    assert r["T"] == pytest.approx((1.46 / 1e-6 + 1.47 / 9e-6)
                                   / (1 / 1e-6 + 1 / 9e-6))
    r2 = phy051.weighted_pair_mean(pc, {(65, 129): 0.0}, 65)
    assert r2["weighted"] is False and r2["T"] == pytest.approx(1.465)
    assert phy051.weighted_pair_mean({(33, 65): None}, None, 65)["T"] is None


def test_mini_production_end_to_end(tmp_path):
    prod = phy051.produce("triangular", ladder=(5, 9), t_grid=(1.3, 1.5),
                          n_seeds=2, n_therm=3, n_meas=8, max_workers=1)
    assert np.asarray(prod["ups_area"]["9"]).shape == (2, 2)
    ratio = (np.asarray(prod["ups_site"]["9"])
             / np.asarray(prod["ups_area"]["9"]))
    assert np.allclose(ratio, phy051.AREA_PER_SITE["triangular"])
    assert prod["unmeasured_L"] == []
    rep = phy051._clean(phy051.build_report(prod))
    json.dumps(rep, allow_nan=False)
    assert rep["o1_label"] == "NEGATIVE_RESULT"       # Leiter != Plan -> S0
    phy051.write_text_report(rep, tmp_path / "r.txt")
    txt = (tmp_path / "r.txt").read_text(encoding="utf-8")
    assert "O1-Diskriminator" in txt and "Gates" in txt


def test_budget_stop_marks_remaining_L_unmeasured():
    prod = phy051.produce("kagome", ladder=(4, 8), t_grid=(0.8,), n_seeds=1,
                          n_therm=1, n_meas=2, max_workers=1,
                          wall_budget_h=0.0)
    # erster Block (L=8) laeuft immer, danach greift das Budget
    assert prod["unmeasured_L"] == [4]
    assert np.all(np.isnan(np.asarray(prod["ups_area"]["4"])))
    out = phy051.analyse(prod, "ups_area")
    assert out["verdict"]["rule"] == "S0"


@pytest.mark.parametrize("lattice", ["triangular", "kagome"])
def test_committed_report_if_present(lattice):
    p = ROOT / "results" / (phy051.report_stem(lattice) + ".json")
    if not p.exists():
        pytest.skip("PHY051-Report noch nicht committet")
    rep = json.loads(p.read_text(encoding="utf-8"))
    plan = phy051.PHY051_PLAN[lattice]
    pr = rep["production"]
    assert pr["ladder"] == list(plan["ladder"])
    assert pr["t_grid"] == pytest.approx(list(plan["t_grid"]))
    assert pr["n_seeds"] == phy051.PHY051_N_SEEDS
    assert pr["n_meas"] == phy051.PHY051_N_MEAS
    assert rep["o1_label"] in ("O1_PER_AREA_CORROBORATED", "O1_CHALLENGED",
                               "NON_DISCRIMINATING", "TENSION_BOTH_CHANNELS",
                               "NEGATIVE_RESULT")
    # Report ist deterministisch aus den gespeicherten Rohdaten reproduzierbar
    prod = {"lattice": lattice, "t_grid": pr["t_grid"], "ladder": pr["ladder"],
            "n_seeds": pr["n_seeds"], "n_therm": pr["n_therm"],
            "n_meas": pr["n_meas"], "unmeasured_L": pr["unmeasured_L"],
            "wall_s": pr["wall_s"], "cpu_s": pr["cpu_s"],
            "max_workers": pr["max_workers"], "have_numba": pr["have_numba"],
            "ups_area": rep["data"]["ups_area"],
            "ups_site": rep["data"]["ups_site"]}
    again = phy051._clean(phy051.build_report(prod))
    assert again["o1_label"] == rep["o1_label"]
    assert again["primary_per_area"]["T_P"] == pytest.approx(
        rep["primary_per_area"]["T_P"], abs=1e-9)
    assert again["pass_gates"] == rep["pass_gates"]
