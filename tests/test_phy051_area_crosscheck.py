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


def _channel_label(rep):
    """Ungegatete Wahrheitstafel aus den im Report stehenden Kanal-Verdikten."""
    return phy051.o1_label(rep["primary_per_area"]["verdict"]["verdict"],
                           rep["convention_crosscheck_per_site"]["verdict"]["verdict"])


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
    # Wahrheitstafel auf Kanal-Ebene; das veroeffentlichte Label ist an die
    # Integritaets-Gates gebunden (Synthetik: 4 Seeds, Nicht-Plan-Gitter ->
    # PASS_INPUT_COMPLETE False -> NEGATIVE_RESULT, Codex #54 Z.420/463)
    assert _channel_label(rep) == "O1_PER_AREA_CORROBORATED"
    assert rep["pass_gates"]["PASS_INPUT_COMPLETE"] is False
    assert rep["o1_label"] == "NEGATIVE_RESULT"
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
    assert _channel_label(rep) == "O1_CHALLENGED"
    assert rep["o1_label"] == "NEGATIVE_RESULT"    # Gates (Nicht-Plan-Eingabe)


def test_stop_rules_fail_closed():
    band = (1.45, 1.48)
    assert phy051.verdict(1.46, 0.001, 0.001, 5, 3, band,
                          ladder_complete=False)["rule"] == "S0"
    assert phy051.verdict(1.46, 0.001, 0.001, 2, 2, band)["rule"] == "S1"
    assert phy051.verdict(None, 0.001, 0.001, 5, 3, band)["rule"] == "S1b"
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


def test_v2_is_the_largest_pair_alone_and_diagnostic_is_not_decisive():
    """Spec 5 (v2): groesstes Paar ALLEIN - fehlt dessen Crossing, ist v2 None
    (nicht das naechstkleinere Paar). Nachtrag v01a: die Diagnostik ohne v3
    traegt decisive=False und aendert das Verdikt nicht."""
    prod = _synthetic_prod("triangular", T0=1.465,
                           t_grid=phy051._grid(1.30, 1.60))
    # groesstes Paar (129,257) durch Rauschen unbrauchbar machen: NaN-freie,
    # aber konstant zu hohe Kurve fuer L=257 -> kein WM-Crossing
    prod["ups_area"]["257"] = (np.asarray(prod["ups_area"]["257"]) + 5.0).tolist()
    out = phy051.analyse(prod, "ups_area")
    e = out["estimates"]
    assert e["pairs"]["129_257"] is None
    assert e["v2_largest_pair"] is None
    assert e["n_pairs_with_crossing"] == 4
    assert out["diagnostic_without_hks3par"]["decisive"] is False
    assert out["verdict"]["verdict"] in ("CONSISTENT", "INCONSISTENT",
                                         "NEGATIVE_RESULT")


def test_stop_rule_s1b_when_primary_missing_but_three_pairs():
    band = (1.45, 1.48)
    assert phy051.verdict(None, 0.001, 0.001, 3, 0, band)["rule"] == "S1b"
    assert phy051.verdict(1.46, 0.001, 0.001, 3, 1, band)["rule"] == "S1b"
    assert phy051.verdict(1.46, 0.001, 0.001, 2, 2, band)["rule"] == "S1"


def test_regenerate_roundtrip(tmp_path):
    """Report -> prod_from_report -> build_report ist idempotent."""
    prod = _synthetic_prod("kagome", T0=0.825, t_grid=phy051._grid(0.70, 0.95))
    rep = phy051._clean(phy051.build_report(prod))
    stem = phy051.report_stem("kagome")
    (tmp_path / f"{stem}.json").write_text(json.dumps(rep), encoding="utf-8")
    again = phy051.regenerate("kagome", out_dir=tmp_path)
    assert again["o1_label"] == rep["o1_label"]
    assert again["primary_per_area"]["T_P"] == pytest.approx(
        rep["primary_per_area"]["T_P"], abs=1e-12)
    assert again["pass_gates"] == rep["pass_gates"]
    assert (tmp_path / f"{stem}.txt").exists()


def test_pool_is_fork_only_and_sequential_otherwise(monkeypatch):
    import multiprocessing
    assert phy051._pool(1) is None
    ex = phy051._pool(2)
    if ex is not None:                       # Plattform mit fork
        assert ex._mp_context.get_start_method() == "fork"
        ex.shutdown()

    def _no_fork(method):
        raise ValueError("cannot find context for 'fork'")
    monkeypatch.setattr(multiprocessing, "get_context", _no_fork)
    assert phy051._pool(4) is None           # fail-safe: sequenziell


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
    assert prod["run_utc_start"] <= prod["run_utc_end"]
    rep = phy051._clean(phy051.build_report(prod))
    json.dumps(rep, allow_nan=False)
    assert rep["production"]["run_utc_end"] == prod["run_utc_end"]
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
def test_committed_report_if_present(lattice, tmp_path):
    p = ROOT / "results" / (phy051.report_stem(lattice) + ".json")
    if not p.exists():
        pytest.skip("PHY051-Report noch nicht committet")
    rep = json.loads(p.read_text(encoding="utf-8"))
    plan = phy051.PHY051_PLAN[lattice]
    pr = rep["production"]
    assert pr["ladder"] == list(plan["ladder"])
    assert pr["t_grid"] == pytest.approx(list(plan["t_grid"]))
    assert pr["n_seeds"] == phy051.PHY051_N_SEEDS
    assert pr["n_therm"] == phy051.PHY051_N_THERM
    assert pr["n_meas"] == phy051.PHY051_N_MEAS
    assert rep["o1_label"] in ("O1_PER_AREA_CORROBORATED", "O1_CHALLENGED",
                               "NON_DISCRIMINATING", "TENSION_BOTH_CHANNELS",
                               "NEGATIVE_RESULT")
    # Report ist deterministisch aus den gespeicherten Rohdaten reproduzierbar:
    # VOLLSTAENDIGER Vergleich (Codex #54 P2, Test-Z.313) - inkl. Lauf-
    # Zeitstempel, Crossings, Unsicherheiten, Verdikt-Details, Referenzen und
    # Rohdaten; Gleitkomma bis rel 1e-9 (plattformabhaengige Fits), alles
    # andere exakt.
    again = phy051._clean(phy051.build_report(phy051.prod_from_report(rep)))
    _assert_report_equal(again, rep)
    # Text-Report: aus dem regenerierten JSON byte-gleich (modulo Zeilenende)
    out = tmp_path / "r.txt"
    phy051.write_text_report(again, out)
    committed_txt = p.with_suffix(".txt").read_bytes().replace(b"\r\n", b"\n")
    assert out.read_bytes().replace(b"\r\n", b"\n") == committed_txt


def _assert_report_equal(a, b, path="$"):
    if isinstance(b, float) and not isinstance(a, bool) and isinstance(a, (int, float)):
        assert a == pytest.approx(b, rel=1e-9, abs=1e-12), path
        return
    assert type(a) is type(b), (path, type(a), type(b))
    if isinstance(b, dict):
        assert set(a) == set(b), (path, set(a) ^ set(b))
        for k in b:
            _assert_report_equal(a[k], b[k], f"{path}.{k}")
    elif isinstance(b, list):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            _assert_report_equal(x, y, f"{path}[{i}]")
    else:
        assert a == b, path


def test_report_comparison_catches_a_stale_field():
    """Positivkontrolle des Vollvergleichs: eine einzige veraltete Zahl tief
    im Report (Jackknife-SE eines Paars) oder ein fehlender Zeitstempel muss
    auffallen."""
    p = ROOT / "results" / (phy051.report_stem("kagome") + ".json")
    rep = json.loads(p.read_text(encoding="utf-8"))
    stale = json.loads(json.dumps(rep))
    k = next(iter(stale["primary_per_area"]["pair_jackknife_se"]))
    stale["primary_per_area"]["pair_jackknife_se"][k] *= 1.001
    with pytest.raises(AssertionError):
        _assert_report_equal(stale, rep)
    no_ts = json.loads(json.dumps(rep))
    no_ts["production"]["run_utc_start"] = None
    with pytest.raises(AssertionError):
        _assert_report_equal(no_ts, rep)


# --- Codex #54 (8c92301): eingebaute Befunde -------------------------------

def _plan_prod(lattice="kagome"):
    """Plan-konforme Produktions-Metadaten + Kanalformen (Werte synthetisch)."""
    plan = phy051.PHY051_PLAN[lattice]
    nt, ns = len(plan["t_grid"]), phy051.PHY051_N_SEEDS
    rows = {str(L): np.ones((nt, ns)).tolist() for L in plan["ladder"]}
    return {"lattice": lattice, "t_grid": [float(t) for t in plan["t_grid"]],
            "ladder": list(plan["ladder"]), "n_seeds": ns,
            "n_therm": phy051.PHY051_N_THERM, "n_meas": phy051.PHY051_N_MEAS,
            "ups_area": rows, "ups_site": json.loads(json.dumps(rows)),
            "unmeasured_L": []}


@pytest.mark.parametrize("lattice", ["triangular", "kagome"])
def test_input_matches_plan_positive_control(lattice):
    assert phy051._input_matches_plan(_plan_prod(lattice)) is True


@pytest.mark.parametrize("mutate", [
    "lattice", "ladder", "ladder_float", "t_grid_shift", "t_grid_nan",
    "t_grid_short", "n_seeds", "n_seeds_bool", "n_therm", "n_meas",
    "site_shape", "site_missing_L", "area_extra_L"])
def test_input_matches_plan_rejects_deviation(mutate):
    """Codex #54 P2 (Z.420): Abweichung vom registrierten Plan ->
    PASS_INPUT_COMPLETE darf nicht gelten, auch wenn prod in sich stimmig ist."""
    prod = _plan_prod("kagome")
    nt = len(prod["t_grid"])
    if mutate == "lattice":
        prod["lattice"] = "square"
    elif mutate == "ladder":
        prod["ladder"] = prod["ladder"][:-1]
    elif mutate == "ladder_float":
        prod["ladder"] = [float(L) for L in prod["ladder"]]
    elif mutate == "t_grid_shift":
        prod["t_grid"] = [t + 0.005 for t in prod["t_grid"]]
    elif mutate == "t_grid_nan":
        prod["t_grid"][3] = math.nan
    elif mutate == "t_grid_short":
        prod["t_grid"] = prod["t_grid"][:-1]
    elif mutate == "n_seeds":
        prod["n_seeds"] = 4
    elif mutate == "n_seeds_bool":
        prod["n_seeds"] = True
    elif mutate == "n_therm":
        prod["n_therm"] = 499
    elif mutate == "n_meas":
        prod["n_meas"] = 1500.0
    elif mutate == "site_shape":
        prod["ups_site"]["32"] = np.ones((nt, 4)).tolist()
    elif mutate == "site_missing_L":
        del prod["ups_site"]["192"]
    elif mutate == "area_extra_L":
        prod["ups_area"]["256"] = np.ones((nt, 8)).tolist()
    assert phy051._input_matches_plan(prod) is False


def test_gated_o1_label_suppresses_affirmative_labels_on_gate_fail():
    """Codex #54 P2 (Z.463): faellt ein Integritaets-Gate, kein bejahendes
    O1-Label. Positivkontrolle: bei overall_pass True gilt die Tafel."""
    pairs = [("CONSISTENT", "INCONSISTENT"), ("INCONSISTENT", "CONSISTENT"),
             ("CONSISTENT", "CONSISTENT"), ("INCONSISTENT", "INCONSISTENT")]
    for a, s in pairs:
        assert phy051.gated_o1_label(a, s, True) == phy051.o1_label(a, s)
        assert phy051.gated_o1_label(a, s, True) != "NEGATIVE_RESULT"
        for bad in (False, None, 1, "True", np.bool_(True)):
            assert phy051.gated_o1_label(a, s, bad) == "NEGATIVE_RESULT", bad


def test_build_report_label_follows_gates(monkeypatch):
    """Verdrahtung: build_report nutzt das gegatete Label. Mit erzwungen
    bestandenen Gates erscheint das bejahende Label (Positivkontrolle),
    mit einem gefallenen Gate NEGATIVE_RESULT; Kanal-Verdikte bleiben."""
    prod = _synthetic_prod("triangular", T0=1.465,
                           t_grid=phy051._grid(1.30, 1.60))
    real = phy051.gates

    def _force(passed):
        def g(*a, **k):
            out = real(*a, **k)
            out["gates"] = {n: passed for n in out["gates"]}
            out["overall_pass"] = all(out["gates"].values())
            return out
        return g
    monkeypatch.setattr(phy051, "gates", _force(True))
    ok = phy051._clean(phy051.build_report(prod))
    assert ok["o1_label"] == "O1_PER_AREA_CORROBORATED"
    monkeypatch.setattr(phy051, "gates", _force(False))
    bad = phy051._clean(phy051.build_report(prod))
    assert bad["o1_label"] == "NEGATIVE_RESULT"
    assert _channel_label(bad) == "O1_PER_AREA_CORROBORATED"


@pytest.mark.parametrize("pos", [0, 3, -1])
def test_jackknife_se_fails_closed_on_missing_replicate(pos):
    """Codex #54 P2 (Z.285): ein fehlendes/nicht-endliches Replikat ->
    keine SE (None), statt SE aus den Ueberlebenden."""
    vals = [1.46, 1.47, 1.45, 1.465, 1.455, 1.468, 1.462, 1.459]
    assert phy051._jk_se(vals) is not None and phy051._jk_se(vals) > 0.0
    for bad in (None, math.nan, math.inf):
        v = list(vals)
        v[pos] = bad
        assert phy051._jk_se(v) is None, (pos, bad)
    assert phy051._jk_se([1.46]) is None


def test_blas_threads_are_forced_not_defaulted():
    """Codex #54 P2 (Z.30): geerbte BLAS-Thread-Variablen werden auf 1
    UEBERSCHRIEBEN (Quelltext-Pruefung: der Import geschieht vor numpy,
    ein Laufzeit-Test im selben Prozess saehe nur den Endzustand)."""
    import ast
    src = (ROOT / "src" / ("260928 PHY051 triangular kagome area-helicity "
                           "crosscheck v01.py")).read_text(encoding="utf-8")
    tree = ast.parse(src)
    loops = [n for n in tree.body if isinstance(n, ast.For)
             and isinstance(n.iter, ast.Tuple)
             and {getattr(e, "value", None) for e in n.iter.elts}
             == {"OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"}]
    assert len(loops) == 1
    body = loops[0].body
    assert len(body) == 1 and isinstance(body[0], ast.Assign)
    assert ast.unparse(body[0]) == "os.environ[_v] = '1'"
    first_numpy = min(n.lineno for n in tree.body if isinstance(n, ast.Import)
                      and any(a.name == "numpy" for a in n.names))
    assert loops[0].lineno < first_numpy
    import os
    for v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        assert os.environ[v] == "1"
