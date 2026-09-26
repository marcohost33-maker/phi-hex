"""Gates fuer PHY045 (Audit O1: Helicity-Normierung im NK-Kriterium).

(a) Geometrie: Flaeche/Site und Upsilon_site(0) der Repo-Builder.
(b) Teil A exakt (MC-frei, klein): die langreichweitige Steifigkeit des
    harmonischen Gitters ist Upsilon pro FLAECHE, nicht pro Site.
(c) Numba-Wolff: T->0-Orakel je Gitter (kurzer Lauf).
(d) Entscheid-Logik und committeter Report.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

import conftest

phy045 = conftest._load(
    "phy045_normalization", "260926 PHY045 helicity normalization O1 test v01.py")

ROOT = Path(__file__).resolve().parents[1]
LATTICES = ("square", "triangular", "honeycomb", "kagome")


def _L(name: str, L: int) -> int:
    return L + 1 if name == "triangular" else L


@pytest.mark.parametrize("name", LATTICES)
def test_builder_geometry_matches_constants(name):
    lat = phy045.build(name, _L(name, 6))
    assert lat.area_per_site == pytest.approx(phy045.AREA_PER_SITE[name],
                                              rel=1e-12)
    assert phy045.upsilon0_site_from_geometry(lat) == pytest.approx(
        phy045.UPSILON0_SITE[name], rel=1e-12)


def test_area_factors_are_the_known_closed_forms():
    a = phy045.AREA_PER_SITE
    assert a["triangular"] == pytest.approx(math.sqrt(3) / 2)
    assert a["honeycomb"] == pytest.approx(3 * math.sqrt(3) / 4)
    assert a["kagome"] == pytest.approx(2 / math.sqrt(3))


def test_part_a_harmonic_stiffness_is_per_area_not_per_site():
    """Exakter Kern von Audit O1 (klein, CI-schnell): K_fit / Upsilon_area(0)
    = 1 auf <= 1 %, K_fit / Upsilon_site(0) = 1/a_s auf allen Gittern."""
    out = phy045.part_a(L=36)
    for name, v in out.items():
        assert v["ratio_K_fit_over_area"] == pytest.approx(1.0, abs=0.01), name
        assert v["ratio_K_fit_over_site"] == pytest.approx(
            1.0 / phy045.AREA_PER_SITE[name], rel=0.01), name


@pytest.mark.parametrize("name", LATTICES)
def test_numba_wolff_t0_oracle(name):
    r = phy045.mc_point((name, _L(name, 6), 0.002, 11, 20, 100))
    assert r["upsilon_site"] == pytest.approx(phy045.UPSILON0_SITE[name],
                                              rel=2e-3)


def test_eta_fit_recovers_power_law():
    Ls = [16, 32, 64, 128]
    m2 = [0.9 * L ** -0.08 * (1 + 0.3 / L ** 2) for L in Ls]
    err = [1e-5] * 4
    assert phy045.eta_fit_corrected(Ls, m2, err)["eta"] == pytest.approx(
        0.08, abs=1e-4)
    assert phy045.eta_fit(Ls[1:], [0.9 * L ** -0.08 for L in Ls[1:]],
                          err[1:])["eta"] == pytest.approx(0.08, abs=1e-9)


def test_evaluate_gates_and_verdict_logic():
    a = {n: {"ratio_K_fit_over_area": 1.001} for n in LATTICES}

    def row(R, a_s, z_site, z_area):
        return {"R_measured": R, "R_if_per_area_correct": a_s,
                "z_vs_per_site": z_site, "z_vs_per_area": z_area}
    b = {"square": row(0.99, 1.0, 0.0, 0.0),
         "triangular": row(0.865, 0.866, -30, -0.2),
         "honeycomb": row(1.289, 1.299, 54, -1.8),
         "kagome": row(1.152, 1.155, 22, -0.4)}
    rep = phy045.evaluate({"part_a_harmonic": a, "part_b_mc": b}, None)
    assert rep["verdict"] == "per_area" and rep["overall_pass"]
    b["honeycomb"] = row(1.01, 1.299, 1.0, -40)       # sah aus wie per Site
    rep = phy045.evaluate({"part_a_harmonic": a, "part_b_mc": b}, None)
    assert rep["verdict"] == "undecided"


def test_committed_report_if_present():
    p = ROOT / "results" / "260926 PHY045 helicity normalization O1 report.json"
    if not p.exists():
        pytest.skip("PHY045-Report noch nicht committet")
    rep = json.loads(p.read_text(encoding="utf-8"))
    assert rep["verdict"] == "per_area" and rep["overall_pass"] is True
    for name in ("triangular", "honeycomb", "kagome"):
        assert abs(rep["part_b_mc"][name]["z_vs_per_site"]) > 10
