"""Regression fuer die sechs Codex-Befunde am PR #75 (Kopf 1e21714,
Review 2026-11-19): 2x P1, 4x P2.

  P1-a  BLAS-Pin erzwungen statt setdefault      (PHY044-WL/045/046/048)
  P1-b  Numba-Backend nur nach Bit-Identitaet    (PHY045 -> 046/048/051)
  P2-a  Jackknife fail-closed bei fehlendem Replikat          (PHY046)
  P2-b  Wall-Budget-Stop S0 in der Produktion + Adjudikation  (PHY046)
  P2-c  v3 = groesstes REGISTRIERTES Paar (und v2 analog)     (PHY046)
  P2-d  PHY031 (honeycomb-Builder) im Provenienz-Fingerprint  (PHY049/050)

Selbstpruefung der eigenen Lieferung (2026-11-21, zweite Runde):

  S-1   Prozess-Pool unter forkserver (Python-3.14-Standard) brach mit
        BrokenProcessPool ab; der parallele Pfad war ungetestet
        (PHY044-WL/045/046/048/049 -> `job_pool` / `map_jobs`)
  S-2   das VAL-BIT-Gate verstellte den globalen numpy-RNG des Aufrufers
  S-3   das Produkt nannte den Mess-Kernel nicht (PHY046 `backend`)

Jeder Guard hat eine Negativprobe, die bei zurueckgenommenem Guard per
Assertion rot wird, und eine Positivkontrolle. Kein Messlauf: die gepinnte
Evidenz in results/ wird nur gelesen; der letzte Abschnitt belegt, dass die
committeten PHY046-/PHY048-Reports unter dem gehaerteten Code unveraendert
regenerieren.
"""
from __future__ import annotations

import ast
import json
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import conftest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"

F44 = "260926 PHY044 honeycomb wl calibration v01.py"
F45 = "260926 PHY045 helicity normalization O1 test v01.py"
F46 = "260926 PHY046 honeycomb w4 wolff area-helicity v01.py"
F48 = "260926 PHY048 square pipeline validation v01.py"
F49 = "260927 PHY049 honeycomb correlation ratio v01.py"
F50 = "260927 PHY050 correlation ratio fss recovery v01.py"
F51 = "260928 PHY051 triangular kagome area-helicity crosscheck v01.py"
F31 = "260607 PHY031 honeycomb tbkt per-site v01.py"

phy045 = conftest._load("phy045_normalization", F45)
phy046 = conftest._load("phy046_w4_v02", F46)
phy048 = conftest._load("phy048_square_validation", F48)
phy049 = conftest._load("phy049_corr_ratio", F49)
phy050 = conftest._load("phy050_corr_ratio_fss", F50)
phy051 = conftest._load("phy051_area_crosscheck", F51)

THREAD_VARS = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
needs_numba = pytest.mark.skipif(not phy045.HAVE_NUMBA,
                                 reason="Numba nicht installiert")


# ============================================================================
# P1-a  BLAS-Pin: erzwungen, vor dem numpy-Import
# ============================================================================

@pytest.mark.parametrize("fname", [F44, F45, F46, F48, F49, F51])
def test_blas_pin_is_an_unconditional_assignment_before_numpy(fname):
    tree = ast.parse((SRC / fname).read_text(encoding="utf-8"))
    loops = [n for n in tree.body if isinstance(n, ast.For)
             and isinstance(n.iter, ast.Tuple)
             and {getattr(e, "value", None) for e in n.iter.elts}
             == set(THREAD_VARS)]
    assert len(loops) == 1
    body = loops[0].body
    assert len(body) == 1 and isinstance(body[0], ast.Assign)
    assert ast.unparse(body[0]) == "os.environ[_v] = '1'"
    first_numpy = min(n.lineno for n in tree.body if isinstance(n, ast.Import)
                      and any(a.name == "numpy" for a in n.names))
    assert loops[0].lineno < first_numpy


def test_no_src_module_pins_blas_threads_with_setdefault():
    """Kein Modul darf den Pin als blosse Vorgabe setzen (Kommentare, die
    `setdefault` nur erwaehnen, zaehlen nicht: geprueft wird der AST)."""
    offenders = []
    for path in sorted(SRC.glob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if (isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Attribute)
                    and node.func.attr == "setdefault"
                    and ast.unparse(node.func.value) == "os.environ"):
                offenders.append(path.name)
    assert offenders == []


@pytest.mark.parametrize("fname", [F44, F45, F46, F48])
def test_inherited_multithread_env_is_overridden_in_a_fresh_process(fname):
    """Verhalten statt Quelltext: ein Runner exportiert 8 Threads; nach dem
    Laden des Moduls in einem frischen Interpreter steht der Pin auf 1, und
    numpy war beim Setzen noch nicht importiert."""
    code = (
        "import importlib.util, json, os, sys\n"
        "assert 'numpy' not in sys.modules\n"
        f"spec = importlib.util.spec_from_file_location('m', {str(SRC / fname)!r})\n"
        "mod = importlib.util.module_from_spec(spec)\n"
        "sys.modules['m'] = mod\n"
        "spec.loader.exec_module(mod)\n"
        f"print(json.dumps({{v: os.environ.get(v) for v in {THREAD_VARS!r}}}))\n"
    )
    env = dict(os.environ, **{v: "8" for v in THREAD_VARS})
    out = subprocess.run([sys.executable, "-c", code], env=env, cwd=ROOT,
                         capture_output=True, text=True, timeout=600)
    assert out.returncode == 0, out.stderr[-2000:]
    assert json.loads(out.stdout.strip().splitlines()[-1]) == {
        v: "1" for v in THREAD_VARS}


# ============================================================================
# P1-b  VAL-BIT: Numba nur nach Bit-Identitaet mit dem Python-Kernel
# ============================================================================

@needs_numba
def test_val_bit_full_kernel_is_bit_identical_here():
    assert phy045.backend_bit_identity() is True


@needs_numba
def test_val_bit_covers_two_coordination_numbers_and_measures():
    degs = set()
    for name, L, _T in phy045.VAL_BIT_CASES:
        _nbr, deg = phy045._nbr_arrays(phy045.build(name, L))
        degs.add(int(deg.max()))
    assert len(degs) >= 2
    n_therm, n_meas = phy045.VAL_BIT_SWEEPS
    assert n_therm >= 1 and n_meas >= 2


@needs_numba
@pytest.mark.parametrize("how", ["one_ulp", "negative_zero", "nan_both",
                                 "shape", "dtype"])
def test_val_bit_comparator_detects_any_deviation(monkeypatch, how):
    """Positivkontrolle des Vergleichers: schon 1 ulp in EINEM Messwert,
    ein Vorzeichen der Null, eine andere Form oder ein anderer dtype macht
    das Gate rot; ein beidseitiges NaN bestaetigt keinen Kernel."""
    real_nb, real_py = phy045._nb_run, phy045._py_run

    def bad_nb(*a):
        out = np.array(real_nb(*a))
        if how == "one_ulp":
            out[-1, 1] = np.nextafter(out[-1, 1], np.inf)
        elif how == "negative_zero":
            out[0, 2] = 0.0
            return out
        elif how == "nan_both":
            out[0, 0] = np.nan
        elif how == "shape":
            return out[:-1]
        elif how == "dtype":
            return out.view(np.int64)    # dieselben Bytes, anderer dtype
        return out

    def py(*a):
        out = np.array(real_py(*a))
        if how == "negative_zero":
            out[0, 2] = -0.0
        elif how == "nan_both":
            out[0, 0] = np.nan
        return out

    monkeypatch.setattr(phy045, "_nb_run", bad_nb)
    monkeypatch.setattr(phy045, "_py_run", py)
    assert phy045.backend_bit_identity() is False


@needs_numba
def test_failed_val_bit_aborts_instead_of_falling_back(monkeypatch):
    monkeypatch.setattr(phy045, "_VAL_BIT_OK", None)
    monkeypatch.setattr(phy045, "backend_bit_identity", lambda: False)
    with pytest.raises(RuntimeError, match="VAL-BIT failed"):
        phy045.production_run_backend()
    # der Fehlschlag ist gecacht: kein zweiter Versuch, kein Python-Fallback
    monkeypatch.setattr(phy045, "backend_bit_identity", lambda: True)
    with pytest.raises(RuntimeError, match="VAL-BIT failed"):
        phy045.production_run_backend()


@needs_numba
def test_successful_val_bit_is_cached_per_interpreter(monkeypatch):
    calls = {"n": 0}

    def ok():
        calls["n"] += 1
        return True

    monkeypatch.setattr(phy045, "_VAL_BIT_OK", None)
    monkeypatch.setattr(phy045, "backend_bit_identity", ok)
    assert phy045.production_run_backend() is phy045._nb_run
    assert phy045.production_run_backend() is phy045._nb_run
    assert calls["n"] == 1


def test_without_numba_the_python_kernel_is_used_without_a_gate(monkeypatch):
    monkeypatch.setattr(phy045, "HAVE_NUMBA", False)
    monkeypatch.setattr(phy045, "_VAL_BIT_OK", None)
    monkeypatch.setattr(phy045, "backend_bit_identity",
                        lambda: pytest.fail("kein Vergleich ohne Numba"))
    assert phy045.production_run_backend() is phy045._py_run


@needs_numba
@pytest.mark.parametrize("label", ["phy045_mc_point", "phy045_efficiency",
                                   "phy046_job", "phy048_job", "phy051_job"])
def test_no_measurement_path_can_bypass_a_failed_val_bit(monkeypatch, label):
    """Jeder Produktionspfad, der den PHY045-Kernel benutzt, haengt am Gate."""
    assert phy046._p45 is phy045 and phy048._p45 is phy045
    assert phy051._p45 is phy045
    monkeypatch.setattr(phy045, "_VAL_BIT_OK", None)
    monkeypatch.setattr(phy045, "backend_bit_identity", lambda: False)
    calls = {
        "phy045_mc_point": lambda: phy045.mc_point(
            ("honeycomb", 4, 0.5, 1, 1, 2)),
        "phy045_efficiency": lambda: phy045.wolff_efficiency_point(
            ("honeycomb", 4, 0.5, 1, 1, 40)),
        "phy046_job": lambda: phy046._job((4, 0, 0.57, 0, 1, 2)),
        "phy048_job": lambda: phy048._job((4, 0, 0.89, 0)),
        "phy051_job": lambda: phy051._job(("kagome", 4, 0, 0.85, 0, 1, 2)),
    }
    with pytest.raises(RuntimeError, match="VAL-BIT failed"):
        calls[label]()


def test_numba_kernel_is_only_handed_out_by_the_gate():
    """Quelltext-Schranke: `_nb_run` wird nirgends mehr direkt aus der
    blossen Verfuegbarkeit gewaehlt; der einzige Leser ausserhalb der
    Definition ist `production_run_backend` / `backend_bit_identity`."""
    for path in sorted(SRC.glob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.IfExp):
                txt = ast.unparse(node)
                assert not ("_nb_run" in txt and "HAVE_NUMBA" in txt
                            and "_nb_run_corr" not in txt), (path.name, txt)
    tree = ast.parse((SRC / F45).read_text(encoding="utf-8"))
    readers = set()
    for fn in ast.walk(tree):
        if isinstance(fn, ast.FunctionDef):
            for node in ast.walk(fn):
                if (isinstance(node, ast.Name) and node.id == "_nb_run"
                        and isinstance(node.ctx, ast.Load)):
                    readers.add(fn.name)
    assert readers == {"backend_bit_identity", "production_run_backend"}


# ============================================================================
# P2-a  PHY046 Jackknife fail-closed
# ============================================================================

@pytest.mark.parametrize("bad", [None, math.nan, math.inf, -math.inf])
@pytest.mark.parametrize("pos", [0, 3, 7])
def test_jk_se_rejects_an_incomplete_replica_set(bad, pos):
    vals = [0.57] * 8
    vals[pos] = bad
    assert phy046._jk_se(vals) is None


def test_jk_se_positive_control_matches_the_jackknife_formula():
    vals = [0.570, 0.571, 0.569, 0.572, 0.568, 0.570, 0.571, 0.569]
    n = len(vals)
    m = sum(vals) / n
    want = math.sqrt((n - 1) / n * sum((x - m) ** 2 for x in vals))
    assert phy046._jk_se(vals) == pytest.approx(want, rel=1e-12)
    assert phy046._jk_se([0.57]) is None


def _wm_curves(T0, C, t_grid, ladder, slope=4.0):
    return {L: np.array([(2 * T0 / math.pi) * (1 + 1 / (2 * math.log(L) + C))
                         - slope * (T - T0) for T in t_grid]) for L in ladder}


def _synthetic_prod(ladder, T0=0.575, noise=1e-4, seed=3, n_seeds=4):
    t = list(phy046.W4V2_T_GRID)
    rng = np.random.default_rng(seed)
    base = _wm_curves(T0, 2.0, t, ladder)
    area = {str(L): (base[L][:, None] + noise * rng.standard_normal(
        (len(t), n_seeds))).tolist() for L in ladder}
    return {"t_grid": t, "ladder": list(ladder), "n_seeds": n_seeds,
            "ups_area": area, "ups_site": area}


def test_one_failed_leave_one_out_fit_ends_in_s2_not_in_a_verdict(monkeypatch):
    """Codex-Beispiel auf Analyse-Ebene: faellt der Primaer-Fit in EINEM
    Leave-one-seed-out-Replikat aus, gibt es kein sigma_sampler und damit
    Regel S2 - vorher wurde die SE aus den Ueberlebenden gerechnet und ein
    C/I-Verdikt ausgesprochen."""
    prod = _synthetic_prod(phy046.W4V2_LADDER)
    base = phy046.analyse(prod, "ups_area")
    assert base["verdict"]["rule"] in ("C", "I")     # Positivkontrolle
    real = phy046.estimate
    n = {"calls": 0}

    def flaky(*a, **k):
        out = real(*a, **k)
        n["calls"] += 1
        if n["calls"] == 3:          # 1 = Vollsatz, 2.. = Jackknife-Replikate
            out = dict(out, primary_hks=None)
        return out

    monkeypatch.setattr(phy046, "estimate", flaky)
    out = phy046.analyse(prod, "ups_area")
    assert out["T_w4"] == base["T_w4"]
    assert out["sigma_sampler"] is None
    assert out["verdict"] == {"verdict": "NEGATIVE_RESULT", "rule": "S2",
                              "sigma_tot": None}


# ============================================================================
# P2-c  PHY046 v2/v3 an die REGISTRIERTEN Paare gebunden
# ============================================================================

_PAIRS = [(32, 64), (48, 96), (64, 128), (96, 192), (128, 256)]


def _estimate_with(monkeypatch, crossings: dict):
    ladder = list(phy046.W4V2_LADDER)
    t = list(phy046.W4V2_T_GRID)
    monkeypatch.setattr(phy046, "pair_crossings",
                        lambda c, tg, lad: dict(crossings))
    curves = _wm_curves(0.575, 2.0, t, ladder)
    se = {L: np.full(len(t), 1e-3) for L in ladder}
    return phy046.estimate(curves, se, t, ladder)


def test_v3_is_none_when_the_largest_registered_pair_does_not_cross(
        monkeypatch):
    pc = dict(zip(_PAIRS, [0.5660, 0.5655, 0.5657, 0.5649, None]))
    out = _estimate_with(monkeypatch, pc)
    assert out["n_pairs_with_crossing"] == 4
    assert out["v3_largest_pair"] is None            # nicht 0.5649 = (96,192)


def test_v3_positive_control_is_the_largest_registered_pair(monkeypatch):
    pc = dict(zip(_PAIRS, [0.5660, None, 0.5657, 0.5649, 0.5663]))
    out = _estimate_with(monkeypatch, pc)
    assert out["v3_largest_pair"] == 0.5663


def test_v2_drops_only_the_smallest_registered_pair(monkeypatch):
    """Kreuzt (32,64) nicht, darf v2 nicht stattdessen (48,96) verwerfen:
    'ohne das kleinste Paar' ist dann der Fit ueber alle Ueberlebenden."""
    pc = dict(zip(_PAIRS, [None, 0.5655, 0.5657, 0.5649, 0.5663]))
    out = _estimate_with(monkeypatch, pc)
    want = phy046.fit_hks([48, 64, 96, 128],
                          [0.5655, 0.5657, 0.5649, 0.5663])
    assert out["v2_hks_drop_smallest"] == want
    assert out["v2_hks_drop_smallest"] == out["primary_hks"]
    # Positivkontrolle: mit allen fuenf Paaren faellt genau (32,64) weg
    pc5 = dict(zip(_PAIRS, [0.5660, 0.5655, 0.5657, 0.5649, 0.5663]))
    out5 = _estimate_with(monkeypatch, pc5)
    assert out5["v2_hks_drop_smallest"] == want
    assert out5["v2_hks_drop_smallest"] != out5["primary_hks"]


def test_missing_largest_pair_is_not_counted_in_sigma_fss(monkeypatch):
    """Auf Analyse-Ebene: Variante v3 ist None und fehlt in der Spannweite."""
    prod = _synthetic_prod(phy046.W4V2_LADDER)
    real = phy046.pair_crossings

    def no_largest(curves, t_grid, ladder):
        out = real(curves, t_grid, ladder)
        out[(128, 256)] = None
        return out

    monkeypatch.setattr(phy046, "pair_crossings", no_largest)
    out = phy046.analyse(prod, "ups_area")
    assert out["estimates"]["pairs"]["128_256"] is None
    assert out["variants"][3] is None
    vv = [v for v in out["variants"] if v is not None]
    assert out["sigma_fss_raw"] == pytest.approx(0.5 * (max(vv) - min(vv)))


# ============================================================================
# P2-b  PHY046 Wall-Budget-Stop S0
# ============================================================================

def _fake_job(args):
    L, t_idx, T, s, _n_therm, _n_meas = args
    v = 0.5 + 1e-3 * L + 1e-2 * t_idx + 1e-4 * s
    return {"L": L, "t_idx": t_idx, "s": s, "ups_area": v,
            "ups_site": 1.299 * v, "wall_s": 0.01}


def _hour_clock():
    """Jeder Aufruf ist eine Stunde spaeter (0, 1, 2, ... h in Sekunden)."""
    state = {"k": -1}

    def clock():
        state["k"] += 1
        return 3600.0 * state["k"]
    return clock


def test_budget_stop_marks_remaining_sizes_unmeasured(monkeypatch):
    seen = []
    monkeypatch.setattr(phy046, "_job",
                        lambda a: (seen.append(a[0]), _fake_job(a))[1])
    prod = phy046.produce(ladder=(4, 8, 16), t_grid=(0.55, 0.56),
                          n_seeds=2, max_workers=1, wall_budget_h=0.5,
                          _clock=_hour_clock())
    assert set(seen) == {16}                 # groesstes L zuerst, dann Stop
    assert prod["unmeasured_L"] == [4, 8]
    assert prod["wall_budget_h"] == 0.5
    assert np.all(np.isfinite(prod["ups_area"]["16"]))
    assert np.all(np.isnan(prod["ups_area"]["8"]))
    assert np.all(np.isnan(prod["ups_site"]["4"]))


def test_within_budget_every_size_is_measured_and_numbers_are_unchanged(
        monkeypatch):
    monkeypatch.setattr(phy046, "_job", _fake_job)
    prod = phy046.produce(ladder=(4, 8, 16), t_grid=(0.55, 0.56),
                          n_seeds=2, max_workers=1, wall_budget_h=1e9,
                          _clock=_hour_clock())
    assert prod["unmeasured_L"] == []
    for L in (4, 8, 16):
        for k in range(2):
            for s in range(2):
                assert prod["ups_area"][str(L)][k][s] == _fake_job(
                    (L, k, 0.0, s, 0, 0))["ups_area"]


def test_zero_budget_still_measures_the_first_block(monkeypatch):
    monkeypatch.setattr(phy046, "_job", _fake_job)
    prod = phy046.produce(ladder=(4, 8), t_grid=(0.55,), n_seeds=2,
                          max_workers=1, wall_budget_h=0.0,
                          _clock=_hour_clock())
    assert prod["unmeasured_L"] == [4]


@pytest.mark.parametrize("bad", [math.nan, math.inf, -1.0, True, "24", None])
def test_invalid_wall_budget_is_rejected_before_any_job(monkeypatch, bad):
    monkeypatch.setattr(phy046, "_job",
                        lambda a: pytest.fail("Job trotz ungueltigem Budget"))
    with pytest.raises(ValueError, match="wall_budget_h"):
        phy046.produce(ladder=(4, 8), t_grid=(0.55,), n_seeds=2,
                       max_workers=1, wall_budget_h=bad)


def test_default_budget_is_the_preregistered_24_hours():
    import inspect
    sig = inspect.signature(phy046.produce)
    assert sig.parameters["wall_budget_h"].default == 24.0
    assert phy046.W4V2_WALL_BUDGET_H == 24.0
    spec = (ROOT / "spec" / "260926 PHI HEX w4 honeycomb preregistration "
            "v02.md").read_text(encoding="utf-8")
    assert "24 h Wall auf 4 Kernen" in spec


def _stopped_prod(drop=(32,)):
    prod = _synthetic_prod(phy046.W4V2_LADDER)
    for L in drop:
        for key in ("ups_area", "ups_site"):
            prod[key][str(L)] = np.full(
                (len(prod["t_grid"]), prod["n_seeds"]), np.nan).tolist()
    prod["unmeasured_L"] = sorted(drop)
    prod["wall_budget_h"] = 24.0
    return prod


def test_stopped_run_gets_s0_even_if_measured_pairs_would_allow_a_verdict():
    """Ohne L=32 bleiben vier Paare mit Crossing (>= 3): ohne S0 gaebe es
    ein gewoehnliches C/I-Verdikt aus einer abgebrochenen Leiter."""
    complete = phy046.analyse(_synthetic_prod(phy046.W4V2_LADDER), "ups_area")
    assert complete["ladder_complete"] is True
    assert complete["verdict"]["rule"] in ("C", "I")  # Positivkontrolle
    out = phy046.analyse(_stopped_prod(), "ups_area")
    assert out["estimates"]["n_pairs_with_crossing"] == 4
    assert out["estimates"]["pairs"]["32_64"] is None
    assert out["ladder_used"] == [48, 64, 96, 128, 192, 256]
    assert out["ladder_complete"] is False
    assert out["verdict"] == {"verdict": "NEGATIVE_RESULT", "rule": "S0",
                              "sigma_tot": None}


def test_s0_cannot_be_bypassed_by_dropping_the_declaration():
    """Fail-closed: NaN-Daten ohne `unmeasured_L` sind trotzdem S0."""
    prod = _stopped_prod()
    del prod["unmeasured_L"]
    assert phy046.analyse(prod, "ups_area")["verdict"]["rule"] == "S0"
    # ... und ein einzelnes fehlendes (T, Seed)-Feld genuegt
    prod = _synthetic_prod(phy046.W4V2_LADDER)
    prod["ups_area"]["96"][5][1] = math.nan
    out = phy046.analyse(prod, "ups_area")
    assert out["verdict"]["rule"] == "S0" and 96 not in out["ladder_used"]


def test_stopped_run_survives_the_json_round_trip_and_is_reported(tmp_path):
    rep = phy046._clean(phy046.build_report(
        dict(_stopped_prod(), n_therm=500, n_meas=1500, wall_s=90000.0,
             cpu_s=1.0, max_workers=4)))
    text = json.dumps(rep, allow_nan=False)          # NaN -> null, kein Crash
    back = json.loads(text)
    assert back["production"]["unmeasured_L"] == [32]
    assert back["primary_per_area"]["verdict"]["rule"] == "S0"
    assert back["convention_crosscheck_per_site"]["verdict"]["rule"] == "S0"
    # Re-Analyse aus dem gelesenen Report (null statt NaN) bleibt S0
    prod = dict(back["production"], **back["data"])
    assert phy046.analyse(prod, "ups_area")["verdict"]["rule"] == "S0"
    phy046.write_text_report(back, tmp_path / "r.txt")
    txt = (tmp_path / "r.txt").read_text(encoding="utf-8")
    assert "S0 BUDGET-STOP" in txt and "L = [32]" in txt
    assert "'rule': 'S0'" in txt


# ============================================================================
# P2-d  PHY031 im Provenienz-Fingerprint (Produzent + Adjudikation)
# ============================================================================

def test_honeycomb_builder_is_fingerprinted_by_producer_and_adjudicator():
    rel = f"src/{F31}"
    assert rel in phy049._RUNTIME_SOURCE_PATHS
    assert phy049._RUNTIME_SOURCE_PATHS == phy050._RUNTIME_SOURCE_PATHS
    assert all((ROOT / p).is_file() for p in phy049._RUNTIME_SOURCE_PATHS)
    # der Fingerprint deckt genau die Datei, an die PHY045.build delegiert
    assert Path(phy045._p31.__file__).resolve() == (ROOT / rel).resolve()
    prov = phy049._runtime_provenance()["source_sha256"]
    assert prov == phy050._current_runtime_file_hashes()
    assert len(prov[rel]) == 64


def _valid_provenance() -> dict:
    prov = phy049._runtime_provenance()
    assert phy050._valid_runtime_provenance(prov) is True   # Positivkontrolle
    return json.loads(json.dumps(prov))


def test_adjudication_rejects_a_changed_honeycomb_builder():
    prov = _valid_provenance()
    prov["source_sha256"][f"src/{F31}"] = "0" * 64
    assert phy050._valid_runtime_provenance(prov) is False


def test_adjudication_rejects_a_fingerprint_without_the_builder():
    """Ein Checkpoint mit dem alten Fuenf-Dateien-Fingerprint ist nicht mehr
    adjudizierbar (es gibt keinen solchen Produktions-Checkpoint)."""
    prov = _valid_provenance()
    del prov["source_sha256"][f"src/{F31}"]
    assert phy050._valid_runtime_provenance(prov) is False


# ============================================================================
# Gepinnte Evidenz regeneriert unter dem gehaerteten Code unveraendert
# ============================================================================

def _diff(a, b, path="") -> list:
    """Jeder Schluessel/Wert aus `a` muss in `b` stehen (Floats rel 1e-9)."""
    if isinstance(a, dict):
        if not isinstance(b, dict):
            return [f"{path}: Typ"]
        out = []
        for k, v in a.items():
            out += [f"{path}/{k}: fehlt"] if k not in b else _diff(
                v, b[k], f"{path}/{k}")
        return out
    if isinstance(a, list):
        if not isinstance(b, list) or len(a) != len(b):
            return [f"{path}: Laenge"]
        return [d for i, (x, y) in enumerate(zip(a, b))
                for d in _diff(x, y, f"{path}[{i}]")]
    if isinstance(a, float) and isinstance(b, (int, float)) \
            and not isinstance(b, bool):
        return [] if math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12) \
            else [f"{path}: {a!r} != {b!r}"]
    return [] if a == b and type(a) is type(b) else [f"{path}: {a!r} != {b!r}"]


def test_diff_helper_positive_control():
    a = {"x": [1.0, {"y": 0.5648514935621091}], "v": None}
    assert _diff(a, json.loads(json.dumps(a))) == []
    assert _diff(a, {"x": [1.0, {"y": 0.5648514935621091}]}) != []
    assert _diff(a, {"x": [1.0, {"y": 0.56485}], "v": None}) != []
    assert _diff(a, {"x": [1.0, {"y": None}], "v": None}) != []
    assert _diff({"n": 5}, {"n": 5.0}) != []


def test_committed_phy046_report_regenerates_unchanged(tmp_path):
    stem = ROOT / "results" / phy046.REPORT_STEM
    rep = json.loads(stem.with_suffix(".json").read_text(encoding="utf-8"))
    prod = dict(rep["production"], **rep["data"])
    new = phy046._clean(phy046.build_report(prod))
    assert _diff(rep, new) == []
    # die einzigen neuen Felder sind die S0-/Leiter-Angaben
    assert set(new["production"]) - set(rep["production"]) == {
        "unmeasured_L", "wall_budget_h", "backend"}
    # der Lauf vom 2026-09-26 hat seinen Kernel nicht protokolliert: das
    # Feld bleibt leer, statt nachtraeglich "numba" zu behaupten
    assert new["production"]["backend"] is None
    for block in ("primary_per_area", "convention_crosscheck_per_site"):
        assert set(new[block]) - set(rep[block]) == {
            "ladder_used", "ladder_complete"}
        assert new[block]["ladder_complete"] is True
    assert new["production"]["unmeasured_L"] == []
    assert new["primary_per_area"]["verdict"]["rule"] == "C"
    assert new["convention_crosscheck_per_site"]["verdict"]["rule"] == "S2"
    phy046.write_text_report(new, tmp_path / "r.txt")
    assert (tmp_path / "r.txt").read_bytes() == stem.with_suffix(
        ".txt").read_bytes()


def test_committed_phy048_analysis_regenerates_unchanged():
    rep = json.loads((ROOT / "results" / f"{phy048.STEM}.json").read_text(
        encoding="utf-8"))
    prod = dict(rep["production"], **rep["data"])
    new = phy046._clean(phy046.analyse(prod, "ups_area"))
    assert _diff(rep["analysis"], new) == []
    assert set(new) - set(rep["analysis"]) == {"ladder_used",
                                               "ladder_complete"}


# ============================================================================
# S-1  Prozess-Pool: unabhaengig von der Start-Methode, parallel == sequenziell
# ============================================================================

class _Recorder:
    """Ersetzt ProcessPoolExecutor: merkt sich die Argumente, startet nichts."""
    calls: list = []

    def __init__(self, **kw):
        type(self).calls.append(kw)


@pytest.fixture
def recorder(monkeypatch):
    _Recorder.calls = []
    monkeypatch.setattr(phy045, "ProcessPoolExecutor", _Recorder)
    return _Recorder


def _a_job(x):
    return x


def test_job_pool_is_sequential_for_a_single_worker(recorder):
    assert phy045.job_pool(_a_job, 1) is None
    assert phy045.job_pool(_a_job, 0) is None
    assert recorder.calls == []


def test_job_pool_requests_fork_explicitly_where_it_exists(recorder,
                                                           monkeypatch):
    """Nicht den Standard-Kontext erben: der ist seit 3.14 forkserver."""
    import multiprocessing
    monkeypatch.setattr(multiprocessing, "get_all_start_methods",
                        lambda: ["forkserver", "fork", "spawn"])
    sentinel = object()
    monkeypatch.setattr(multiprocessing, "get_context",
                        lambda m: (m, sentinel))
    assert isinstance(phy045.job_pool(_a_job, 3), _Recorder)
    assert recorder.calls == [{"max_workers": 3,
                               "mp_context": ("fork", sentinel)}]


def test_job_pool_without_fork_depends_on_where_the_job_lives(
        recorder, monkeypatch, capsys):
    import multiprocessing
    monkeypatch.setattr(multiprocessing, "get_all_start_methods",
                        lambda: ["spawn"])
    # per Kunstnamen geladen: ein spawn-Kind koennte den Job nicht
    # importieren -> sequenziell, mit Hinweis
    assert _a_job.__module__ != "__main__"
    assert phy045.job_pool(_a_job, 4) is None
    assert recorder.calls == []
    assert "sequenziell" in capsys.readouterr().out

    def main_job(x):
        return x
    main_job.__module__ = "__main__"       # als Skript gestartet
    assert isinstance(phy045.job_pool(main_job, 4), _Recorder)
    assert recorder.calls == [{"max_workers": 4}]


def test_map_jobs_keeps_input_order_sequentially():
    assert phy045.map_jobs(_a_job, [3, 1, 2], 1) == [3, 1, 2]


_POOL_SCRIPT = r"""
import json, multiprocessing, sys
import numpy as np
if __name__ == "__main__":
    # der Python-3.14-Standard, auf jeder Version erzwungen: so reproduziert
    # der Test den BrokenProcessPool auch unter 3.12/3.13
    multiprocessing.set_start_method("forkserver", force=True)
    sys.path.insert(0, "tests")
    import conftest
    p45 = conftest._load("phy045_normalization",
                         "260926 PHY045 helicity normalization O1 test v01.py")
    p46 = conftest._load("phy046_w4_v02",
                         "260926 PHY046 honeycomb w4 wolff area-helicity v01.py")
    p48 = conftest._load("phy048_square_validation",
                         "260926 PHY048 square pipeline validation v01.py")
    p49 = conftest._load("phy049_corr_ratio",
                         "260927 PHY049 honeycomb correlation ratio v01.py")
    p44 = conftest._load("phy044_wl_calibration",
                         "260926 PHY044 honeycomb wl calibration v01.py")
    out = {"default_start_method": multiprocessing.get_start_method()}

    kw = dict(ladder=(4, 8), t_grid=(0.55, 0.60), n_seeds=2, n_therm=3,
              n_meas=6)
    a, b = p46.produce(max_workers=2, **kw), p46.produce(max_workers=1, **kw)
    out["phy046"] = all(a[k] == b[k] for k in ("ups_area", "ups_site"))
    out["phy046_finite"] = bool(np.all(np.isfinite(a["ups_area"]["8"])))

    p48.LADDER, p48.T_GRID = (4, 8), (0.88, 0.90)
    p48.N_SEEDS, p48.N_THERM, p48.N_MEAS = 2, 3, 6
    out["phy048"] = (p48.produce(max_workers=2)["ups_area"]
                     == p48.produce(max_workers=1)["ups_area"])

    plan = {"square": (0.45, (4, 6, 8, 10))}
    strip = lambda r: json.dumps(r, sort_keys=True)
    pb = lambda w: p45.part_b(n_seeds=2, n_therm=2, n_meas=40,
                              max_workers=w, plan=plan)
    out["phy045_part_b"] = strip(pb(2)) == strip(pb(1))

    jobs = [(8, k, T, s, 2, 3) for k, T in enumerate((0.56, 0.58))
            for s in range(2)]
    drop = lambda rows: [{k: v for k, v in r.items() if k != "wall_s"}
                         for r in rows]
    out["phy049"] = (drop(p45.map_jobs(p49._job, jobs, 2))
                     == drop([p49._job(j) for j in jobs]))

    def fake(args):
        return {"echo": list(args)}
    p44._job = fake
    j44 = [("x", 4, 0), ("x", 8, 1), ("x", 8, 0)]
    order = sorted(j44, key=lambda j: (-j[1], j[2]))
    out["phy044_wl"] = (p44._run_jobs(j44, 2)
                        == [{"echo": list(j)} for j in order])
    print("RESULT " + json.dumps(out))
"""


def test_parallel_production_is_bit_equal_to_sequential_under_forkserver(
        tmp_path):
    """Ende-zu-Ende in einem frischen Interpreter mit forkserver als
    Standard: jeder parallele Produktionspfad laeuft durch und liefert
    exakt die sequenziellen Zahlen (vorher: BrokenProcessPool)."""
    script = tmp_path / "pool_check.py"
    script.write_text(_POOL_SCRIPT, encoding="utf-8")
    out = subprocess.run([sys.executable, str(script)], cwd=ROOT,
                         capture_output=True, text=True, timeout=900)
    assert out.returncode == 0, out.stderr[-3000:]
    line = [ln for ln in out.stdout.splitlines() if ln.startswith("RESULT ")]
    assert len(line) == 1, out.stdout[-2000:]
    assert json.loads(line[0][7:]) == {
        "default_start_method": "forkserver", "phy046": True,
        "phy046_finite": True, "phy048": True, "phy045_part_b": True,
        "phy049": True, "phy044_wl": True}


def test_w4_modules_never_build_a_pool_with_the_inherited_default_context():
    """Quelltext-Schranke: PHY046/048/049 bauen keinen eigenen Pool mehr;
    PHY045 und PHY044-WL nur ueber den fork-Zweig oder den `__main__`-Zweig."""
    for fname in (F46, F48, F49):
        assert "ProcessPoolExecutor" not in (SRC / fname).read_text(
            encoding="utf-8"), fname
    for fname in (F45, F44):
        tree = ast.parse((SRC / fname).read_text(encoding="utf-8"))
        calls = [n for n in ast.walk(tree) if isinstance(n, ast.Call)
                 and getattr(n.func, "id", None) == "ProcessPoolExecutor"]
        assert len(calls) == 2, fname
        with_ctx = [c for c in calls
                    if any(k.arg == "mp_context" for k in c.keywords)]
        assert len(with_ctx) == 1, fname
        assert "get_context('fork')" in ast.unparse(with_ctx[0])


# ============================================================================
# S-2  VAL-BIT veraendert den globalen numpy-RNG des Aufrufers nicht
# ============================================================================

@needs_numba
def test_val_bit_leaves_the_callers_global_rng_untouched():
    np.random.seed(20261121)
    np.random.random(7)
    before = np.random.get_state()
    want = np.random.RandomState()
    want.set_state(before)
    assert phy045.backend_bit_identity() is True
    after = np.random.get_state()
    assert before[0] == after[0] and before[2:] == after[2:]
    assert np.array_equal(before[1], after[1])
    assert np.random.random() == want.random_sample()


@needs_numba
def test_val_bit_restores_the_rng_even_when_it_fails(monkeypatch):
    real = phy045._nb_run
    monkeypatch.setattr(phy045, "_nb_run", lambda *a: np.array(real(*a))[:-1])
    np.random.seed(5)
    before = np.random.get_state()
    assert phy045.backend_bit_identity() is False
    assert np.array_equal(before[1], np.random.get_state()[1])
    assert before[2] == np.random.get_state()[2]


# ============================================================================
# S-3  Das Produkt nennt den Kernel, mit dem gemessen wurde
# ============================================================================

def test_production_names_the_kernel_that_measured(monkeypatch):
    monkeypatch.setattr(phy046, "_job", _fake_job)
    prod = phy046.produce(ladder=(4, 8), t_grid=(0.55,), n_seeds=2,
                          max_workers=1)
    assert prod["backend"] == ("numba" if phy045.HAVE_NUMBA else "python")
    rep = phy046.build_report(dict(prod, n_therm=1, n_meas=1))
    assert rep["production"]["backend"] == prod["backend"]


def test_backend_name_follows_the_gate_not_the_mere_availability(monkeypatch):
    monkeypatch.setattr(phy045, "HAVE_NUMBA", False)
    monkeypatch.setattr(phy045, "_VAL_BIT_OK", None)
    assert phy045.production_backend_name() == "python"


@needs_numba
def test_production_stops_before_any_job_when_val_bit_fails(monkeypatch):
    monkeypatch.setattr(phy045, "_VAL_BIT_OK", None)
    monkeypatch.setattr(phy045, "backend_bit_identity", lambda: False)
    monkeypatch.setattr(phy046, "_job",
                        lambda a: pytest.fail("Job trotz VAL-BIT-Fehlschlag"))
    with pytest.raises(RuntimeError, match="VAL-BIT failed"):
        phy046.produce(ladder=(4, 8), t_grid=(0.55,), n_seeds=2,
                       max_workers=1)
