"""Pytest-Fixtures + Modul-Loader fuer Phi-Hex.

Die Engine-Dateien tragen historisch Leerzeichen im Namen
("260602 PHI HEX core v2 2 hardened.py") und die PHY0xx-Module importieren
sie unter sauberen Modulnamen (phi_hex_core_v2, phy026_wolff_cluster, ...)
ueber einen Dev-Pfad /home/claude. Dieser Loader registriert die Engine-
Dateien unter genau diesen sauberen Namen in sys.modules, sodass sowohl die
Tests als auch die bestehenden PHY0xx-Imports plattformunabhaengig aufloesen.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parent.parent / "src"

# Logischer Name -> Dateiname (mit Leerzeichen). Reihenfolge = Ladereihenfolge:
# Kern zuerst, dann Module, die ihn importieren.
_MODULE_FILES = {
    "phi_hex_core_v2": "260602 PHI HEX core v2 2 hardened.py",
    "phy026_wolff_cluster": "260602 PHY026 wolff cluster v01.py",
    "phy028_square_validation": "260602 PHY028 square validation v01.py",
    "phy030_triangular_tbkt": "260604 PHY030 triangular tbkt per-site v01.py",
    "phy030_tbkt_wm_logfit": "260604 PHY030 triangular tbkt per-site v02 wm-logfit.py",
    "phy031_honeycomb": "260607 PHY031 honeycomb tbkt per-site v01.py",
    "phy032_honeycomb_wm_bootstrap":
        "260607 PHY032 honeycomb wm-logfit bootstrap v01.py",
    "phy033_kagome": "260609 PHY033 kagome tbkt per-site v01.py",
    "phy034_hks_extrapolation":
        "260611 PHY034 honeycomb hks extrapolation v01.py",
    "phy035_hks_dense_ladder":
        "260611 PHY035 honeycomb hks dense ladder v01.py",
    "phy036_hks_large_ladder":
        "260611 PHY036 honeycomb hks large ladder v01.py",
    "phy037_hks_multilattice":
        "260611 PHY037 hks multilattice extrapolation v01.py",
    "phy038_hks_highstat":
        "260612 PHY038 honeycomb hks highstat v01.py",
    # Coverage-Schliessung 2026-07-10 (Code-Audit): diese fuenf Module hatten
    # bis dahin KEIN Test-Gate (nur compile+lint). Import + Orakel-Smokes in
    # tests/test_module_smoke_gates.py.
    "phy017_validation_kernel": "phy017 validation kernel v1 3 1.py",
    "phy027_precision_fss": "260602 PHY027 precision fss v01.py",
    "phy029_triangular_sandvik": "260602 PHY029 triangular sandvik v01.py",
    "phy024_r0_reference": str(Path("phy024")
                               / "PHY024_R0_XY_PhiHex_Reference_v01.py"),
    "phy024_r2_honeycomb": str(Path("phy024")
                               / "PHY024_R2_PeriodicHoneycomb_XY_PhiHex_v01.py"),
}


def _load(name: str, filename: str):
    if name in sys.modules:
        return sys.modules[name]
    path = SRC / filename
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    # VOR exec registrieren (dataclass-Resolver + interne Imports brauchen das).
    sys.modules[name] = mod
    try:
        spec.loader.exec_module(mod)
    except BaseException:
        # Haertung 2026-07-10 (Code-Audit): halb-initialisiertes Modul nicht
        # in sys.modules zuruecklassen - ein Retry bekaeme sonst still das
        # kaputte Modul statt eines erneuten Fehlers.
        sys.modules.pop(name, None)
        raise
    return mod


def load_all():
    for name, filename in _MODULE_FILES.items():
        _load(name, filename)
    return {name: sys.modules[name] for name in _MODULE_FILES}


# Beim Import von conftest sofort registrieren, damit auch die internen
# `sys.path.insert(0, "/home/claude"); from phi_hex_core_v2 import ...`-Zeilen
# der PHY0xx-Module greifen (sie finden den Namen dann bereits in sys.modules).
load_all()
