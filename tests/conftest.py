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
}


def _load(name: str, filename: str):
    if name in sys.modules:
        return sys.modules[name]
    path = SRC / filename
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    # VOR exec registrieren (dataclass-Resolver + interne Imports brauchen das).
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def load_all():
    for name, filename in _MODULE_FILES.items():
        _load(name, filename)
    return {name: sys.modules[name] for name in _MODULE_FILES}


# Beim Import von conftest sofort registrieren, damit auch die internen
# `sys.path.insert(0, "/home/claude"); from phi_hex_core_v2 import ...`-Zeilen
# der PHY0xx-Module greifen (sie finden den Namen dann bereits in sys.modules).
load_all()
