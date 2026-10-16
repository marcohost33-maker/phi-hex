"""README figure: honeycomb helicity modulus from the committed PHY042 gate report.

Every plotted number is read from
``results/260707 PHY042 honeycomb wl-fss L24-32-48 gate report.json``; nothing
is fitted or smoothed here. Curves are drawn solid inside each L's quotable
domain (``domain_tmax_spread004``: walker spread <= 0.04) and faded beyond it.

    python docs/figures/make_readme_figure.py

writes ``docs/figures/phy042_helicity_honeycomb_{light,dark}.png``.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "results" / "260707 PHY042 honeycomb wl-fss L24-32-48 gate report.json"
OUT = Path(__file__).resolve().parent

# Categorical slots 1-3 of the validated dataviz reference palette
# (validate_palette.js: all checks pass in light and dark mode).
THEMES = {
    "light": {"surface": "#fcfcfb", "ink": "#1f1f1e", "muted": "#6b6a64",
              "grid": "#e4e3dd", "series": ("#2a78d6", "#eb6834", "#1baf7a")},
    "dark": {"surface": "#1a1a19", "ink": "#f0efec", "muted": "#a8a79d",
             "grid": "#34342f", "series": ("#3987e5", "#d95926", "#199e70")},
}
LINESTYLES = ("-", (0, (6, 2)), (0, (2, 1.5)))  # secondary encoding besides color


def load() -> tuple[dict, str]:
    raw = REPORT.read_bytes()
    return json.loads(raw), hashlib.sha256(raw).hexdigest()[:16].upper()


def draw(report: dict, sha16: str, theme: str) -> Path:
    c = THEMES[theme]
    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=150)
    fig.patch.set_facecolor(c["surface"])
    ax.set_facecolor(c["surface"])
    curves = report["main_curves"]
    domains = report["domain_tmax_spread004"]
    t_all = curves[next(iter(curves))]["T"]
    ax.plot(t_all, [2.0 * t / math.pi for t in t_all], color=c["muted"],
            lw=1.2, ls=(0, (1, 2)), zorder=1)
    ax.annotate("Υ = 2T/π  (BKT criterion)", xy=(t_all[-1], 2 * t_all[-1] / math.pi),
                xytext=(-4, 6), textcoords="offset points", ha="right",
                color=c["muted"], fontsize=8.5)
    ref = report["reference_band"]["multi_lattice"]
    ax.axvline(ref, color=c["muted"], lw=0.8, zorder=1)
    ax.annotate(f"reference T_BKT ≈ {ref}", xy=(ref, 1), xycoords=("data", "axes fraction"),
                xytext=(4, -12), textcoords="offset points", color=c["muted"], fontsize=8.5)
    for i, L in enumerate(sorted(curves, key=int)):
        T, y = curves[L]["T"], curves[L]["y2"]
        tmax = domains.get(L)
        color, ls = c["series"][i], LINESTYLES[i]
        inside = [k for k, t in enumerate(T) if tmax is not None and t <= tmax + 1e-12]
        ax.plot(T, y, color=color, lw=2, ls=ls, alpha=0.28, zorder=2)
        if inside:
            Ti = [T[k] for k in inside]
            yi = [y[k] for k in inside]
            ax.plot(Ti, yi, color=color, lw=2, ls=ls, zorder=3, solid_capstyle="round")
            # label at the end of the quotable segment, where the curves separate
            at_edge = Ti[-1] >= t_all[-1] - 1e-9
            ax.annotate(f"L = {L}", xy=(Ti[-1], yi[-1]),
                        xytext=(-4, 9) if at_edge else (7, 2),
                        textcoords="offset points", ha="right" if at_edge else "left",
                        va="center", color=c["ink"], fontsize=9,
                        bbox={"boxstyle": "round,pad=0.15", "fc": c["surface"], "ec": "none"})
    ax.set_xlabel("temperature T", color=c["ink"])
    ax.set_ylabel("helicity modulus Υ", color=c["ink"])
    ax.set_title("Phi-Hex · XY model on the honeycomb lattice (Wang–Landau, PHY042)",
                 color=c["ink"], fontsize=10.5, loc="left")
    ax.grid(True, color=c["grid"], lw=0.6)
    ax.set_axisbelow(True)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(colors=c["muted"], labelsize=8.5)
    ax.set_xlim(t_all[0] - 0.012, t_all[-1] + 0.002)
    fig.text(0.01, 0.01,
             f"Data: results/260707 PHY042 … gate report.json (SHA-256 {sha16}…). "
             "Faded: outside the quotable domain (walker spread > 0.04).",
             color=c["muted"], fontsize=7)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    out = OUT / f"phy042_helicity_honeycomb_{theme}.png"
    fig.savefig(out, facecolor=c["surface"])
    plt.close(fig)
    return out


if __name__ == "__main__":
    data, sha = load()
    for name in THEMES:
        print(draw(data, sha, name))
