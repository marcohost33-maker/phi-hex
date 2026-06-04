#!/usr/bin/env python3
"""
PHY024-R0 XY-Phi-Hex Reference Implementation v0.1
Coworker Research / Coworkerz — 2026-06-01

Purpose
-------
This script is a minimal, reproducible R0 reference for testing whether explicit
XY phase coupling can repair the PHY023 negative result in Phi-Hex v1.x.

It does NOT prove BKT. It only computes finite-size, BKT-compatible indicators:
- g1(r) phase correlation by graph distance
- algebraic vs exponential decay fit comparison
- vortex/antivortex winding on hex-ring plaquettes
- vortex pairing fraction
- helicity analog
- bootstrap confidence intervals
- Benjamini-Hochberg FDR helper
- optional degree-preserving edge shuffle null model

Dependencies
------------
Python >= 3.10, NumPy.

Example
-------
python PHY024_R0_XY_PhiHex_Reference_v01.py --radius 3 --jtheta 0.4 --sigma 0.05 --seeds 16 --steps 800 --burn-in 300
python PHY024_R0_XY_PhiHex_Reference_v01.py --radius 3 --jtheta 0.4 --sigma 0.05 --seeds 16 --null shuffled
"""

from __future__ import annotations

import argparse
import json
import math
from collections import deque
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np


Axial = Tuple[int, int]
Edge = Tuple[int, int]


DIRECTIONS: List[Axial] = [
    (1, 0),
    (1, -1),
    (0, -1),
    (-1, 0),
    (-1, 1),
    (0, 1),
]


def wrap_angle(x: np.ndarray | float) -> np.ndarray | float:
    """Wrap angle(s) to (-pi, pi]."""
    return (x + math.pi) % (2.0 * math.pi) - math.pi


@dataclass(frozen=True)
class HexPatch:
    radius: int
    coords: List[Axial]
    index: Dict[Axial, int]
    edges: List[Edge]
    neighbors: List[List[int]]
    plaquette_rings: List[List[int]]
    dist_matrix: np.ndarray


def axial_distance(a: Axial, b: Axial) -> int:
    """Hex-grid distance between two axial coordinates."""
    aq, ar = a
    bq, br = b
    return int((abs(aq - bq) + abs(ar - br) + abs((aq + ar) - (bq + br))) // 2)


def build_hex_patch(radius: int) -> HexPatch:
    """
    Build a finite axial hex patch.

    Radius convention:
    radius=0 -> 1 node
    radius=1 -> 7 nodes
    radius=2 -> 19 nodes
    radius=3 -> 37 nodes
    """
    if radius < 1:
        raise ValueError("radius must be >= 1")

    coords: List[Axial] = []
    for q in range(-radius, radius + 1):
        for r in range(-radius, radius + 1):
            if max(abs(q), abs(r), abs(q + r)) <= radius:
                coords.append((q, r))

    coords.sort()
    index = {coord: i for i, coord in enumerate(coords)}
    n = len(coords)

    edge_set: set[Edge] = set()
    neighbors: List[List[int]] = [[] for _ in range(n)]

    for i, (q, r) in enumerate(coords):
        for dq, dr in DIRECTIONS:
            nb = (q + dq, r + dr)
            j = index.get(nb)
            if j is not None:
                a, b = sorted((i, j))
                edge_set.add((a, b))
                neighbors[i].append(j)

    edges = sorted(edge_set)

    # Plaquette ring: for each interior center with all six neighbor cells present,
    # use the six neighbors around the center as a closed ring.
    plaquette_rings: List[List[int]] = []
    for q, r in coords:
        ring: List[int] = []
        ok = True
        for dq, dr in DIRECTIONS:
            j = index.get((q + dq, r + dr))
            if j is None:
                ok = False
                break
            ring.append(j)
        if ok:
            plaquette_rings.append(ring)

    dist_matrix = np.zeros((n, n), dtype=np.int16)
    for i, a in enumerate(coords):
        for j, b in enumerate(coords):
            dist_matrix[i, j] = axial_distance(a, b)

    return HexPatch(
        radius=radius,
        coords=coords,
        index=index,
        edges=edges,
        neighbors=neighbors,
        plaquette_rings=plaquette_rings,
        dist_matrix=dist_matrix,
    )


def double_edge_swap(
    edges: List[Edge],
    n_nodes: int,
    n_swaps: int,
    rng: np.random.Generator,
) -> List[Edge]:
    """
    Degree-preserving randomization using simple double-edge swaps.

    This is a lightweight null model. It preserves degrees approximately exactly
    when swaps are accepted, avoids self-loops and duplicate edges, and keeps
    the graph undirected.
    """
    edge_set = {tuple(sorted(e)) for e in edges}
    edge_list = list(edge_set)

    for _ in range(n_swaps):
        if len(edge_list) < 2:
            break
        idx1, idx2 = rng.choice(len(edge_list), size=2, replace=False)
        a, b = edge_list[idx1]
        c, d = edge_list[idx2]
        if len({a, b, c, d}) < 4:
            continue

        if rng.random() < 0.5:
            new1 = tuple(sorted((a, d)))
            new2 = tuple(sorted((c, b)))
        else:
            new1 = tuple(sorted((a, c)))
            new2 = tuple(sorted((b, d)))

        if new1[0] == new1[1] or new2[0] == new2[1]:
            continue
        if new1 in edge_set or new2 in edge_set:
            continue

        edge_set.remove((a, b))
        edge_set.remove((c, d))
        edge_set.add(new1)
        edge_set.add(new2)
        edge_list[idx1] = new1
        edge_list[idx2] = new2

    return sorted(edge_set)


def neighbors_from_edges(n: int, edges: Sequence[Edge]) -> List[List[int]]:
    neighbors: List[List[int]] = [[] for _ in range(n)]
    for i, j in edges:
        neighbors[i].append(j)
        neighbors[j].append(i)
    return neighbors


def graph_distances(n: int, neighbors: Sequence[Sequence[int]]) -> np.ndarray:
    dist = np.full((n, n), fill_value=32767, dtype=np.int16)
    for src in range(n):
        dist[src, src] = 0
        q: deque[int] = deque([src])
        while q:
            u = q.popleft()
            for v in neighbors[u]:
                if dist[src, v] > dist[src, u] + 1:
                    dist[src, v] = dist[src, u] + 1
                    q.append(v)
    return dist


def local_coherence_feedback(theta: np.ndarray, neighbors: Sequence[Sequence[int]]) -> np.ndarray:
    """
    Experimental local-coherence feedback term.

    This is NOT claimed as a final Wasserstein/Madelung operator. In R0, lambda
    should usually be 0. The term is included only to test whether a secondary
    feedback can add anything beyond primary XY phase coupling.
    """
    n = len(theta)
    c = np.zeros(n, dtype=float)
    for i in range(n):
        if neighbors[i]:
            c[i] = float(np.mean(np.cos(theta[np.array(neighbors[i])] - theta[i])))
    q = np.zeros(n, dtype=float)
    for i in range(n):
        if neighbors[i]:
            q[i] = float(np.mean(c[np.array(neighbors[i])] - c[i]))
    return q


def simulate_xy_phihex(
    *,
    neighbors: Sequence[Sequence[int]],
    jtheta: float,
    sigma: float,
    lam: float,
    steps: int,
    burn_in: int,
    dt: float,
    seed: int,
    omega_std: float,
) -> np.ndarray:
    """
    Euler-Maruyama simulation of the XY-Phi-Hex phase dynamics.

    Returns sampled theta states after burn-in.
    """
    if burn_in >= steps:
        raise ValueError("burn_in must be smaller than steps")

    rng = np.random.default_rng(seed)
    n = len(neighbors)
    theta = rng.uniform(0.0, 2.0 * math.pi, size=n)
    omega = rng.normal(0.0, omega_std, size=n)

    samples: List[np.ndarray] = []
    sqrt_dt = math.sqrt(dt)

    for t in range(steps):
        coupling = np.zeros(n, dtype=float)
        for i in range(n):
            if neighbors[i]:
                js = np.array(neighbors[i], dtype=int)
                coupling[i] = float(np.sum(np.sin(theta[js] - theta[i])))

        feedback = local_coherence_feedback(theta, neighbors) if lam != 0.0 else 0.0
        noise = rng.normal(0.0, 1.0, size=n)

        theta = theta + dt * (omega + jtheta * coupling + lam * feedback) + sigma * sqrt_dt * noise
        theta = np.asarray(wrap_angle(theta))

        if t >= burn_in:
            # Sample sparsely to reduce autocorrelation.
            if (t - burn_in) % 10 == 0:
                samples.append(theta.copy())

    return np.stack(samples, axis=0)


def g1_by_distance(theta_samples: np.ndarray, dist_matrix: np.ndarray) -> Dict[int, float]:
    """
    Compute g1(r)=mean cos(theta_i-theta_j) by graph distance.
    """
    n = theta_samples.shape[1]
    max_dist = int(np.max(dist_matrix[dist_matrix < 32767]))
    out: Dict[int, float] = {}

    for d in range(1, max_dist + 1):
        pairs = np.argwhere(dist_matrix == d)
        pairs = pairs[pairs[:, 0] < pairs[:, 1]]
        if len(pairs) == 0:
            continue
        vals = []
        for theta in theta_samples:
            vals.append(np.cos(theta[pairs[:, 0]] - theta[pairs[:, 1]]))
        out[d] = float(np.mean(np.concatenate(vals)))
    return out


def fit_decay(g1: Dict[int, float]) -> Dict[str, float | str]:
    """
    Compare algebraic and exponential decay fits using simple least squares on log(g1).

    Algebraic: log(g) = a - eta * log(r)
    Exponential: log(g) = a - r / xi

    Only positive g1 values are usable. This is an R0 diagnostic, not a final
    inferential model.
    """
    xs = []
    ys = []
    for r, g in sorted(g1.items()):
        if r >= 1 and g > 1e-8:
            xs.append(float(r))
            ys.append(math.log(float(g)))

    if len(xs) < 3:
        return {
            "fit_status": "INSUFFICIENT_POSITIVE_G1",
            "best_model": "none",
            "aic_power": float("nan"),
            "aic_exp": float("nan"),
            "eta_power": float("nan"),
            "xi_exp": float("nan"),
        }

    x = np.array(xs)
    y = np.array(ys)

    # Power law.
    Xp = np.column_stack([np.ones_like(x), -np.log(x)])
    beta_p, *_ = np.linalg.lstsq(Xp, y, rcond=None)
    resid_p = y - Xp @ beta_p
    rss_p = float(np.sum(resid_p ** 2))

    # Exponential.
    Xe = np.column_stack([np.ones_like(x), -x])
    beta_e, *_ = np.linalg.lstsq(Xe, y, rcond=None)
    resid_e = y - Xe @ beta_e
    rss_e = float(np.sum(resid_e ** 2))

    k = 2
    n = len(y)
    eps = 1e-12
    aic_p = n * math.log((rss_p + eps) / n) + 2 * k
    aic_e = n * math.log((rss_e + eps) / n) + 2 * k

    slope_e = float(beta_e[1])
    xi = 1.0 / slope_e if slope_e > 1e-12 else float("inf")

    return {
        "fit_status": "OK",
        "best_model": "power" if aic_p < aic_e else "exponential",
        "aic_power": float(aic_p),
        "aic_exp": float(aic_e),
        "eta_power": float(beta_p[1]),
        "xi_exp": float(xi),
    }


def plaquette_windings(theta: np.ndarray, plaquette_rings: Sequence[Sequence[int]]) -> np.ndarray:
    charges: List[int] = []
    for ring in plaquette_rings:
        if len(ring) < 3:
            continue
        total = 0.0
        for a, b in zip(ring, ring[1:] + ring[:1]):
            total += float(wrap_angle(theta[b] - theta[a]))
        charges.append(int(round(total / (2.0 * math.pi))))
    return np.array(charges, dtype=int)


def vortex_metrics(theta_samples: np.ndarray, plaquette_rings: Sequence[Sequence[int]]) -> Dict[str, float]:
    if not plaquette_rings:
        return {
            "vortex_density": float("nan"),
            "mean_abs_charge": float("nan"),
            "charge_balance_abs": float("nan"),
        }

    densities = []
    abs_charges = []
    balance = []
    for theta in theta_samples:
        q = plaquette_windings(theta, plaquette_rings)
        densities.append(float(np.mean(q != 0)))
        abs_charges.append(float(np.mean(np.abs(q))))
        balance.append(float(abs(np.sum(q))))

    return {
        "vortex_density": float(np.mean(densities)),
        "mean_abs_charge": float(np.mean(abs_charges)),
        "charge_balance_abs": float(np.mean(balance)),
    }


def helicity_analog(theta_samples: np.ndarray, edges: Sequence[Edge]) -> float:
    if not edges:
        return float("nan")
    vals = []
    ii = np.array([e[0] for e in edges], dtype=int)
    jj = np.array([e[1] for e in edges], dtype=int)
    for theta in theta_samples:
        vals.append(np.cos(theta[ii] - theta[jj]))
    return float(np.mean(np.concatenate(vals)))


def binder_like(theta_samples: np.ndarray) -> float:
    """
    Binder-like cumulant based on complex magnetization magnitude.
    U = 1 - <m^4> / (3 <m^2>^2)

    For XY models, alternative Binder definitions exist. This R0 indicator is
    used only as a finite-size diagnostic, not a standalone order proof.
    """
    m = np.abs(np.mean(np.exp(1j * theta_samples), axis=1))
    m2 = float(np.mean(m ** 2))
    m4 = float(np.mean(m ** 4))
    if m2 <= 1e-12:
        return float("nan")
    return float(1.0 - m4 / (3.0 * m2 * m2))


def bootstrap_ci(values: Sequence[float], seed: int, n_boot: int = 1000, alpha: float = 0.05) -> Tuple[float, float]:
    arr = np.array(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    means = []
    for _ in range(n_boot):
        sample = rng.choice(arr, size=len(arr), replace=True)
        means.append(float(np.mean(sample)))
    lo = float(np.quantile(means, alpha / 2.0))
    hi = float(np.quantile(means, 1.0 - alpha / 2.0))
    return lo, hi


def benjamini_hochberg(pvalues: Sequence[float], q: float = 0.05) -> List[bool]:
    """
    Benjamini-Hochberg FDR decision helper.
    Returns Boolean reject decisions in original order.
    """
    p = np.array(pvalues, dtype=float)
    m = len(p)
    order = np.argsort(p)
    sorted_p = p[order]

    k_max = -1
    for k, pv in enumerate(sorted_p, start=1):
        if pv <= (k / m) * q:
            k_max = k

    reject = np.zeros(m, dtype=bool)
    if k_max > 0:
        reject[order[:k_max]] = True
    return reject.tolist()


def run_condition(
    *,
    radius: int,
    jtheta: float,
    sigma: float,
    lam: float,
    seeds: int,
    steps: int,
    burn_in: int,
    dt: float,
    omega_std: float,
    null: str,
) -> Dict[str, object]:
    patch = build_hex_patch(radius)
    edges = patch.edges
    neighbors = patch.neighbors
    plaquettes = patch.plaquette_rings
    dist_matrix = patch.dist_matrix

    if null == "shuffled":
        rng = np.random.default_rng(12345 + radius)
        edges = double_edge_swap(edges, len(patch.coords), n_swaps=max(100, 10 * len(edges)), rng=rng)
        neighbors = neighbors_from_edges(len(patch.coords), edges)
        dist_matrix = graph_distances(len(patch.coords), neighbors)
        # Plaquettes are geometry-specific; shuffled graph has no meaningful hex plaquettes.
        plaquettes = []
    elif null == "hex":
        pass
    else:
        raise ValueError("null must be 'hex' or 'shuffled'")

    per_seed = []
    g1_accum: Dict[int, List[float]] = {}

    for seed in range(seeds):
        samples = simulate_xy_phihex(
            neighbors=neighbors,
            jtheta=jtheta,
            sigma=sigma,
            lam=lam,
            steps=steps,
            burn_in=burn_in,
            dt=dt,
            seed=seed,
            omega_std=omega_std,
        )
        g1 = g1_by_distance(samples, dist_matrix)
        vort = vortex_metrics(samples, plaquettes)
        hel = helicity_analog(samples, edges)
        bind = binder_like(samples)
        fit = fit_decay(g1)

        for d, val in g1.items():
            g1_accum.setdefault(d, []).append(val)

        per_seed.append({
            "seed": seed,
            "helicity_analog": hel,
            "binder_like": bind,
            "vortex_density": vort["vortex_density"],
            "mean_abs_charge": vort["mean_abs_charge"],
            "fit_best_model": fit["best_model"],
            "fit_status": fit["fit_status"],
        })

    g1_mean = {str(d): float(np.mean(vals)) for d, vals in sorted(g1_accum.items())}
    g1_ci = {
        str(d): {
            "lo": bootstrap_ci(vals, seed=999 + d)[0],
            "hi": bootstrap_ci(vals, seed=999 + d)[1],
        }
        for d, vals in sorted(g1_accum.items())
    }

    helicities = [x["helicity_analog"] for x in per_seed]
    binders = [x["binder_like"] for x in per_seed]
    vortex_densities = [x["vortex_density"] for x in per_seed if math.isfinite(float(x["vortex_density"]))]

    aggregate_fit = fit_decay({int(k): v for k, v in g1_mean.items()})

    return {
        "module": "PHY024_R0_XY_PhiHex_reference",
        "radius": radius,
        "n_nodes": len(patch.coords),
        "n_edges": len(edges),
        "n_plaquettes": len(plaquettes),
        "null_model": null,
        "parameters": {
            "jtheta": jtheta,
            "sigma": sigma,
            "lambda": lam,
            "seeds": seeds,
            "steps": steps,
            "burn_in": burn_in,
            "dt": dt,
            "omega_std": omega_std,
        },
        "metrics": {
            "g1_mean_by_distance": g1_mean,
            "g1_bootstrap_ci_by_distance": g1_ci,
            "helicity_analog_mean": float(np.mean(helicities)),
            "helicity_analog_ci": {
                "lo": bootstrap_ci(helicities, seed=555)[0],
                "hi": bootstrap_ci(helicities, seed=555)[1],
            },
            "binder_like_mean": float(np.mean(binders)),
            "binder_like_ci": {
                "lo": bootstrap_ci(binders, seed=556)[0],
                "hi": bootstrap_ci(binders, seed=556)[1],
            },
            "vortex_density_mean": float(np.mean(vortex_densities)) if vortex_densities else float("nan"),
            "vortex_density_ci": {
                "lo": bootstrap_ci(vortex_densities, seed=557)[0] if vortex_densities else float("nan"),
                "hi": bootstrap_ci(vortex_densities, seed=557)[1] if vortex_densities else float("nan"),
            },
            "aggregate_decay_fit": aggregate_fit,
        },
        "per_seed_summary": per_seed,
        "claim_guard": [
            "This script tests finite-size BKT-compatible indicators only.",
            "It does not prove a BKT transition.",
            "It does not make cosmological HQST claims.",
            "J_theta is primary; lambda/backreaction is secondary.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="PHY024-R0 XY-Phi-Hex reference run")
    parser.add_argument("--radius", type=int, default=3)
    parser.add_argument("--jtheta", type=float, default=0.4)
    parser.add_argument("--sigma", type=float, default=0.05)
    parser.add_argument("--lambda", dest="lam", type=float, default=0.0)
    parser.add_argument("--seeds", type=int, default=16)
    parser.add_argument("--steps", type=int, default=800)
    parser.add_argument("--burn-in", type=int, default=300)
    parser.add_argument("--dt", type=float, default=0.05)
    parser.add_argument("--omega-std", type=float, default=0.0)
    parser.add_argument("--null", choices=["hex", "shuffled"], default="hex")
    parser.add_argument("--out", type=str, default="")
    args = parser.parse_args()

    report = run_condition(
        radius=args.radius,
        jtheta=args.jtheta,
        sigma=args.sigma,
        lam=args.lam,
        seeds=args.seeds,
        steps=args.steps,
        burn_in=args.burn_in,
        dt=args.dt,
        omega_std=args.omega_std,
        null=args.null,
    )

    text = json.dumps(report, indent=2, sort_keys=True)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
