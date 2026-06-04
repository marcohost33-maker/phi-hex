#!/usr/bin/env python3
"""
PHY024-R2 Periodic Honeycomb XY-PhiHex Reference v0.1

This script improves R0/R1 by using periodic honeycomb topology and explicit
hex plaquette windings.

It does not prove BKT. It computes finite-size indicators and null-model checks.

Dependencies: Python >=3.10, NumPy.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import deque
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

import numpy as np

Edge = Tuple[int, int]


def wrap_angle(x):
    return (x + math.pi) % (2.0 * math.pi) - math.pi


@dataclass(frozen=True)
class HoneycombTorus:
    lx: int
    ly: int
    n: int
    edges: List[Edge]
    neighbors: List[List[int]]
    plaquettes: List[List[int]]
    x_boundary_edges: List[Edge]


def idx(x: int, y: int, sub: int, lx: int, ly: int) -> int:
    x %= lx
    y %= ly
    return 2 * (y * lx + x) + sub


def build_honeycomb_torus(lx: int, ly: int) -> HoneycombTorus:
    if lx < 3 or ly < 3:
        raise ValueError("lx and ly must be >= 3 for meaningful periodic honeycomb plaquettes")

    edges_set = set()
    x_boundary_edges = []

    # Sublattice convention:
    # A(x,y) connects to B(x,y), B(x-1,y), B(x,y-1).
    for x in range(lx):
        for y in range(ly):
            a = idx(x, y, 0, lx, ly)
            candidates = [
                (x, y, False),
                (x - 1, y, x == 0),
                (x, y - 1, False),
            ]
            for bx, by, crosses_x in candidates:
                b = idx(bx, by, 1, lx, ly)
                e = tuple(sorted((a, b)))
                edges_set.add(e)
                if crosses_x:
                    x_boundary_edges.append(e)

    edges = sorted(edges_set)
    n = 2 * lx * ly
    neighbors = [[] for _ in range(n)]
    for i, j in edges:
        neighbors[i].append(j)
        neighbors[j].append(i)

    # One ordered hexagonal plaquette per unit cell.
    plaquettes = []
    for x in range(lx):
        for y in range(ly):
            ring = [
                idx(x, y, 0, lx, ly),          # A(x,y)
                idx(x, y, 1, lx, ly),          # B(x,y)
                idx(x + 1, y, 0, lx, ly),      # A(x+1,y)
                idx(x + 1, y - 1, 1, lx, ly),  # B(x+1,y-1)
                idx(x + 1, y - 1, 0, lx, ly),  # A(x+1,y-1)
                idx(x, y - 1, 1, lx, ly),      # B(x,y-1)
            ]
            plaquettes.append(ring)

    return HoneycombTorus(
        lx=lx,
        ly=ly,
        n=n,
        edges=edges,
        neighbors=neighbors,
        plaquettes=plaquettes,
        x_boundary_edges=sorted(set(x_boundary_edges)),
    )


def graph_distances(n: int, neighbors: Sequence[Sequence[int]]) -> np.ndarray:
    dist = np.full((n, n), 32767, dtype=np.int16)
    for s in range(n):
        dist[s, s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            for v in neighbors[u]:
                if dist[s, v] > dist[s, u] + 1:
                    dist[s, v] = dist[s, u] + 1
                    q.append(v)
    return dist


def double_edge_swap(edges: List[Edge], n: int, swaps: int, rng: np.random.Generator) -> List[Edge]:
    edge_set = {tuple(sorted(e)) for e in edges}
    edge_list = list(edge_set)
    for _ in range(swaps):
        i1, i2 = rng.choice(len(edge_list), 2, replace=False)
        a, b = edge_list[i1]
        c, d = edge_list[i2]
        if len({a, b, c, d}) < 4:
            continue
        new1 = tuple(sorted((a, d)))
        new2 = tuple(sorted((c, b)))
        if new1[0] == new1[1] or new2[0] == new2[1]:
            continue
        if new1 in edge_set or new2 in edge_set:
            continue
        edge_set.remove((a, b))
        edge_set.remove((c, d))
        edge_set.add(new1)
        edge_set.add(new2)
        edge_list[i1] = new1
        edge_list[i2] = new2
    return sorted(edge_set)


def neighbors_from_edges(n: int, edges: Sequence[Edge]) -> List[List[int]]:
    out = [[] for _ in range(n)]
    for i, j in edges:
        out[i].append(j)
        out[j].append(i)
    return out


def simulate(
    neighbors: Sequence[Sequence[int]],
    jtheta: float,
    sigma: float,
    steps: int,
    burn_in: int,
    dt: float,
    seed: int,
    omega_std: float = 0.0,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n = len(neighbors)
    theta = rng.uniform(0.0, 2.0 * math.pi, n)
    omega = rng.normal(0.0, omega_std, n)
    samples = []
    sqrt_dt = math.sqrt(dt)

    for t in range(steps):
        coupling = np.zeros(n, dtype=float)
        for i in range(n):
            if neighbors[i]:
                js = np.array(neighbors[i], dtype=int)
                coupling[i] = np.sum(np.sin(theta[js] - theta[i]))
        theta = theta + dt * (omega + jtheta * coupling) + sigma * sqrt_dt * rng.normal(0.0, 1.0, n)
        theta = wrap_angle(theta)
        if t >= burn_in and (t - burn_in) % 10 == 0:
            samples.append(theta.copy())
    return np.stack(samples, axis=0)


def g1_by_distance(samples: np.ndarray, dist: np.ndarray) -> Dict[int, float]:
    out = {}
    maxd = int(np.max(dist[dist < 32767]))
    for d in range(1, maxd + 1):
        pairs = np.argwhere(dist == d)
        pairs = pairs[pairs[:, 0] < pairs[:, 1]]
        if len(pairs) == 0:
            continue
        vals = []
        for th in samples:
            vals.append(np.cos(th[pairs[:, 0]] - th[pairs[:, 1]]))
        out[d] = float(np.mean(np.concatenate(vals)))
    return out


def plaquette_winding(theta: np.ndarray, ring: Sequence[int]) -> int:
    total = 0.0
    for a, b in zip(ring, list(ring[1:]) + [ring[0]]):
        total += float(wrap_angle(theta[b] - theta[a]))
    return int(round(total / (2.0 * math.pi)))


def vortex_metrics(samples: np.ndarray, plaquettes: Sequence[Sequence[int]]) -> Dict[str, float]:
    if not plaquettes:
        return {"vortex_density": float("nan"), "pair_fraction": float("nan"), "mean_abs_charge": float("nan")}

    dens = []
    absq = []
    pair_fracs = []

    # Crude plaquette-center distance on plaquette index torus is omitted in R2 v0.1.
    # Pair fraction here is charge-balance proxy: min(Nplus,Nminus)/(Nplus+Nminus).
    for th in samples:
        q = np.array([plaquette_winding(th, p) for p in plaquettes], dtype=int)
        nonzero = q[q != 0]
        dens.append(float(len(nonzero) / len(q)))
        absq.append(float(np.mean(np.abs(q))))
        nplus = int(np.sum(q > 0))
        nminus = int(np.sum(q < 0))
        denom = nplus + nminus
        pair_fracs.append(float(min(nplus, nminus) / denom) if denom else 1.0)

    return {
        "vortex_density": float(np.mean(dens)),
        "pair_fraction": float(np.mean(pair_fracs)),
        "mean_abs_charge": float(np.mean(absq)),
    }


def energy(theta: np.ndarray, edges: Sequence[Edge], jtheta: float, twist_edges: Sequence[Edge] = (), phi: float = 0.0) -> float:
    twist_set = {tuple(sorted(e)) for e in twist_edges}
    e = 0.0
    for i, j in edges:
        d = theta[i] - theta[j]
        if tuple(sorted((i, j))) in twist_set:
            d += phi
        e += -jtheta * math.cos(float(d))
    return e


def twist_energy_curvature(samples: np.ndarray, edges: Sequence[Edge], twist_edges: Sequence[Edge], jtheta: float, delta: float = 1e-3) -> float:
    if not twist_edges or jtheta == 0:
        return float("nan")
    vals = []
    n = samples.shape[1]
    for th in samples:
        e0 = energy(th, edges, jtheta)
        ep = energy(th, edges, jtheta, twist_edges, +delta)
        em = energy(th, edges, jtheta, twist_edges, -delta)
        vals.append((ep + em - 2.0 * e0) / (delta * delta) / n)
    return float(np.mean(vals))


def edge_coherence(samples: np.ndarray, edges: Sequence[Edge]) -> float:
    vals = []
    ii = np.array([e[0] for e in edges], dtype=int)
    jj = np.array([e[1] for e in edges], dtype=int)
    for th in samples:
        vals.append(np.cos(th[ii] - th[jj]))
    return float(np.mean(np.concatenate(vals)))


def magnetization(samples: np.ndarray) -> float:
    return float(np.mean(np.abs(np.mean(np.exp(1j * samples), axis=1))))


def run_condition(lx: int, ly: int, jtheta: float, sigma: float, seeds: int, steps: int, burn_in: int, dt: float, null: str) -> Dict[str, object]:
    torus = build_honeycomb_torus(lx, ly)
    edges = torus.edges
    neighbors = torus.neighbors
    plaquettes = torus.plaquettes
    x_boundary_edges = torus.x_boundary_edges

    if null == "shuffled":
        rng = np.random.default_rng(1000 + lx * 17 + ly)
        edges = double_edge_swap(edges, torus.n, max(100, 20 * len(edges)), rng)
        neighbors = neighbors_from_edges(torus.n, edges)
        plaquettes = []
        x_boundary_edges = []
    elif null != "honeycomb":
        raise ValueError("null must be honeycomb or shuffled")

    dist = graph_distances(torus.n, neighbors)
    per_seed = []
    g1_acc = {}

    for seed in range(seeds):
        samples = simulate(neighbors, jtheta, sigma, steps, burn_in, dt, seed)
        g1 = g1_by_distance(samples, dist)
        for k, v in g1.items():
            g1_acc.setdefault(k, []).append(v)
        vort = vortex_metrics(samples, plaquettes)
        per_seed.append({
            "seed": seed,
            "edge_coherence": edge_coherence(samples, edges),
            "magnetization": magnetization(samples),
            "twist_energy_curvature_proxy": twist_energy_curvature(samples, edges, x_boundary_edges, jtheta),
            "vortex_density": vort["vortex_density"],
            "vortex_pair_fraction": vort["pair_fraction"],
            "mean_abs_charge": vort["mean_abs_charge"],
        })

    def mean_of(name: str) -> float:
        vals = [x[name] for x in per_seed if math.isfinite(float(x[name]))]
        return float(np.mean(vals)) if vals else float("nan")

    return {
        "module": "PHY024_R2_periodic_honeycomb_xy_phihex",
        "status": "R2_METHOD_HARDENING_NOT_FINAL_BKT_RESULT",
        "lx": lx,
        "ly": ly,
        "n_nodes": torus.n,
        "n_edges": len(edges),
        "n_plaquettes": len(plaquettes),
        "null_model": null,
        "parameters": {
            "jtheta": jtheta,
            "sigma": sigma,
            "seeds": seeds,
            "steps": steps,
            "burn_in": burn_in,
            "dt": dt,
        },
        "metrics": {
            "g1_mean_by_distance": {str(k): float(np.mean(v)) for k, v in sorted(g1_acc.items())},
            "edge_coherence_mean": mean_of("edge_coherence"),
            "magnetization_mean": mean_of("magnetization"),
            "twist_energy_curvature_proxy_mean": mean_of("twist_energy_curvature_proxy"),
            "vortex_density_mean": mean_of("vortex_density"),
            "vortex_pair_fraction_mean": mean_of("vortex_pair_fraction"),
            "mean_abs_charge_mean": mean_of("mean_abs_charge"),
        },
        "per_seed": per_seed,
        "claim_guard": [
            "R2 tests indicators only; no BKT proof.",
            "twist_energy_curvature_proxy is not full thermodynamic helicity modulus.",
            "Backreaction/lambda remains disabled.",
        ],
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--lx", type=int, default=4)
    p.add_argument("--ly", type=int, default=4)
    p.add_argument("--jtheta", type=float, default=0.4)
    p.add_argument("--sigma", type=float, default=0.05)
    p.add_argument("--seeds", type=int, default=8)
    p.add_argument("--steps", type=int, default=600)
    p.add_argument("--burn-in", type=int, default=250)
    p.add_argument("--dt", type=float, default=0.05)
    p.add_argument("--null", choices=["honeycomb", "shuffled"], default="honeycomb")
    p.add_argument("--out", type=str, default="")
    args = p.parse_args()

    report = run_condition(
        lx=args.lx,
        ly=args.ly,
        jtheta=args.jtheta,
        sigma=args.sigma,
        seeds=args.seeds,
        steps=args.steps,
        burn_in=args.burn_in,
        dt=args.dt,
        null=args.null,
    )
    text = json.dumps(report, indent=2, sort_keys=True)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
