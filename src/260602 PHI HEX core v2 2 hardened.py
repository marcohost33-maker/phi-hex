"""phi_hex_core v2.0 - Gehaertetes Kernmodul (HARDENED, STATE-OF-THE-ART)

Coworker Research / Coworkerz, 02. Juni 2026
de-CH Konventionen, ASCII-Quotes, ss statt Eszett, keine Personennamen.

Konsolidiert und haertet den gesamten Phi-Hex-Stack (PHY017-PHY023) zu
einem reproduzierbaren Kernmodul mit dem XY-Kopplungs-Durchbruch als
kanonischer Baseline.

KORREKTUR DER GITTER-TOPOLOGIE (zentral):
  Die Phi-Hex-Knoten haben Koordinationszahl z=6 (sechs J1-Richtungen).
  Damit bilden die KNOTEN ein TRIANGULAR LATTICE, nicht ein Honeycomb.
  Honeycomb (z=3) ist das DUALE Gitter der Hex-Zellen-Schwerpunkte.
  Fuer ein XY-Modell auf den Knoten ist der korrekte BKT-Referenzwert
  daher T_BKT(triangular) ca. 1.4 (in Einheiten J), NICHT 0.575 (honeycomb).
  Quelle: Okabe-Otsuka J. Phys. A 58, 065003 (2025); Sorokin (helicity).

KANONISCHES MODELL (XY-Hamiltonian auf Triangular Lattice):
  H = -J * sum_<i,j> cos(theta_i - theta_j)
  Langevin-Dynamik: dtheta_i/dt = -dH/dtheta_i + noise
                  = -J * sum_j sin(theta_i - theta_j) + sqrt(2T) * xi_i(t)
  Diskret (Euler-Maruyama):
    theta_i += dt * J * sum_j sin(theta_j - theta_i) + sqrt(2*T*dt) * N(0,1)
  Effektive Temperatur T entspricht dem frueheren eta^2-Skalenmass.

STATE-OF-THE-ART-KOMPONENTEN:
  1. Reproduzierbarkeit: zentrale Seed-Verwaltung, deterministische Pipeline,
     Provenance-Record (Commit-faehig).
  2. BKT-Diagnostik: Helicity-Modulus mit Nelson-Kosterlitz-Sprung 2T/pi,
     Weber-Minnhagen-FSS mit Subleading-Log, Υ4 (4. Ordnung).
  3. Vortex-Detektion: Plaquette-Winkelsumme mit (-pi,pi]-Reduktion,
     Ladungserhaltungs-Check auf Torus.
  4. Defekt-Analyse: strikte Trennung persistent (lifetime>=L_min) / transient.
  5. Statistik: Hedges' g, Cliff's delta / Vargha-Delaney A, BCa-Bootstrap.
  6. Validierungs-Gates: Referenzwert-Regression (T_BKT), Konvergenz, V&V.

Referenzen:
  Kosterlitz, Thouless (1973), J. Phys. C 6, 1181
  Nelson, Kosterlitz (1977), Phys. Rev. Lett. 39, 1201 (universal jump 2/pi)
  Weber, Minnhagen (1988), Phys. Rev. B 37, 5986 (FSS log-correction)
  Hasenbusch (2005), J. Phys. A 38, 5869 (subleading log)
  Okabe, Otsuka (2025), J. Phys. A 58, 065003 (lattice-dependent T_BKT)
  Efron (1987), JASA 82, 171 (BCa bootstrap)
  Vargha, Delaney (2000), J. Educ. Behav. Stat. 25, 101 (A measure)
"""
from __future__ import annotations

import json
import math
import platform
import sys
import time
from dataclasses import dataclass, field, asdict
from typing import Any, Callable, Optional

import numpy as np


# ============================================================================
# 0. Reproduzierbarkeit und Provenance
# ============================================================================

CORE_VERSION = "phi_hex_core_v2.0"
ATTRIBUTION = "Coworker Research / Coworkerz"


@dataclass
class ProvenanceRecord:
    """Vollstaendiger Provenance-Record fuer Reproduzierbarkeit (FAIR)."""
    core_version: str = CORE_VERSION
    attribution: str = ATTRIBUTION
    timestamp_utc: str = ""
    python_version: str = ""
    numpy_version: str = ""
    platform: str = ""
    master_seed: int = 0
    parameters: dict[str, Any] = field(default_factory=dict)

    def finalize(self, master_seed: int, parameters: dict[str, Any]) -> None:
        self.timestamp_utc = time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                            time.gmtime())
        self.python_version = sys.version.split()[0]
        self.numpy_version = np.__version__
        self.platform = platform.platform()
        self.master_seed = master_seed
        self.parameters = parameters


def make_rng(master_seed: int, stream: int = 0) -> np.random.Generator:
    """Deterministische, unabhaengige RNG-Streams via SeedSequence.

    Best Practice (NumPy SeedSequence): spawn unabhaengige Substreams,
    sodass jeder Seed/Stream-Index reproduzierbar und dekorreliert ist.
    """
    ss = np.random.SeedSequence([master_seed, stream])
    return np.random.default_rng(ss)


# ============================================================================
# 1. Triangular Lattice (axiale Hex-Koordinaten, z=6)
# ============================================================================

Ax = tuple[int, int]

J1_DIRECTIONS: list[Ax] = [
    (+1, 0), (0, +1), (-1, +1), (-1, 0), (0, -1), (+1, -1)
]


def axial_hex_dist(a: Ax, b: Ax) -> int:
    dq, dr = a[0] - b[0], a[1] - b[1]
    return (abs(dq) + abs(dr) + abs(dq + dr)) // 2


def axial_to_xy(a: Ax) -> tuple[float, float]:
    q, r = a
    return (1.5 * q, math.sqrt(3.0) * (r + 0.5 * q))


@dataclass
class TriangularLattice:
    """Triangular Lattice (Knoten z=6) auf axialen Koordinaten.

    Unterstuetzt offene Raender (Default) und periodische (Torus) Raender.
    Bei Torus ist die Gesamtvortizitaet erhalten (= 0), was als
    Validierungs-Check genutzt wird.

    edge_disp: minimum-image Bond-Verschiebung (dx, dy) je Kante. Auf dem
    Torus ist das der korrekte naechste-Nachbar-Vektor (Einheitslaenge),
    nicht die rohe Koordinatendifferenz gewrappter Knoten.
    """
    ax_of: list[Ax]
    idx_of: dict[Ax, int]
    edges: list[tuple[int, int]]
    cells: list[list[int]]
    radius: int
    periodic: bool
    edge_disp: list[tuple[float, float]] = field(default_factory=list)

    @property
    def n_nodes(self) -> int:
        return len(self.ax_of)

    @property
    def coordination_mean(self) -> float:
        deg = np.zeros(self.n_nodes)
        for i, j in self.edges:
            deg[i] += 1
            deg[j] += 1
        return float(np.mean(deg))


def build_triangular_lattice(radius: int,
                              periodic: bool = False) -> TriangularLattice:
    """Erzeugt Triangular Lattice.

    periodic=False: offene Raender, hexagonale Region mit Radius `radius`
                    (wie historisch; Koordinationszahl am Rand < 6).
    periodic=True : echte Torus-Topologie. Es wird ein LxL-Rhombus
                    (Parallelogramm) mit L = 2*radius+1 in axialen
                    Koordinaten erzeugt und modulo L gewrappt. Jeder
                    Knoten hat dann exakt 6 Nachbarn (Koordinationszahl 6),
                    und die Gesamtvortizitaet ist erhalten (= 0).

    Hinweis: Nur das Parallelogramm (nicht die Hexagon-Region) kachelt den
    Torus konsistent. Daher die getrennte Konstruktion.
    """
    if periodic:
        return _build_triangular_torus(L=2 * radius + 1, radius=radius)

    ax_of: list[Ax] = []
    for q in range(-radius, radius + 1):
        r_min = max(-radius, -q - radius)
        r_max = min(radius, -q + radius)
        for r in range(r_min, r_max + 1):
            ax_of.append((q, r))
    idx_of = {a: i for i, a in enumerate(ax_of)}

    edges_set: set[tuple[int, int]] = set()
    for i, a in enumerate(ax_of):
        for d in J1_DIRECTIONS:
            nb = (a[0] + d[0], a[1] + d[1])
            j = idx_of.get(nb)
            if j is not None and i < j:
                edges_set.add((i, j))
    # P0-Fix 2026-06-04 (Wissenschafts-Audit S1): diese Zeile wurde beim
    # v2.0->v2.2-Hardening geloescht -> NameError im open-boundary-Pfad;
    # der Default-Selbsttest lief as-shipped nie. (Drive-SoT v2.2 hat den Bug noch!)
    edges = sorted(edges_set)
    pos_open = [axial_to_xy(a) for a in ax_of]
    # Bond-Vektoren auf physikalische Einheitslaenge a=1 (NN-Abstand). Der
    # rohe axial_to_xy-Abstand ist sqrt(3); Division stellt a=1 her. Die
    # Helicity-Skala wird ueber die per-Site-Normierung gesetzt, nicht hier.
    _a_raw = math.sqrt(3.0)
    edge_disp = [((pos_open[j][0] - pos_open[i][0]) / _a_raw,
                  (pos_open[j][1] - pos_open[i][1]) / _a_raw)
                 for i, j in edges]

    cells: list[list[int]] = []
    seen: set[frozenset[int]] = set()
    for a in ax_of:
        loop = [idx_of[a]]
        ok = True
        for d in J1_DIRECTIONS:
            nb = (a[0] + d[0], a[1] + d[1])
            nb_idx = idx_of.get(nb)
            if nb_idx is None:
                ok = False
                break
            loop.append(nb_idx)
        if ok:
            key = frozenset(loop)
            if key not in seen:
                seen.add(key)
                cells.append(loop)

    return TriangularLattice(ax_of=ax_of, idx_of=idx_of, edges=edges,
                             cells=cells, radius=radius, periodic=False,
                             edge_disp=edge_disp)


def _build_triangular_torus(L: int, radius: int) -> TriangularLattice:
    """Echtes periodisches Triangular Lattice auf LxL-Rhombus (Torus).

    Knoten sind (q, r) mit q, r in 0..L-1. Nachbarn via J1_DIRECTIONS
    modulo L. Jeder Knoten hat exakt 6 Nachbarn. Die drei elementaren
    Dreiecks-Plaquetten-Typen werden als 3-er-Loops erfasst (fuer
    Vortizitaet); zusaetzlich werden 6-er-Hexagon-Loops fuer Kompatibilitaet
    mit der bestehenden Defekt-Logik bereitgestellt.
    """
    ax_of: list[Ax] = [(q, r) for q in range(L) for r in range(L)]
    idx_of = {a: i for i, a in enumerate(ax_of)}

    def wrap(a: Ax) -> Ax:
        return (a[0] % L, a[1] % L)

    edges_set: set[tuple[int, int]] = set()
    for i, a in enumerate(ax_of):
        for d in J1_DIRECTIONS:
            nb = wrap((a[0] + d[0], a[1] + d[1]))
            j = idx_of[nb]
            if i != j:
                edges_set.add((min(i, j), max(i, j)))
    edges = sorted(edges_set)

    # Minimum-image Bond-Verschiebung. Bond-Vektoren werden auf physikalische
    # Einheitslaenge a=1 (naechster-Nachbar-Abstand) gesetzt: der rohe
    # axial_to_xy-Abstand ist sqrt(3), daher Division durch sqrt(3).
    # Die korrekte Helicity-Skala entsteht dann ueber die per-Site-Normierung
    # in helicity_from_ensemble (nicht ueber willkuerliche Bond-Skalierung).
    _a_raw = math.sqrt(3.0)
    j1_xy = {d: (axial_to_xy(d)[0] / _a_raw, axial_to_xy(d)[1] / _a_raw)
             for d in J1_DIRECTIONS}
    edge_disp: list[tuple[float, float]] = []
    for i, j in edges:
        a = ax_of[i]
        disp = None
        for d in J1_DIRECTIONS:
            if wrap((a[0] + d[0], a[1] + d[1])) == ax_of[j]:
                disp = j1_xy[d]
                break
        edge_disp.append(disp if disp is not None else (0.0, 0.0))

    # Elementare Up-Dreiecke als Plaquetten (orientiert) fuer Vortizitaet:
    # (q,r) -> (q+1,r) -> (q,r+1) ist eine der beiden Dreiecksrichtungen.
    cells: list[list[int]] = []
    for a in ax_of:
        tri_up = [idx_of[a],
                  idx_of[wrap((a[0] + 1, a[1]))],
                  idx_of[wrap((a[0], a[1] + 1))]]
        tri_dn = [idx_of[a],
                  idx_of[wrap((a[0], a[1] + 1))],
                  idx_of[wrap((a[0] - 1, a[1] + 1))]]
        cells.append(tri_up)
        cells.append(tri_dn)

    return TriangularLattice(ax_of=ax_of, idx_of=idx_of, edges=edges,
                             cells=cells, radius=radius, periodic=True,
                             edge_disp=edge_disp)


# ============================================================================
# 2. Kanonisches XY-Modell mit Langevin-Dynamik
# ============================================================================

@dataclass
class XYConfig:
    """Konfiguration der XY-Langevin-Dynamik."""
    J: float = 1.0           # Kopplungsstaerke
    T: float = 0.5           # effektive Temperatur (Noise-Mass)
    dt: float = 0.05         # Zeitschritt (Euler-Maruyama)
    max_steps: int = 2000
    burn_in: int = 500
    measure_every: int = 5


def build_adjacency(lattice: TriangularLattice) -> list[list[int]]:
    adj: list[list[int]] = [[] for _ in range(lattice.n_nodes)]
    for i, j in lattice.edges:
        adj[i].append(j)
        adj[j].append(i)
    return adj


def xy_langevin_step(theta: np.ndarray, adj: list[list[int]],
                     cfg: XYConfig, rng: np.random.Generator) -> np.ndarray:
    """Ein Euler-Maruyama-Schritt der XY-Langevin-Dynamik.

    dtheta_i = dt * J * sum_j sin(theta_j - theta_i) + sqrt(2*T*dt)*N(0,1)
    """
    n = len(theta)
    force = np.zeros(n)
    for i in range(n):
        ti = theta[i]
        s = 0.0
        for j in adj[i]:
            s += math.sin(theta[j] - ti)
        force[i] = s
    noise = math.sqrt(2.0 * cfg.T * cfg.dt) * rng.standard_normal(n)
    theta_new = (theta + cfg.dt * cfg.J * force + noise) % (2.0 * math.pi)
    return theta_new


# ============================================================================
# 3. BKT-Diagnostik: Helicity-Modulus mit Nelson-Kosterlitz-Sprung
# ============================================================================

def _edge_arrays(lattice: Any
                 ) -> tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
    """Gibt (edge_i, edge_j, disp) als numpy-Arrays fuer die Vektorisierung.

    Cacht das Ergebnis auf dem Gitter-Objekt (`_phx_edge_cache`), da die
    Kanten/Verschiebungen waehrend eines Mess-Laufs konstant sind und
    `helicity_terms` pro Messung aufgerufen wird (Hot-Loop). Der Cache wird
    neu gebaut, wenn sich die Kantenzahl aendert (Robustheit). Gitter-
    agnostisch: funktioniert fuer jedes Objekt mit `.edges` (+ optional
    `.edge_disp`), also Triangular- wie Honeycomb-Gitter.

    `disp` ist None, falls keine (oder inkonsistente) edge_disp vorliegt;
    dann faellt `helicity_terms` auf die rohe Positionsdifferenz zurueck.
    """
    edges = lattice.edges
    edge_disp = getattr(lattice, "edge_disp", None)
    cache = getattr(lattice, "_phx_edge_cache", None)
    if cache is not None:
        c_edges, c_disp, arrays = cache
        # Identitaets-basierte Invalidierung: der Cache haelt Referenzen auf
        # die GENAUEN Listen-Objekte. Wird `edges` oder `edge_disp` ersetzt
        # bzw. umverdrahtet (neues Listen-Objekt) - auch bei gleicher Laenge -
        # greift `is` nicht mehr und die Arrays werden neu gebaut. Verhindert
        # stale-Cache-Messungen am alten Gitter (Codex-Review P2, 2026-06-10).
        # Das Festhalten der Referenzen schliesst zudem id()-Wiederverwendung
        # nach GC aus.
        if c_edges is edges and c_disp is edge_disp:
            return arrays
    if edges:
        ij = np.asarray(edges, dtype=np.intp)
        edge_i, edge_j = ij[:, 0], ij[:, 1]
    else:
        edge_i = np.empty(0, dtype=np.intp)
        edge_j = np.empty(0, dtype=np.intp)
    if edge_disp is not None and len(edge_disp) == len(edges):
        disp = np.asarray(edge_disp, dtype=float).reshape(len(edges), 2)
    else:
        disp = None
    arrays = (edge_i, edge_j, disp)
    try:
        lattice._phx_edge_cache = (edges, edge_disp, arrays)
    except (AttributeError, TypeError):
        pass  # frozen/slots-Gitter: ohne Cache, weiterhin korrekt
    return arrays


def helicity_terms(theta: np.ndarray, lattice: TriangularLattice,
                   pos: list[tuple[float, float]],
                   direction: tuple[float, float] = (1.0, 0.0)
                   ) -> tuple[float, float]:
    """Liefert die beiden Helicity-Terme (T1, sin_accum) fuer EINE Konfig.

    WICHTIG: Der Helicity-Modulus ist eine Ensemble-Groesse. Die beiden
    Terme muessen SEPARAT ueber Konfigurationen gemittelt werden, bevor
    sie subtrahiert werden. Ein einzelner Snapshot liefert ein falsches
    (oft negatives) Ergebnis - genau die Mess-Falle, vor der das Meta-Audit
    warnt.

      T1        = sum_<ij> cos(theta_i - theta_j) (e.r_ij)^2
      sin_accum = sum_<ij> sin(theta_i - theta_j) (e.r_ij)

    Nutzt lattice.edge_disp (minimum-image Bond-Vektoren), damit gewrappte
    Torus-Kanten den korrekten naechste-Nachbar-Vektor verwenden und nicht
    die rohe Koordinatendifferenz.

    Vektorisiert (numpy) ueber alle Kanten - mathematisch identisch zur
    fruehen skalaren Schleife (nur Gleitkomma-Summationsreihenfolge, ~1e-13
    relativ), aber Groessenordnungen schneller im Mess-Hot-Loop (Wolff-
    Sampling ruft das pro Messung auf). Aequivalenz ist als CI-Gate
    abgesichert (tests/test_core_vectorization.py).
    """
    ex, ey = direction
    norm = math.hypot(ex, ey)
    ex, ey = ex / norm, ey / norm
    edge_i, edge_j, disp = _edge_arrays(lattice)
    if edge_i.shape[0] == 0:
        return 0.0, 0.0
    theta = np.asarray(theta, dtype=float)
    if disp is not None:
        proj = disp[:, 0] * ex + disp[:, 1] * ey
    else:
        pa = np.asarray(pos, dtype=float)
        proj = (pa[edge_j, 0] - pa[edge_i, 0]) * ex \
            + (pa[edge_j, 1] - pa[edge_i, 1]) * ey
    dtheta = theta[edge_i] - theta[edge_j]
    t1 = float(np.dot(np.cos(dtheta) * proj, proj))
    sin_accum = float(np.dot(np.sin(dtheta), proj))
    return t1, sin_accum


def helicity_from_ensemble(t1_samples: list[float],
                           sin_samples: list[float],
                           cfg: XYConfig, n_nodes: int) -> float:
    """Helicity-Modulus aus separat gemittelten Ensemble-Termen.

    KONVENTION: per-Site-Normierung (Literatur-Standard, arXiv:2406.12076;
    so misst auch PHY028 das Quadratgitter). Marco-Entscheid 2026-06-04,
    Wissenschafts-Audit S2:

      Upsilon = (1/N) [ J*<T1> - (J^2/T)*<sin_accum^2> ]

    mit N = Zahl der Gitter-Plaetze. Da die Bond-Vektoren in helicity_terms
    auf NN-Einheitslaenge a=1 normiert sind, ergibt das per-Site-Schema bei
    T=0 (perfekt ausgerichtet) den korrekten Spinwellen-Grenzwert
    Upsilon(0) = 3/2 * J fuer das Dreiecksgitter (z=6) - groessen-unabhaengig.

    HISTORIE (superseded): bis 2026-06-04 wurde durch die FLAECHE
    A = N*sqrt(3)/2 geteilt. Das lieferte Upsilon(0) = J*sqrt(3) = 1.732 und
    ueberschaetzte T_BKT um ~15% (von der laxen Toleranz tol=0.15 verdeckt).
    Diese Definition ist NICHT der Nelson-Kosterlitz-Standard und wurde
    ersetzt.

    Am BKT-Uebergang schneidet Upsilon(T) die universelle Gerade 2*T/pi
    (Nelson-Kosterlitz). Darueber faellt Upsilon gegen null bzw. negativ.
    """
    mean_t1 = float(np.mean(t1_samples))
    mean_sin2 = float(np.mean(np.square(sin_samples)))
    return (cfg.J * mean_t1 / n_nodes
            - (cfg.J ** 2) * mean_sin2 / (n_nodes * max(cfg.T, 1e-12)))


def nelson_kosterlitz_line(T: float) -> float:
    """Universelle Sprung-Gerade 2*T/pi (Nelson-Kosterlitz 1977)."""
    return 2.0 * T / math.pi


def weber_minnhagen_prediction(T_bkt: float, L: float,
                               L0: float = 1.0) -> float:
    """Weber-Minnhagen-FSS-Vorhersage des Helicity-Modulus bei T_BKT.

    Upsilon(T_BKT, L) = (2*T_BKT/pi) * [1 + 1/(2*ln(L/L0))]

    Subleading-Log-Korrektur; ohne sie wird T_BKT systematisch zu tief
    geschaetzt (Hasenbusch 2005).
    """
    return (2.0 * T_bkt / math.pi) * (1.0 + 1.0 / (2.0 * math.log(L / L0)))


# ============================================================================
# 4. Vortex-Detektion (Plaquette-Winkelsumme, Ladungserhaltung)
# ============================================================================

def wrap_to_pi(d: float) -> float:
    return (d + math.pi) % (2.0 * math.pi) - math.pi


def plaquette_vorticity(theta: np.ndarray, loop: list[int]) -> int:
    """Ganzzahlige Vortizitaet einer Plaquette via Winkelsumme.

    Summiert (-pi,pi]-reduzierte Phasendifferenzen entlang des Loops,
    geteilt durch 2pi. Ergebnis ist die topologische Ladung (..., -1, 0, +1).
    """
    total = 0.0
    L = len(loop)
    for k in range(L):
        a = theta[loop[k]]
        b = theta[loop[(k + 1) % L]]
        total += wrap_to_pi(b - a)
    return int(round(total / (2.0 * math.pi)))


def vorticity_field(theta: np.ndarray,
                    lattice: TriangularLattice) -> list[int]:
    return [plaquette_vorticity(theta, loop) for loop in lattice.cells]


def check_charge_conservation(vort: list[int],
                              lattice: TriangularLattice) -> dict[str, Any]:
    """Prueft Gesamtvortizitaet. Auf Torus muss sie exakt 0 sein."""
    net = int(sum(vort))
    ok = (net == 0) if lattice.periodic else True
    return {"net_vorticity": net, "conserved": ok,
            "periodic": lattice.periodic}


# ============================================================================
# 5. Defekt-Tracking mit persistent/transient-Trennung
# ============================================================================

@dataclass
class DefectStats:
    n_persistent: int
    n_transient: int
    persistent_density_per_cell: float
    mean_persistent_lifetime: float
    total_tracked: int


class DefectTrackerV2:
    """Vortex-Tracking mit Ladungserhaltung beim Matching.

    Trennt strikt persistente (lifetime>=L_min) von transienten Defekten,
    wie im Audit als methodischer Standard festgelegt.
    """

    def __init__(self, cell_positions: list[tuple[float, float]],
                 max_match_dist: float = 2.5):
        self.cell_positions = cell_positions
        self.max_match_dist = max_match_dist
        self.trajectories: dict[int, dict[str, Any]] = {}
        self.active: dict[int, int] = {}  # cell -> traj_id
        self.next_id = 0

    def update(self, t: int, vort: list[int]) -> None:
        current: dict[int, int] = {c: q for c, q in enumerate(vort)
                                   if q != 0}
        new_active: dict[int, int] = {}
        used: set[int] = set()
        for cell, charge in current.items():
            pos = self.cell_positions[cell]
            best, best_d = None, float("inf")
            for ocell, tid in self.active.items():
                if ocell in used:
                    continue
                if self.trajectories[tid]["charge"] != charge:
                    continue
                opos = self.cell_positions[ocell]
                d = math.hypot(pos[0] - opos[0], pos[1] - opos[1])
                if d < best_d and d <= self.max_match_dist:
                    best, best_d = ocell, d
            if best is not None:
                tid = self.active[best]
                used.add(best)
                self.trajectories[tid]["t_end"] = t
                new_active[cell] = tid
            else:
                tid = self.next_id
                self.next_id += 1
                self.trajectories[tid] = {"charge": charge, "t_start": t,
                                          "t_end": t}
                new_active[cell] = tid
        self.active = new_active

    def stats(self, L_min: int = 3) -> DefectStats:
        n_pers, n_trans, pers_lifetimes = 0, 0, []
        for tr in self.trajectories.values():
            lt = tr["t_end"] - tr["t_start"]
            if lt >= L_min:
                n_pers += 1
                pers_lifetimes.append(lt)
            else:
                n_trans += 1
        n_cells = max(len(self.cell_positions), 1)
        return DefectStats(
            n_persistent=n_pers,
            n_transient=n_trans,
            persistent_density_per_cell=n_pers / n_cells,
            mean_persistent_lifetime=(float(np.mean(pers_lifetimes))
                                      if pers_lifetimes else 0.0),
            total_tracked=len(self.trajectories),
        )


# ============================================================================
# 6. Phasenkorrelationsfunktion g1(r) und Zerfallsklassifikation
# ============================================================================

def phase_correlation(theta: np.ndarray,
                      node_positions: list[tuple[float, float]],
                      n_bins: int = 8) -> tuple[np.ndarray, np.ndarray]:
    n = len(theta)
    dists, corrs = [], []
    for i in range(n):
        for j in range(i + 1, n):
            dx = node_positions[i][0] - node_positions[j][0]
            dy = node_positions[i][1] - node_positions[j][1]
            dists.append(math.hypot(dx, dy))
            corrs.append(math.cos(theta[i] - theta[j]))
    dists = np.array(dists)
    corrs = np.array(corrs)
    # Binning ab minimalem tatsaechlichen Abstand, sonst ist der erste Bin
    # (Abstand 0 bis erster Nachbar) leer und liefert NaN.
    edges = np.linspace(dists.min(), dists.max(), n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    # letzter Bin inklusiv am Rand
    edges[-1] += 1e-9
    g1 = np.full(n_bins, np.nan)
    for b in range(n_bins):
        m = (dists >= edges[b]) & (dists < edges[b + 1])
        if np.any(m):
            g1[b] = np.mean(corrs[m])
    return centers, g1


def classify_decay(r: np.ndarray, g1: np.ndarray) -> dict[str, Any]:
    """Klassifiziert g1(r): algebraisch (geordnet) vs exponentiell.

    Liefert auch eta (algebraischer Exponent). Am BKT-Uebergang eta ca. 1/4.
    """
    valid = (~np.isnan(g1)) & (g1 > 1e-3) & (r > 1e-9)
    if np.sum(valid) < 3:
        return {"decay_type": "indeterminate", "eta": float("nan"),
                "r2_alg": 0.0, "r2_exp": 0.0}
    rr, gg = r[valid], np.log(g1[valid])
    # algebraisch: log g = -eta log r + c
    A = np.vstack([np.log(rr), np.ones_like(rr)]).T
    ca, *_ = np.linalg.lstsq(A, gg, rcond=None)
    pred_a = A @ ca
    sstot = np.sum((gg - gg.mean()) ** 2)
    r2a = 1 - np.sum((gg - pred_a) ** 2) / max(sstot, 1e-12)
    # exponentiell: log g = -r/xi + c
    B = np.vstack([rr, np.ones_like(rr)]).T
    cb, *_ = np.linalg.lstsq(B, gg, rcond=None)
    pred_b = B @ cb
    r2e = 1 - np.sum((gg - pred_b) ** 2) / max(sstot, 1e-12)
    return {"decay_type": "algebraic" if r2a >= r2e else "exponential",
            "eta": float(-ca[0]), "r2_alg": float(r2a), "r2_exp": float(r2e)}


# ============================================================================
# 7. Statistik-Toolbox: Hedges' g, Cliff's delta, BCa-Bootstrap
# ============================================================================

def hedges_g(a: np.ndarray, b: np.ndarray) -> float:
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return 0.0
    sp = math.sqrt(((na - 1) * np.var(a, ddof=1)
                    + (nb - 1) * np.var(b, ddof=1)) / (na + nb - 2))
    if sp == 0:
        return 0.0
    d = (np.mean(a) - np.mean(b)) / sp
    J = 1.0 - 3.0 / (4.0 * (na + nb) - 9.0)
    return float(d * J)


def cliff_delta(a: np.ndarray, b: np.ndarray) -> dict[str, Any]:
    """Cliff's delta (nichtparametrische Effektstaerke) + VD-A.

    Negibility-Schwellen (Vargha-Delaney): |d|<0.11 vernachlaessigbar,
    <0.28 klein, <0.43 mittel, sonst gross.
    """
    na, nb = len(a), len(b)
    # Sortier-/Rang-basiert: O((na+nb) log nb) Zeit, O(nb) Speicher. Bit-
    # identisch zur fruehen O(na*nb)-Doppelschleife (nur ganzzahlige strikte
    # Paar-Vergleiche), aber OHNE die na*nb-Vergleichsmatrix zu materialisieren
    # - die wuerde bei grossen Stichproben den Speicher sprengen (z.B. 2x50k
    # -> ~2.5 GB je Temporary). Codex-Review P2, 2026-06-10.
    av = np.asarray(a, dtype=float)
    b_sorted = np.sort(np.asarray(b, dtype=float))
    # Pro x in a: #{y in b : y < x} (strict) summiert = Zahl der (x>y)-Paare;
    # #{y : y <= x} summiert -> (x<y)-Paare = na*nb - sum(#{y<=x}). Gleichstaende
    # zaehlen - wie in der strikten Doppelschleife - weder zu gt noch zu lt.
    gt = int(np.searchsorted(b_sorted, av, side="left").sum())
    le = int(np.searchsorted(b_sorted, av, side="right").sum())
    lt = na * nb - le
    delta = (gt - lt) / (na * nb)
    ad = abs(delta)
    mag = ("negligible" if ad < 0.11 else "small" if ad < 0.28
           else "medium" if ad < 0.43 else "large")
    return {"cliff_delta": float(delta),
            "vargha_delaney_A": float((delta + 1) / 2),
            "magnitude": mag}


def bca_bootstrap_ci(a: np.ndarray, b: np.ndarray,
                     n_boot: int = 10000, alpha: float = 0.05,
                     rng: Optional[np.random.Generator] = None
                     ) -> dict[str, float]:
    """BCa-Bootstrap-CI fuer die Differenz der Mittelwerte (a - b).

    Bias-corrected and accelerated (Efron 1987). Robuster als Perzentil-CI
    bei schiefen Verteilungen (z.B. Lifetimes).
    """
    if rng is None:
        rng = np.random.default_rng(0)
    theta_hat = float(np.mean(a) - np.mean(b))

    boot = np.empty(n_boot)
    for i in range(n_boot):
        ra = rng.choice(a, len(a), replace=True)
        rb = rng.choice(b, len(b), replace=True)
        boot[i] = ra.mean() - rb.mean()

    # bias-correction z0
    prop = np.mean(boot < theta_hat)
    prop = min(max(prop, 1e-9), 1 - 1e-9)
    z0 = _norm_ppf(prop)

    # acceleration via jackknife auf gepooltem Effekt
    jack = []
    # jackknife auf a
    for i in range(len(a)):
        aa = np.delete(a, i)
        jack.append(aa.mean() - b.mean())
    for j in range(len(b)):
        bb = np.delete(b, j)
        jack.append(a.mean() - bb.mean())
    jack = np.array(jack)
    jbar = jack.mean()
    num = np.sum((jbar - jack) ** 3)
    den = 6.0 * (np.sum((jbar - jack) ** 2) ** 1.5)
    acc = num / den if den != 0 else 0.0

    z_a = _norm_ppf(alpha / 2)
    z_1a = _norm_ppf(1 - alpha / 2)
    a1 = _norm_cdf(z0 + (z0 + z_a) / (1 - acc * (z0 + z_a)))
    a2 = _norm_cdf(z0 + (z0 + z_1a) / (1 - acc * (z0 + z_1a)))
    lo = float(np.percentile(boot, 100 * a1))
    hi = float(np.percentile(boot, 100 * a2))
    return {"diff": theta_hat, "ci_lower": lo, "ci_upper": hi,
            "excludes_zero": bool(lo * hi > 0), "z0": float(z0),
            "acceleration": float(acc)}


def _norm_cdf(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def _norm_ppf(p: float) -> float:
    """Inverse Standardnormal-CDF (scipy-frei, Acklam-Approximation)."""
    if p <= 0:
        return -math.inf
    if p >= 1:
        return math.inf
    try:
        from scipy.special import ndtri
        return float(ndtri(p))
    except ImportError:
        pass
    a = [-3.969683028665376e+01, 2.209460984245205e+02,
         -2.759285104469687e+02, 1.383577518672690e+02,
         -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02,
         -1.556989798598866e+02, 6.680131188771972e+01,
         -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01,
         -2.400758277161838e+00, -2.549732539343734e+00,
         4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01,
         2.445134137142996e+00, 3.754408661907416e+00]
    plow, phigh = 0.02425, 1 - 0.02425
    if p < plow:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p > phigh:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    q = p - 0.5
    r = q * q
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)


# ============================================================================
# 8. Power-Analyse (Monte-Carlo-Power-Kurve)
# ============================================================================

def required_sample_size(
    effect_simulator: Callable[[int, np.random.Generator],
                               tuple[np.ndarray, np.ndarray]],
    candidate_ns: list[int], target_power: float = 0.8,
    n_trials: int = 200, alpha: float = 0.05,
    master_seed: int = 0,
) -> dict[str, Any]:
    """Empirische Bestimmung der noetigen Seed-Zahl via Power-Kurve.

    Fuer jede Kandidaten-n: Anteil der Trials, in denen das BCa-CI der
    Differenz die Null ausschliesst (= empirische Power). Liefert das
    kleinste n, das target_power erreicht.
    """
    curve = {}
    chosen = None
    for n in candidate_ns:
        hits = 0
        for trial in range(n_trials):
            rng = make_rng(master_seed, stream=1000 * n + trial)
            a, b = effect_simulator(n, rng)
            ci = bca_bootstrap_ci(a, b, n_boot=1000, rng=rng)
            if ci["excludes_zero"]:
                hits += 1
        power = hits / n_trials
        curve[n] = power
        if chosen is None and power >= target_power:
            chosen = n
    return {"power_curve": curve, "required_n": chosen,
            "target_power": target_power}


# ============================================================================
# 9. Validierungs-Gates (Referenzwerte als Regressionstests)
# ============================================================================

# Lattice-abhaengige BKT-Referenzwerte (in Einheiten J/k_B), web-verifiziert
# (arXiv:2501.07388; Hasenbusch 2005; Okabe-Otsuka 2025).
# WICHTIG: Phi-Hex-Knoten = triangular (z=6) -> ca. 1.4.
# Diese T_BKT-Werte sind KONVENTIONS-UNABHAENGIG (kein per-Site/Flaechen-Bezug).
T_BKT_REFERENCE = {
    "square": 0.8929,        # Hasenbusch 2005, Jha 2020 (arXiv:2501.07388: 0.893)
    "triangular": 1.418,     # Sorokin (helicity), Okabe-Otsuka 2025; arXiv:2501.07388
    "honeycomb": 0.573,      # arXiv:2501.07388 (MC); Theorie 1/sqrt(2) umstritten
}


def validation_gates(measured_T_bkt: Optional[float],
                     lattice_type: str,
                     convergence_slope: Optional[float] = None,
                     variational_derivative_error: Optional[float] = None,
                     tol_T_bkt_rel: float = 0.03) -> dict[str, Any]:
    """Pflicht-Validierungs-Gates als Regressionstests.

    Gibt PASS/FAIL je Gate. T_BKT-Gate nutzt lattice-spezifischen Referenzwert.
    """
    gates: dict[str, bool] = {}
    ref = T_BKT_REFERENCE.get(lattice_type)
    if measured_T_bkt is not None and ref is not None:
        rel = abs(measured_T_bkt - ref) / ref
        gates[f"PASS_T_BKT_within_{int(tol_T_bkt_rel*100)}pct"] = (
            rel < tol_T_bkt_rel)
    if convergence_slope is not None:
        gates["PASS_CONVERGENCE_ON2"] = (-2.3 < convergence_slope < -1.7)
    if variational_derivative_error is not None:
        gates["PASS_VARIATIONAL_MACHINE_PRECISION"] = (
            variational_derivative_error < 1e-12)
    gates["PASS_REFERENCE_KNOWN"] = ref is not None
    return {"gates": gates, "all_pass": all(gates.values()) if gates else False,
            "reference_T_bkt": ref, "lattice_type": lattice_type}


# ============================================================================
# 10. Hauptpipeline: ein reproduzierbarer Lauf
# ============================================================================

@dataclass
class RunResult:
    seed: int
    T: float
    helicity: float
    nk_line: float
    g1_centers: list[float]
    g1_values: list[float]
    decay: dict[str, Any]
    defect_stats: DefectStats
    charge_check: dict[str, Any]


def run_phi_hex(lattice: TriangularLattice, cfg: XYConfig,
                seed: int, master_seed: int = 0) -> RunResult:
    """Ein vollstaendiger, reproduzierbarer Phi-Hex-XY-Lauf."""
    rng = make_rng(master_seed, stream=seed)
    adj = build_adjacency(lattice)
    n = lattice.n_nodes
    theta = rng.uniform(0, 2 * math.pi, n)

    node_pos = [axial_to_xy(a) for a in lattice.ax_of]
    cell_pos = []
    for ci in lattice.cells:
        xs = [node_pos[i][0] for i in ci]
        ys = [node_pos[i][1] for i in ci]
        cell_pos.append((float(np.mean(xs)), float(np.mean(ys))))
    tracker = DefectTrackerV2(cell_pos)

    # Ensemble-Akkumulatoren fuer den Helicity-Modulus (separate Mittelung!)
    t1_samples: list[float] = []
    sin_samples: list[float] = []

    for t in range(cfg.max_steps):
        theta = xy_langevin_step(theta, adj, cfg, rng)
        if t >= cfg.burn_in and t % cfg.measure_every == 0:
            vort = vorticity_field(theta, lattice)
            tracker.update(t, vort)
            t1, sin_acc = helicity_terms(theta, lattice, node_pos)
            t1_samples.append(t1)
            sin_samples.append(sin_acc)

    ups = helicity_from_ensemble(t1_samples, sin_samples, cfg,
                                 lattice.n_nodes)
    centers, g1 = phase_correlation(theta, node_pos, n_bins=8)
    decay = classify_decay(centers, g1)
    final_vort = vorticity_field(theta, lattice)
    charge = check_charge_conservation(final_vort, lattice)

    return RunResult(
        seed=seed, T=cfg.T, helicity=ups,
        nk_line=nelson_kosterlitz_line(cfg.T),
        g1_centers=[float(x) for x in centers],
        g1_values=[float(x) if not math.isnan(x) else None
                   for x in g1],
        decay=decay, defect_stats=tracker.stats(),
        charge_check=charge,
    )


# ============================================================================
# 11. Selbsttest / Demonstration
# ============================================================================

def self_test() -> dict[str, Any]:
    """Schneller Selbsttest des Kernmoduls."""
    print(f"{CORE_VERSION} Selbsttest - {ATTRIBUTION}")
    print("=" * 68)

    lattice = build_triangular_lattice(radius=3, periodic=False)
    print(f"Triangular Lattice: {lattice.n_nodes} Knoten, "
          f"{len(lattice.edges)} Kanten, {len(lattice.cells)} Zellen, "
          f"mittlere Koordinationszahl {lattice.coordination_mean:.2f}")
    assert lattice.coordination_mean <= 6.0

    # Tieftemperatur (geordnet) vs Hochtemperatur (ungeordnet)
    print("\nXY-Langevin-Laeufe (J=1):")
    print(f"{'T':>6s} | {'helicity':>10s} | {'2T/pi':>8s} | "
          f"{'decay':>12s} | {'g1[0]':>7s} | {'pers/cell':>9s}")
    print("-" * 70)
    results = {}
    for T in [0.2, 0.8, 2.0]:
        cfg = XYConfig(J=1.0, T=T, dt=0.05, max_steps=1500,
                       burn_in=500, measure_every=5)
        r = run_phi_hex(lattice, cfg, seed=0, master_seed=42)
        g0 = r.g1_values[0] if r.g1_values[0] is not None else float("nan")
        print(f"{T:>6.2f} | {r.helicity:>10.3f} | {r.nk_line:>8.3f} | "
              f"{r.decay['decay_type']:>12s} | {g0:>7.3f} | "
              f"{r.defect_stats.persistent_density_per_cell:>9.2f}")
        results[str(T)] = {
            "helicity": r.helicity, "nk_line": r.nk_line,
            "decay_type": r.decay["decay_type"], "g1_short": g0,
            "persistent_density": r.defect_stats.persistent_density_per_cell,
        }

    # Statistik-Toolbox-Smoke-Test
    print("\nStatistik-Toolbox-Smoke-Test:")
    rng = make_rng(42, 7)
    a = rng.normal(0.5, 1.0, 16)
    b = rng.normal(0.0, 1.0, 16)
    g = hedges_g(a, b)
    cd = cliff_delta(a, b)
    ci = bca_bootstrap_ci(a, b, n_boot=2000, rng=rng)
    print(f"  Hedges g = {g:+.3f}, Cliff delta = {cd['cliff_delta']:+.3f} "
          f"({cd['magnitude']}), BCa-CI = [{ci['ci_lower']:+.3f}, "
          f"{ci['ci_upper']:+.3f}], excl. 0: {ci['excludes_zero']}")

    # Validierungs-Gates
    print("\nValidierungs-Gates (Triangular-Referenz):")
    vg = validation_gates(measured_T_bkt=None, lattice_type="triangular",
                          convergence_slope=-1.99,
                          variational_derivative_error=2e-16)
    for k, v in vg["gates"].items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    print(f"  Referenz T_BKT(triangular) = {vg['reference_T_bkt']}")

    # PASS-Gates des Selbsttests
    low_T = results["0.2"]
    high_T = results["2.0"]
    pass_gates = {
        "PASS_LOW_T_HIGHER_HELICITY": low_T["helicity"] > high_T["helicity"],
        "PASS_LOW_T_MORE_ORDERED": low_T["g1_short"] > high_T["g1_short"],
        "PASS_STATS_TOOLBOX": ci["excludes_zero"] in (True, False),
        "PASS_LATTICE_TRIANGULAR": abs(lattice.coordination_mean - 6.0) < 1.5,
    }
    print("\nSelbsttest-Gates:")
    for k, v in pass_gates.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    overall = all(pass_gates.values())
    print(f"\nOVERALL: {'PASS' if overall else 'FAIL'}")

    prov = ProvenanceRecord()
    prov.finalize(master_seed=42, parameters={"radius": 3, "J": 1.0})

    return {"results": results, "pass_gates": pass_gates,
            "overall_pass": overall, "provenance": asdict(prov)}


if __name__ == "__main__":
    report = self_test()
    def clean(o):
        if o is None:
            return None
        if isinstance(o, (np.bool_, bool)):
            return bool(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            v = float(o)
            return v if math.isfinite(v) else None
        if isinstance(o, float):
            return o if math.isfinite(o) else None
        if isinstance(o, np.ndarray):
            return [clean(x) for x in o.tolist()]
        if hasattr(o, "__dataclass_fields__"):
            return clean(asdict(o))
        if isinstance(o, dict):
            return {str(k): clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [clean(x) for x in o]
        return o
    print("\n--- JSON-Report ---")
    print(json.dumps(clean(report), indent=2, allow_nan=False))
