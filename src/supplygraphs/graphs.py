"""Threshold graphs, the distance from strong regularity, and benchmark graphs.

The central object is :func:`regularity_profile`, which implements
Proposition `prop:residual` of the paper: for a graph with adjacency matrix A,
the least-squares fit of the strongly regular identity

    A^2 = k I + lambda A + mu (J - I - A)

has minimum squared residual  n Var(D) + 2m Var(L) + 2m' Var(M),  where D is the
list of degrees, L the common-neighbour counts over ordered adjacent pairs, and
M the counts over ordered distinct non-adjacent pairs.
"""

from __future__ import annotations

import itertools

import networkx as nx
import numpy as np


# ---------------------------------------------------------------- construction
def link_weights(a: np.ndarray) -> np.ndarray:
    """w_ij = max(a_ij, a_ji) for i != j, and 0 on the diagonal.

    w_ij >= theta exactly when one of the two industries buys at least theta
    dollars from the other per dollar of its own output."""
    w = np.maximum(a, a.T).astype(float)
    np.fill_diagonal(w, 0.0)
    return w


def threshold_adjacency(w: np.ndarray, theta: float) -> np.ndarray:
    """A_theta: i ~ j iff i != j and w_ij >= theta (Definition `def:threshold`)."""
    a = (w >= theta).astype(np.int64)
    np.fill_diagonal(a, 0)
    return a


def from_edges(n: int, edges) -> np.ndarray:
    a = np.zeros((n, n), dtype=np.int64)
    for u, v in edges:
        a[u, v] = a[v, u] = 1
    return a


def hamming_graph(d: int, q: int) -> np.ndarray:
    """H(d, q): words of length d over q symbols, adjacent iff Hamming distance 1."""
    words = list(itertools.product(range(q), repeat=d))
    w = np.array(words)
    dist = (w[:, None, :] != w[None, :, :]).sum(axis=2)
    return (dist == 1).astype(np.int64)


def petersen_graph() -> np.ndarray:
    """Kneser graph K(5,2): 2-subsets of a 5-set, adjacent iff disjoint."""
    pairs = list(itertools.combinations(range(5), 2))
    return np.array([[int(not set(a) & set(b)) for b in pairs] for a in pairs])


def paley_graph(p: int) -> np.ndarray:
    """Paley graph on Z_p for a prime p = 1 mod 4."""
    squares = {(x * x) % p for x in range(1, p)}
    return np.array([[int(i != j and (i - j) % p in squares) for j in range(p)] for i in range(p)])


def bowtie_graph() -> np.ndarray:
    """Two triangles sharing v0, ordered v0..v4 as in Example `ex:bowtie`."""
    return from_edges(5, [(0, 1), (0, 2), (1, 2), (0, 3), (0, 4), (3, 4)])


# ------------------------------------------------------------ the instrument
def square(a: np.ndarray) -> np.ndarray:
    """A^2 for a 0/1 matrix, computed in floating point (fast BLAS) and rounded
    back to exact integers; exact because every entry is at most n < 2^53."""
    af = np.asarray(a, dtype=float)
    return np.rint(af @ af).astype(np.int64)


def common_neighbour_lists(a: np.ndarray):
    """Return (D, L, M): degrees, counts on ordered edges, counts on ordered non-edges."""
    a = np.asarray(a)
    n = a.shape[0]
    a2 = square(a)
    off = ~np.eye(n, dtype=bool)
    adj = (a == 1) & off
    non = (a == 0) & off
    return np.diag(a2).copy(), a2[adj], a2[non]


def regularity_profile(a: np.ndarray) -> dict:
    """Fitted (k, lambda, mu), the three residual terms and delta(G).

    delta(G) = min ||A^2 - kI - lambda A - mu Abar||_F / ||A^2||_F,
    which is 0 exactly for strongly regular graphs (Corollary `cor:delta-zero`).
    """
    a = np.asarray(a)
    n = a.shape[0]
    d, lam, mu = common_neighbour_lists(a)
    out = {"n": n, "m": int(len(lam) // 2), "m_bar": int(len(mu) // 2)}
    if len(lam) == 0 or len(mu) == 0:
        out.update({k: np.nan for k in (
            "k_fit", "lambda_fit", "mu_fit", "gap", "term_degree", "term_edge",
            "term_nonedge", "residual_sq", "delta", "frac_degree", "frac_edge",
            "frac_nonedge", "sd_lambda", "sd_mu")})
        return out
    t_d = n * d.var()
    t_l = len(lam) * lam.var()
    t_m = len(mu) * mu.var()
    res = t_d + t_l + t_m
    a2 = square(a)
    out.update(
        k_fit=float(d.mean()),
        lambda_fit=float(lam.mean()),
        mu_fit=float(mu.mean()),
        gap=float(lam.mean() - mu.mean()),
        sd_lambda=float(lam.std()),
        sd_mu=float(mu.std()),
        term_degree=float(t_d),
        term_edge=float(t_l),
        term_nonedge=float(t_m),
        residual_sq=float(res),
        delta=float(np.sqrt(res) / np.linalg.norm(a2)),
        frac_degree=float(t_d / res) if res > 0 else np.nan,
        frac_edge=float(t_l / res) if res > 0 else np.nan,
        frac_nonedge=float(t_m / res) if res > 0 else np.nan,
    )
    return out


def residual_norm_sq(a: np.ndarray, k: float, lam: float, mu: float) -> float:
    """||A^2 - kI - lam A - mu (J - I - A)||_F^2 for arbitrary (k, lam, mu)."""
    n = a.shape[0]
    i = np.eye(n)
    abar = np.ones((n, n)) - i - a
    r = a @ a - k * i - lam * a - mu * abar
    return float((r * r).sum())


# ------------------------------------------------------------ null models
def erdos_renyi(n: int, p: float, rng: np.random.Generator) -> np.ndarray:
    """G(n, p): each of the C(n,2) edges independently with probability p."""
    upper = np.triu(rng.random((n, n)) < p, k=1)
    return (upper | upper.T).astype(np.int64)


def degree_preserving(a: np.ndarray, rng: np.random.Generator, swaps_per_edge: int = 10) -> np.ndarray:
    """Random graph with the same degree sequence, by double-edge swaps.

    Repeatedly pick edges {a,b}, {c,d} and replace them with {a,d}, {c,b} when
    that creates no loop or multi-edge. Every degree is unchanged by each swap.
    """
    adj = np.array(a, dtype=bool, copy=True)
    iu = np.triu_indices_from(adj, k=1)
    mask = adj[iu]
    edges = np.column_stack([iu[0][mask], iu[1][mask]])
    m = len(edges)
    if m < 2:
        return adj.astype(np.int64)
    tries = swaps_per_edge * m
    picks = rng.integers(0, m, size=(tries, 2))
    flips = rng.random(tries) < 0.5
    for (i, j), flip in zip(picks, flips):
        if i == j:
            continue
        a_, b_ = edges[i]
        c_, d_ = edges[j]
        if flip:
            c_, d_ = d_, c_
        if a_ == d_ or c_ == b_ or a_ == c_ or b_ == d_:
            continue
        if adj[a_, d_] or adj[c_, b_]:
            continue
        adj[a_, b_] = adj[b_, a_] = adj[c_, d_] = adj[d_, c_] = False
        adj[a_, d_] = adj[d_, a_] = adj[c_, b_] = adj[b_, c_] = True
        edges[i] = (a_, d_)
        edges[j] = (c_, b_)
    return adj.astype(np.int64)


# ------------------------------------------------------------ descriptors
def spectrum(a: np.ndarray) -> np.ndarray:
    return np.sort(np.linalg.eigvalsh(np.asarray(a, dtype=float)))[::-1]


def n_distinct(values: np.ndarray, tol: float = 1e-6) -> int:
    v = np.sort(np.asarray(values))
    return int(1 + (np.diff(v) > tol).sum()) if len(v) else 0


def structure_profile(a: np.ndarray, clique_max_density: float = 0.15) -> dict:
    """Degrees, components, clique number and spectral quantities of one graph."""
    a = np.asarray(a)
    n = a.shape[0]
    g = nx.from_numpy_array(a)
    comps = sorted((len(c) for c in nx.connected_components(g)), reverse=True)
    deg = a.sum(axis=1)
    m = int(deg.sum() // 2)
    density = m / (n * (n - 1) / 2)
    ev = spectrum(a)
    if m == 0:
        omega = 1
    elif density <= clique_max_density:
        omega = max(len(c) for c in nx.find_cliques(g))
    else:
        omega = np.nan  # exact search is expensive for dense graphs; not reported
    return {
        "m": m,
        "density": density,
        "mean_degree": float(deg.mean()),
        "sd_degree": float(deg.std()),
        "max_degree": int(deg.max()),
        "n_isolated": int((deg == 0).sum()),
        "n_components": len(comps),
        "largest_component": comps[0],
        "clique_number": omega,
        "lambda1": float(ev[0]),
        "lambda2": float(ev[1]),
        "spectral_gap": float(ev[0] - ev[1]),
        "lambda_min": float(ev[-1]),
        "n_distinct_eigenvalues": n_distinct(ev),
    }
