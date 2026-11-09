"""The threshold sweep and every number reported in the paper.

All randomness flows from one seed (SEED), so a rerun reproduces every value.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import graphs as gr

SEED = 20261027
# Thresholds in dollars of input per dollar of output, evenly spaced on a log scale
# from 1e-4 (0.01 cents) to 1e-1 (10 cents).
THETAS = np.round(10 ** np.linspace(-4, -1, 49), 8)
N_ER = 200          # Erdos--Renyi samples per threshold
N_DP = 20           # degree-preserving samples per threshold
SWAPS_PER_EDGE = 5
TARGET_DEGREE = 38  # degree of H(2,20), the strongly regular benchmark on 400 vertices
TARGET_DENSITY = 38 / 399


def _band(values) -> tuple[float, float, float]:
    v = np.asarray([x for x in values if np.isfinite(x)])
    if len(v) == 0:
        return np.nan, np.nan, np.nan
    return float(v.mean()), float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))


def null_bands(a: np.ndarray, rng: np.random.Generator, n_er: int = N_ER, n_dp: int = N_DP) -> dict:
    """delta and lambda-mu gap for ER graphs of equal density and for
    degree-preserving rewirings of the same graph."""
    n = a.shape[0]
    m = a.sum() // 2
    p = m / (n * (n - 1) / 2)
    out = {}
    er = [gr.regularity_profile(gr.erdos_renyi(n, p, rng)) for _ in range(n_er)]
    out["er_delta"], out["er_delta_lo"], out["er_delta_hi"] = _band(x["delta"] for x in er)
    out["er_gap"], out["er_gap_lo"], out["er_gap_hi"] = _band(x["gap"] for x in er)
    if n_dp:
        dp = [gr.regularity_profile(gr.degree_preserving(a, rng, SWAPS_PER_EDGE)) for _ in range(n_dp)]
        out["dp_delta"], out["dp_delta_lo"], out["dp_delta_hi"] = _band(x["delta"] for x in dp)
        out["dp_gap"], out["dp_gap_lo"], out["dp_gap_hi"] = _band(x["gap"] for x in dp)
    return out


def sweep(w: np.ndarray, thetas=THETAS, seed: int = SEED, n_er: int = N_ER, n_dp: int = N_DP) -> pd.DataFrame:
    """One row per threshold: structure, regularity profile and null bands."""
    rng = np.random.default_rng(seed)
    rows = []
    for theta in thetas:
        a = gr.threshold_adjacency(w, theta)
        row = {"theta": float(theta)}
        row.update(gr.structure_profile(a))
        prof = gr.regularity_profile(a)
        row.update({k: v for k, v in prof.items() if k not in ("n", "m")})
        if 0 < row["m"] < a.shape[0] * (a.shape[0] - 1) // 2:
            row.update(null_bands(a, rng, n_er, n_dp))
        rows.append(row)
    return pd.DataFrame(rows)


def theta_for_density(w: np.ndarray, density: float) -> float:
    """The threshold whose graph has density closest to `density`: the value of
    the k-th largest off-diagonal weight, k = round(density * C(n,2))."""
    n = w.shape[0]
    vals = np.sort(w[np.triu_indices(n, 1)])[::-1]
    k = int(round(density * len(vals)))
    return float(vals[k - 1])


def group_profile(a: np.ndarray, groups) -> dict:
    """How far edges, and the lambda-mu gap, concentrate inside groups."""
    s = np.asarray(groups)
    n = len(s)
    same = (s[:, None] == s[None, :]) & ~np.eye(n, dtype=bool)
    adj = (a == 1) & ~np.eye(n, dtype=bool)
    non = (a == 0) & ~np.eye(n, dtype=bool)
    a2 = gr.square(a)
    counts = pd.Series(s).value_counts()
    return {
        "frac_edges_within": float(adj[same].sum() / adj.sum()),
        "frac_pairs_within": float((counts * (counts - 1)).sum() / (n * (n - 1))),
        "lambda_within": float(a2[adj & same].mean()) if (adj & same).any() else np.nan,
        "lambda_across": float(a2[adj & ~same].mean()) if (adj & ~same).any() else np.nan,
        "mu_within": float(a2[non & same].mean()) if (non & same).any() else np.nan,
        "mu_across": float(a2[non & ~same].mean()) if (non & ~same).any() else np.nan,
        "n_edges_within": int(adj[same].sum() // 2),
        "n_edges_across": int(adj[~same].sum() // 2),
    }


def at_density(w: np.ndarray, density: float, rng: np.random.Generator,
               n_er: int = N_ER, n_dp: int = 50) -> dict:
    """Regularity profile and null bands at the density-matched threshold."""
    theta = theta_for_density(w, density)
    a = gr.threshold_adjacency(w, theta)
    out = {"theta": theta}
    out.update(gr.structure_profile(a))
    out.update({k: v for k, v in gr.regularity_profile(a).items() if k not in ("n", "m")})
    out.update(null_bands(a, rng, n_er, n_dp))
    return out


def ranges(thetas, mask):
    """Contiguous theta intervals where mask holds, as [(lo, hi), ...]."""
    out, start, prev = [], None, None
    for t, m in zip(thetas, mask):
        if m and start is None:
            start = t
        if not m and start is not None:
            out.append((start, prev))
            start = None
        prev = t
    if start is not None:
        out.append((start, prev))
    return [(float(a), float(b)) for a, b in out]
