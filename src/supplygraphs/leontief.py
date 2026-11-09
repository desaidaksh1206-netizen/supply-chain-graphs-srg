"""Supply chains as weighted walks: the Leontief inverse (paper Section 2--3).

For the direct-requirements matrix A (a_ij = $ of industry i's output used per
$ of industry j's output), the weighted-walk theorem gives

    (A^k)_ij = sum over walks i = v0 -> v1 -> ... -> vk = j of a_{v0 v1} ... a_{v(k-1) vk},

the requirement of j for i's output that travels through supply chains of
exactly k steps. When the spectral radius rho(A) < 1,

    L = (I - A)^-1 = I + A + A^2 + ...      (the Leontief inverse)

collects all of them.
"""

from __future__ import annotations

import numpy as np


def spectral_radius(a: np.ndarray) -> float:
    return float(np.max(np.abs(np.linalg.eigvals(a))))


def perron(a: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    """Perron root and normalised (sum 1) non-negative right and left eigenvectors."""
    ev, vr = np.linalg.eig(a)
    i = int(np.argmax(ev.real))
    rho = float(ev[i].real)
    right = np.abs(vr[:, i].real)
    evl, vl = np.linalg.eig(a.T)
    j = int(np.argmax(evl.real))
    left = np.abs(vl[:, j].real)
    return rho, right / right.sum(), left / left.sum()


def leontief_inverse(a: np.ndarray) -> np.ndarray:
    n = a.shape[0]
    return np.linalg.solve(np.eye(n) - a, np.eye(n))


def neumann_partial_sum(a: np.ndarray, k: int) -> np.ndarray:
    """I + A + ... + A^k."""
    n = a.shape[0]
    s, p = np.eye(n), np.eye(n)
    for _ in range(k):
        p = p @ a
        s = s + p
    return s


def walk_shares(a: np.ndarray, kmax: int = 12) -> np.ndarray:
    """s_k = 1'A^k 1 / 1'(L - I) 1 for k = 1..kmax: the fraction of all indirect
    and direct requirements that travel through supply chains of exactly k steps.
    The s_k over all k >= 1 sum to 1."""
    n = a.shape[0]
    total = (leontief_inverse(a) - np.eye(n)).sum()
    out, p = [], np.eye(n)
    for _ in range(kmax):
        p = p @ a
        out.append(p.sum() / total)
    return np.array(out)


def average_chain_length(a: np.ndarray) -> float:
    """sum_k k s_k = 1'(L - I) L 1 / 1'(L - I) 1, using sum_k k A^k = (L - I) L."""
    n = a.shape[0]
    L = leontief_inverse(a)
    return float(((L - np.eye(n)) @ L).sum() / (L - np.eye(n)).sum())


def output_multipliers(a: np.ndarray) -> np.ndarray:
    """Column sums of L: total output, across all industries, needed to deliver
    $1 of industry j's product to final users (Miller & Blair 2009)."""
    return leontief_inverse(a).sum(axis=0)
