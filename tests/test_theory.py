"""Numerical checks of every theorem and worked example the paper relies on.

These tests do not prove anything; the proofs are in the paper. They guard
against a typo in a formula or a bug in the code that implements it.
"""

import itertools

import numpy as np
import pytest

from supplygraphs import graphs as gr
from supplygraphs import leontief as lf

RNG = np.random.default_rng(20261020)


# ------------------------------------------------ Proposition prop:residual
@pytest.mark.parametrize("seed", range(5))
def test_residual_identity_on_random_graphs(seed):
    rng = np.random.default_rng(seed)
    a = gr.erdos_renyi(30, 0.3, rng)
    prof = gr.regularity_profile(a)
    # value of the minimum equals the variance decomposition
    at_means = gr.residual_norm_sq(a, prof["k_fit"], prof["lambda_fit"], prof["mu_fit"])
    assert at_means == pytest.approx(prof["residual_sq"], rel=1e-10)
    # and no nearby (k, lambda, mu) does better
    for dk, dl, dm in rng.normal(scale=0.5, size=(50, 3)):
        other = gr.residual_norm_sq(
            a, prof["k_fit"] + dk, prof["lambda_fit"] + dl, prof["mu_fit"] + dm
        )
        assert other >= at_means - 1e-9


def test_delta_between_zero_and_one():
    for p in (0.1, 0.5, 0.9):
        d = gr.regularity_profile(gr.erdos_renyi(40, p, RNG))["delta"]
        assert 0 <= d <= 1


# ------------------------------------------------ Corollary cor:delta-zero
@pytest.mark.parametrize(
    "a, params",
    [
        (gr.petersen_graph(), (3, 0, 1)),
        (gr.paley_graph(13), (6, 2, 3)),
        (gr.hamming_graph(2, 3), (4, 1, 2)),
        (gr.hamming_graph(2, 10), (18, 8, 2)),
        (gr.hamming_graph(2, 20), (38, 18, 2)),
    ],
)
def test_strongly_regular_graphs_have_delta_zero(a, params):
    prof = gr.regularity_profile(a)
    assert prof["delta"] == pytest.approx(0, abs=1e-12)
    assert (prof["k_fit"], prof["lambda_fit"], prof["mu_fit"]) == params


def test_hamming_3_is_not_strongly_regular():
    assert gr.regularity_profile(gr.hamming_graph(3, 3))["delta"] > 0


def test_hamming_2_20_spectrum_and_cliques():
    h = gr.hamming_graph(2, 20)
    ev = np.round(gr.spectrum(h)).astype(int)
    vals, counts = np.unique(ev, return_counts=True)
    assert dict(zip(vals.tolist(), counts.tolist())) == {-2: 361, 18: 38, 38: 1}
    prof = gr.structure_profile(h, clique_max_density=0.2)
    assert prof["clique_number"] == 20


# ------------------------------------------------ Example ex:bowtie
def test_bowtie_example_numbers():
    a = gr.bowtie_graph()
    expected_a2 = np.array(
        [[4, 1, 1, 1, 1], [1, 2, 1, 1, 1], [1, 1, 2, 1, 1], [1, 1, 1, 2, 1], [1, 1, 1, 1, 2]]
    )
    assert (a @ a == expected_a2).all()
    prof = gr.regularity_profile(a)
    assert prof["lambda_fit"] == 1 and prof["mu_fit"] == 1
    assert prof["k_fit"] == pytest.approx(12 / 5)
    assert prof["term_edge"] == 0 and prof["term_nonedge"] == 0
    assert prof["residual_sq"] == pytest.approx(16 / 5)


# ------------------------------------------------ Proposition prop:er
def test_er_common_neighbours_binomial_and_independent_of_adjacency():
    n, p, reps = 60, 0.3, 300
    lam, mu = [], []
    for _ in range(reps):
        _, l_, m_ = gr.common_neighbour_lists(gr.erdos_renyi(n, p, RNG))
        lam.append(l_)
        mu.append(m_)
    lam, mu = np.concatenate(lam), np.concatenate(mu)
    mean, var = (n - 2) * p**2, (n - 2) * p**2 * (1 - p**2)
    for x in (lam, mu):
        assert x.mean() == pytest.approx(mean, rel=0.02)
        assert x.var() == pytest.approx(var, rel=0.05)


# ------------------------------------------------ degree-preserving null
def test_degree_preserving_keeps_degrees_and_simplicity():
    a = gr.erdos_renyi(50, 0.2, RNG)
    b = gr.degree_preserving(a, RNG)
    assert (a.sum(0) == b.sum(0)).all()
    assert (b == b.T).all() and np.diag(b).sum() == 0 and set(np.unique(b)) <= {0, 1}
    assert (a != b).any()


# ------------------------------------------------ threshold graphs
def test_threshold_graphs_are_nested_and_symmetric():
    a = RNG.random((30, 30)) * 0.1
    w = gr.link_weights(a)
    assert (w == w.T).all() and np.diag(w).sum() == 0
    prev = None
    for theta in 10 ** np.linspace(-1, -4, 13):
        g = gr.threshold_adjacency(w, theta)
        if prev is not None:
            assert (g >= prev).all()
        prev = g


def test_square_is_exact():
    a = gr.erdos_renyi(120, 0.4, RNG)
    assert (gr.square(a) == a.astype(np.int64) @ a.astype(np.int64)).all()


# ------------------------------------------------ weighted walks and the Leontief inverse
def _random_requirements(n, scale=0.15, seed=0):
    rng = np.random.default_rng(seed)
    return rng.random((n, n)) * scale / n * 3


def test_weighted_walk_theorem_by_enumeration():
    """(A^k)_ij equals the sum over all walks of length k of the product of weights."""
    a = _random_requirements(4, seed=1)
    n, k = 4, 3
    brute = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            for mid in itertools.product(range(n), repeat=k - 1):
                path = (i, *mid, j)
                brute[i, j] += np.prod([a[path[t], path[t + 1]] for t in range(k)])
    assert np.allclose(np.linalg.matrix_power(a, k), brute)


def test_neumann_series_converges_to_leontief_inverse():
    a = _random_requirements(30, seed=2)
    assert lf.spectral_radius(a) < 1
    L = lf.leontief_inverse(a)
    assert np.allclose(L @ (np.eye(30) - a), np.eye(30))
    assert np.abs(lf.neumann_partial_sum(a, 60) - L).max() < 1e-10
    assert (L >= 0).all()


def test_column_sum_bound_on_spectral_radius():
    a = _random_requirements(25, seed=3)
    assert lf.spectral_radius(a) <= a.sum(axis=0).max() + 1e-12


def test_walk_shares_and_average_length():
    a = _random_requirements(20, seed=4)
    s = lf.walk_shares(a, 200)
    assert s.sum() == pytest.approx(1, abs=1e-10)
    assert lf.average_chain_length(a) == pytest.approx((np.arange(1, 201) * s).sum(), rel=1e-9)


def test_perron_vectors():
    a = _random_requirements(15, seed=5)
    rho, right, left = lf.perron(a)
    assert rho == pytest.approx(lf.spectral_radius(a))
    assert np.allclose(a @ right, rho * right)
    assert np.allclose(left @ a, rho * left)
    assert (right >= 0).all() and (left >= 0).all()


def test_output_multipliers_at_least_one():
    a = _random_requirements(15, seed=6)
    assert (lf.output_multipliers(a) >= 1).all()
