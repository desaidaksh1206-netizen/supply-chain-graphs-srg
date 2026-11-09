# Methods: the mathematics, step by step

This document explains every step of the analysis in the order the code runs it. For each step it covers what is computed, why, and which theorem makes the step legitimate. The paper (`paper/sec-supply-chain.tex`) gives the same material more briefly. Proofs are given here in full. Anything that is only computed is labelled as such.

**Notation.** *N* = 400 industries. *A* is the direct-requirements matrix (§1). 1 is the all-ones vector. In §§4–6, *A* also denotes the 0/1 adjacency matrix of a graph; the context makes clear which is meant. *J* is the all-ones matrix and *Ā* = *J* − *I* − *A*.

**Vertex ordering.** Industries are numbered 1–400 in order of their BEA industry code (`data/processed/industries.csv`, column `index`). Every matrix uses this ordering, except where a figure regroups by sector to show block structure. Regrouping replaces *A* by *PAP*ᵀ for a permutation matrix *P*, which describes the same graph.

---

## 1. From BEA tables to the requirements matrix

*Code:* `data.direct_requirements`

The BEA publishes two tables per year, in millions of dollars:

- **Use** *U*, commodity × industry: how much of each commodity each industry buys as an intermediate input.
- **Make** *V*, industry × commodity: how much of each commodity each industry produces.

With industry outputs *g* = *V*1 and commodity outputs *q* = *V*ᵀ1:

- *B* = *U* diag(*g*)⁻¹: commodity inputs per dollar of industry output;
- *D* = *V* diag(*q*)⁻¹: market shares, the fraction of commodity *c* made by industry *i* (a zero column for a commodity no industry produces, where *q* = 0);
- **A = DB**: industry × industry direct requirements.

This is the standard "industry-technology" construction (Miller & Blair 2009, Ch. 5; Horowitz & Planting 2009). Reading *DB* entry by entry: industry *j* buys commodity *c* (*B*), and the dollars spent on *c* are attributed to the industries that make *c* in proportion to their market shares (*D*). So *a*ᵢⱼ is the dollars of industry *i*'s output embodied in one dollar of *j*'s output, as a direct input.

Commodities produced by no industry (noncomparable imports, and used and secondhand goods) have *q* = 0. Their columns of *D* are set to zero, so they act like imports or value added: inputs from outside the network.

## 2. Supply chains are weighted walks

*Code:* `leontief.py`; *tests:* `test_weighted_walk_theorem_by_enumeration`, `test_neumann_series_converges_to_leontief_inverse`

Treat *A* as a weighted directed graph with an arc *i* → *j* of weight *a*ᵢⱼ: from supplier to buyer. A supply chain of length *k* is a walk *i* = *v*₀ → *v*₁ → ⋯ → *v*ₖ = *j*, and its weight is *a*_{v₀v₁} ⋯ *a*_{vₖ₋₁vₖ}.

**Theorem (weighted walks).** (*A*ᵏ)ᵢⱼ is the sum of the weights of all walks of length *k* from *i* to *j*.

*Proof.* By induction on *k*. For *k* = 1 the only walk is the arc itself. A walk of length *k* + 1 is a walk of length *k* to some ℓ followed by the arc ℓ → *j*, so the total weight is Σ_ℓ (*A*ᵏ)ᵢℓ *a*ℓⱼ = (*A*ᵏ⁺¹)ᵢⱼ. ∎

With all weights equal to 1 this is the article's walk-counting theorem. **Economic reading:** making $1 of *j* uses *a*ℓⱼ of each supplier ℓ, which in turn uses *a*ᵢℓ of *i* per dollar, and so on back along the chain. So (*A*ᵏ)ᵢⱼ is the dollars of *i*'s output that $1 of *j* needs through chains of exactly *k* − 1 intermediaries.

**Theorem (total requirements).** If *A* ≥ 0 and ρ(*A*) < 1, then *I* − *A* is invertible, *L* = (*I* − *A*)⁻¹ = *I* + *A* + *A*² + ⋯, and *L* ≥ 0 entrywise.

*Proof.* ρ(*A*) < 1 implies *A*ᵏ → 0 (Horn & Johnson 2013, §5.6). The identity (*I* − *A*)(*I* + ⋯ + *A*ᵏ) = *I* − *A*ᵏ⁺¹ holds for every *k*. If (*I* − *A*)*x* = 0 then *x* = *A*ᵏ⁺¹*x* → 0, so *I* − *A* is invertible. Hence the partial sums equal (*I* − *A*)⁻¹(*I* − *A*ᵏ⁺¹) → (*I* − *A*)⁻¹. Each partial sum is ≥ 0, so the limit is too. ∎

*L* is the **Leontief inverse** (Leontief 1936). Its column sums are **output multipliers**: the total output across the economy needed to deliver $1 of *j*'s product to final users. Its row sums are **forward linkages**: how much of *i*'s output is needed when every industry delivers $1 to final users.

**Proposition (column-sum bound).** For any real square matrix, ρ(*A*) ≤ max_j Σᵢ |*a*ᵢⱼ|.

*Proof.* If *Av* = λ*v* with *v* ≠ 0, then |λ|‖*v*‖₁ = ‖*Av*‖₁ ≤ Σⱼ |*v*ⱼ| Σᵢ |*a*ᵢⱼ| ≤ (maxⱼ Σᵢ |*a*ᵢⱼ|)‖*v*‖₁. ∎

Column *j* of *A* sums to the share of *j*'s output spent on intermediate inputs, so one expects every column sum to be below 1, which would give ρ(*A*) < 1 at once. **On our data that fails for one industry.** State and local government passenger transit has negative value added in the tables: its intermediate purchases are 1.59 times the value of its output. So the bound cannot be used, and the code computes ρ(*A*) directly: 0.489 for 2017 and 0.513 for 2012. The test suite checks the bound itself on random matrices.

**Proposition (average chain length).** Let *s*ₖ = 1ᵀ*A*ᵏ1 / 1ᵀ(*L* − *I*)1. Then Σₖ *s*ₖ = 1 and Σₖ *k s*ₖ = 1ᵀ(*L* − *I*)*L*1 / 1ᵀ(*L* − *I*)1.

*Proof.* *L* − *I* = Σ_{k≥1} *A*ᵏ gives the first claim. For the second, Σ_{k≥1} *k A*ᵏ = Σ_{k≥1} Σ_{j≥k} *A*ʲ = Σ_{k≥1} *A*ᵏ*L* = (*L* − *I*)*L*, because Σ_{j≥k} *A*ʲ = *A*ᵏ*L*. ∎

*s*ₖ is the share of all requirements that travel along chains of exactly *k* steps. Σ *k s*ₖ is the dollar-weighted average chain length, an economy-wide version of the "average propagation length" of Dietzenbacher, Romero & Bosma (2005).

**Perron–Frobenius.** For a non-negative matrix, ρ(*A*) is itself an eigenvalue, with a non-negative eigenvector (Horn & Johnson 2013, §8.3). The largest eigenvalue sets the rate at which *A*ᵏ decays, which is why successive shares *s*ₖ₊₁/*s*ₖ approach ρ(*A*).

## 3. From requirements to a graph

*Code:* `graphs.link_weights`, `graphs.threshold_adjacency`

Set *w*ᵢⱼ = max(*a*ᵢⱼ, *a*ⱼᵢ) for *i* ≠ *j*. The threshold graph *G*_θ joins *i* and *j* when *w*ᵢⱼ ≥ θ: one of the two buys at least θ dollars from the other per dollar of its own output. Taking the maximum makes the graph undirected, as the article's graphs are. Direction is kept in §2, which uses *A* itself.

Lowering θ only adds edges. Because the *a*ᵢⱼ span several orders of magnitude, the 49 thresholds are spaced evenly on a log scale from 10⁻⁴ to 10⁻¹.

## 4. The instrument: distance from strong regularity

*Code:* `graphs.regularity_profile`; *tests:* `test_residual_identity_on_random_graphs`, `test_strongly_regular_graphs_have_delta_zero`, `test_bowtie_example_numbers`

For a 0/1 adjacency matrix, (*A*²)ᵤᵥ = *c*(*u*,*v*) is the number of common neighbours, and (*A*²)ᵤᵤ = deg *u*. A graph that is neither complete nor edgeless is strongly regular exactly when *A*² = *kI* + λ*A* + μ*Ā*. Since *I*, *A* and *Ā* have disjoint supports covering every entry, this says three things: every degree is *k*, every common-neighbour count on an edge is λ, and every count on a non-edge is μ. In general these quantities form three lists: 𝒟 (degrees), ℒ (counts on the 2*m* ordered edges) and ℳ (counts on the 2*m̄* ordered non-edges).

**Proposition.** min over *k*, λ, μ of ‖*A*² − *kI* − λ*A* − μ*Ā*‖²_F = *n*Var(𝒟) + 2*m*Var(ℒ) + 2*m̄*Var(ℳ), attained exactly at the three means.

*Proof.* The disjoint supports split the squared norm into three sums, each involving one unknown. For any list, Σ(*x*ᵢ − *a*)² = ℓVar(*x*) + ℓ(*x̄* − *a*)², which is minimised at *a* = *x̄*. ∎

**δ(*G*)** = (minimum residual norm) / ‖*A*²‖_F lies in [0, 1], and δ = 0 exactly when *G* is strongly regular. The three shares of the squared residual tell which part of the definition fails. **The bow-tie** (two triangles sharing a vertex) has λ = μ = 1 exactly but degrees 4, 2, 2, 2, 2, so only the degree term detects that it is not strongly regular.

`graphs.square` computes *A*² with floating-point BLAS and rounds back to integers. This is exact because every entry is at most *n* < 2⁵³; `test_square_is_exact` checks it.

## 5. The random benchmarks

**Erdős–Rényi** *G*(*n*, *p*), with *p* equal to the density of *G*_θ. **Proposition:** *c*(*u*,*v*) ~ Bin(*n* − 2, *p*²), independent of whether *u* ~ *v*. *Proof:* *w* is a common neighbour iff both *uw* and *vw* are edges (probability *p*²). These events are independent across *w* and never involve the edge *uv*. ∎ So in a random graph λ̄ ≈ μ̄.

**Degree-preserving rewiring** (Maslov & Sneppen 2002): repeatedly swap {*a*,*b*}, {*c*,*d*} → {*a*,*d*}, {*c*,*b*} when no loop or repeated edge results, with 5 attempts per edge. Every degree is preserved. This separates "unequal degrees" from "where the edges are". It is a computation; nothing about it is proved here.

At each θ: 200 Erdős–Rényi samples and 20 rewirings (50 at θ\*). The code reports the mean and central 95% range of δ and of λ̄ − μ̄. These ranges describe the benchmarks; they are not confidence intervals for the industry network, which is one fixed graph. Seed: `analysis.SEED = 20261027`.

## 6. The strongly regular benchmark: H(2,20)

Vertices: pairs (*a*, *b*) with *a*, *b* ∈ {0, …, 19}, adjacent when they differ in exactly one coordinate.

- *k* = 19 + 19 = 38.
- λ: adjacent (*a*,*b*) and (*a*′,*b*) share exactly the pairs (*a*″, *b*) with *a*″ ∉ {*a*, *a*′}, so λ = 18.
- μ: non-adjacent (*a*,*b*) and (*a*′,*b*′) with *a* ≠ *a*′ and *b* ≠ *b*′ share exactly (*a*, *b*′) and (*a*′, *b*), so μ = 2.
- Eigenvalues: 38, and the roots of *x*² − (λ − μ)*x* + (μ − *k*) = *x*² − 16*x* − 36 = (*x* − 18)(*x* + 2) (Brouwer & Haemers 2012). With multiplicities *f* and *g*: *f* + *g* = 399 and 38 + 18*f* − 2*g* = 0, so *f* = 38 and *g* = 361.
- Clique number 20: each "line" {(*a*, *b*) : *b* fixed} (or *a* fixed) is a clique. If three vertices are pairwise adjacent and two of them differ in the first coordinate, the third must agree with both in the second coordinate, so every clique of three or more lies in a line.

*N* = 400 was chosen so that this benchmark has exactly the order of the industry network. θ\* is the threshold at which *G*_θ has exactly the density of *H*(2,20), 38/399, and so mean degree 38 (`analysis.theta_for_density`).

## 7. Robustness

- **2012 at the same detail.** Same 400 industry codes, 2012 tables, full sweep (100 ER and 10 rewirings per threshold) plus the θ\* comparison.
- **Every year 2012–2023 at the 71-industry summary level.** Here the threshold matches the density of *H*(2,20). The code computes ρ(*A*), the average chain length, *s*₁, the gap and δ, with 200 ER samples and 50 rewirings per year.
- **Consistency check** (`test_two_levels_agree_on_spectral_radius`). For 2017, ρ(*A*) at the 71-industry level (0.4889) agrees with the 400-industry level (0.4886) to 0.1%. This is not a theorem. It indicates the two levels were built the same way.

## 8. Glossary for readers new to input-output economics

- **Industry vs commodity:** BEA classifies both what is made (commodities) and who makes it (industries). Most industries make mainly one commodity.
- **Intermediate input:** something bought to make something else, as opposed to final demand (households, investment, government, exports).
- **Value added:** what is left of a dollar of output after paying for intermediate inputs: wages, profits, taxes.
- **Direct requirements *A*:** inputs needed *directly* per dollar of output.
- **Total requirements *L*:** direct plus all indirect inputs, along supply chains of every length.
- **Output multiplier:** column sum of *L*: total economy-wide output generated by $1 of final demand for one industry's product.
