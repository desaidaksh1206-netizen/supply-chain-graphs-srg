# The LaTeX section: structure and how to use it

`paper/sec-supply-chain.tex` is the application section of the article. It replaces the earlier finance section. It follows the advisor's requirements for the next revision, numbered as in the advisor's email:

| Requirement | How the section meets it |
|---|---|
| 1. No code in the body | Only definitions, theorems, proofs, matrices, figures and tables. All code is in this repository. Each of the six figures is numbered, captioned and discussed in the text. |
| 2. Adjacency matrices early, ordering stated | The vertex ordering (by BEA industry code) is stated before the first matrix. The bow-tie example shows its ordering, *A* and *A*². Figure `fig:adjacency` shows *A*<sub>θ\*</sub> before any matrix property of the network is discussed, and explains the regrouped version as *PAP*ᵀ. |
| 3. Scholarly references | The mathematics is cited to books and peer-reviewed papers (Horn & Johnson, Miller & Blair, Brouwer & Haemers, …). The BEA tables are cited as the data source, together with the peer-reviewed paper describing `useeior`. |
| 4. Proof vs computation | Four results about walks and the Leontief inverse are proved in the text: Theorems `thm:weighted-walks` and `thm:neumann`, Propositions `prop:colsum` and `prop:chain-length`. So are the graph results: Propositions `prop:residual` and `prop:er` and Corollary `cor:delta-zero`. Two matrix-analysis facts are cited. Every number about the economy is introduced as a computation on the stated tables. Table `tab:supply-status` classifies each link. |
| 5. Johnson graphs, Hamming graphs and association schemes | The section does not redefine these families; it relies on Section `sec:johnson-hamming`. Where it uses the Hamming graph *H*(2,20), it states the vertex set, the adjacency condition, the parameters (400, 38, 18, 2), the eigenvalues with their multiplicities and the clique number, citing Brouwer & Haemers for the eigenvalues; [`METHODS.md`](METHODS.md) §6 derives them. Association schemes are not used. As the advisor asked, the definition of *H*(2,20) cites their thesis (Joshi 2023, Example 2.2.2), which defines the Hamming scheme and the Hamming graph. |
| 6. Applications | Agreed with the advisor at the November 2026 Zoom meeting: this section is the article's only application. |
| 7. LaTeX article format | Numbered definitions, theorems, equations, figures and tables, with a BibTeX bibliography. The title, abstract and introduction belong to `main.tex`. |

**How it ties to the rest of the article.** The weighted walks theorem is the article's walks-and-matrix-powers theorem with weights. It reduces to the unweighted theorem when every weight is 1. The Leontief inverse is then a direct application of it, not an analogy. The second half reuses the definition of strong regularity, the Hamming graphs, and the spectral facts from earlier sections.

## Files

| File | Purpose |
|---|---|
| `sec-supply-chain.tex` | The section. It contains no hand-typed numbers. |
| `generated/numbers.tex` | About 110 `\newcommand` macros, one per number quoted in the text, written by `scripts/run_analysis.py` |
| `refs-supply-chain.bib` | Bibliography for the section |
| `figures/*.pdf` | Figures (PNG copies are for the README) |
| `preview.tex`, `preview.pdf` | The section compiled on its own |

## Putting it into the article on Overleaf

1. Upload `sec-supply-chain.tex`, `refs-supply-chain.bib`, `generated/numbers.tex`, and these six figures from `figures/`: `fig-walks.pdf`, `fig-adjacency.pdf`, `fig-common-neighbours.pdf`, `fig-delta-sweep.pdf`, `fig-spectra.pdf`, `fig-years.pdf`. Keep the folder names. The other figures (2012 vs 2017, residual shares, heat map) are for the repository only.
2. If the old finance section is in the Overleaf project, remove its `\input{sec-finance}`, its `generated/numbers.tex` and `refs-finance.bib`. The two sections define macros with the same names, so they cannot both be loaded.
3. In `main.tex`, the preamble needs `amsmath, amssymb, amsthm, graphicx, booktabs` and the theorem environments `theorem, proposition, corollary, definition, example, remark`. Then add:
   ```latex
   \input{generated/numbers}
   ```
4. Where the section belongs (after the sections on spectra and on Johnson and Hamming graphs):
   ```latex
   \input{sec-supply-chain}
   ```
5. Add the bibliography file alongside the main one:
   ```latex
   \bibliography{refs,refs-supply-chain}
   ```
   Delete from `refs-supply-chain.bib` any entry that is already in `refs.bib`, for example `brouwer2012` and `mentorthesis`.
6. Make sure the main article defines these four labels:

   | Label | Should point to |
   |---|---|
   | `sec:srg` | the section defining strongly regular graphs |
   | `sec:johnson-hamming` | the section on Johnson and Hamming graphs |
   | `sec:spectra` | the section containing "λ₁ lies between the average and the maximum degree" |
   | `thm:walks` | the theorem that (*A*ᵏ)ᵢⱼ counts walks of length *k* |

7. The footnote at the end of the data subsection links to this repository, https://github.com/desaidaksh1206-netizen/supply-chain-graphs-srg, with `\url`. `main.tex` gets `\url` from `hyperref`, and `preview.tex` from the `url` package.

## Before submission

- **Bibliography.** Every entry in `refs-supply-chain.bib` has been checked against the publication of record. Horn & Johnson is cited by section (§5.6 and §8.3) because the second edition's theorem numbers could not be checked. The article previously cited Theorems 5.6.12 and 8.3.1; check those numbers in a print copy of the second edition before citing theorems again.
- **The BEA handbook** (Horowitz & Planting 2009) is a government methodology report, not a book or paper. At the November 2026 Zoom meeting the advisor accepted it, and the BEA tables as the data source.
- **Length.** The section is about eleven pages in the preview, and the advisor agreed to keep that length (November 2026). If it ever needs to be shorter, cut in this order: the spectra figure (keep the sentence about 399 against 3 distinct eigenvalues), then Proposition `prop:colsum` and its paragraph (keep one sentence saying ρ(*A*) was computed directly), then the year-by-year figure (keep the sentence summarising it).

## Changing a number

Never edit a number in `sec-supply-chain.tex`. Change the code or the data, run `python scripts/run_analysis.py`, and upload the new `generated/numbers.tex` and figures. To rebuild only the figures and macros from saved results, run `python scripts/run_analysis.py --reuse`, which takes under a minute.
