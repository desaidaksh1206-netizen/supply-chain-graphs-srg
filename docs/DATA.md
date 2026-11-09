# Data: sources, construction and cleaning

Only public data is used. Nothing here comes from any employer, client or proprietary source. Every file can be traced to one pinned public commit, and every change made to the data is listed below.

## 1. Source

The **U.S. Bureau of Economic Analysis (BEA)** publishes the input-output accounts of the US economy: for each year, a **Make** table (which industry produces which commodity) and a **Use** table (which industry buys which commodity), in millions of dollars. We take these tables as redistributed, unchanged, by the U.S. Environmental Protection Agency's open-source R package **`useeior`** (Li et al. 2022, *Applied Sciences* 12(9), 4469). That package stores each BEA table as an R data file.

| | |
|---|---|
| Repository | `github.com/USEPA/useeior` |
| Pinned commit | `1565a70b7da2e50ff7ab5fc344fc0ba95f307a8a` (v1.7.2, 28 July 2025) |
| Files | 34 `.rda` files, listed with their SHA-256 in `data/SHA256SUMS.json` |
| Tables | Make and Use, **before redefinitions**, Use at **producer prices**, **2017 BEA industry classification** |
| Detail level | 402 industries; years 2012 and 2017 |
| Summary level | 71 industries; every year 2012–2023 |
| Also used | `MasterCrosswalk2017` (detail → summary → sector codes), industry names |

`scripts/fetch_data.py` downloads each file from `raw.githubusercontent.com` at the pinned commit and stops if any checksum differs from `data/SHA256SUMS.json`. Reading `.rda` files needs the `pyreadr` package.

**Why these choices.** *Before redefinitions* keeps each establishment's secondary products where BEA recorded them, which is the form the Make–Use construction below expects. *Producer prices* exclude trade and transport margins from each purchase, which BEA records separately as purchases from the trade and transportation industries. *The 2017 classification* is available for both 2012 and 2017 at the detail level, so the two years have identical industry codes and can be compared vertex by vertex.

## 2. Building the requirements matrix

*Code:* `data.direct_requirements`

From Use *U* (commodity × industry, intermediate purchases only) and Make *V* (industry × commodity), with *g* = *V*1 and *q* = *V*ᵀ1:

  *A* = *D B*, with *B* = *U* diag(*g*)⁻¹ and *D* = *V* diag(*q*)⁻¹.

This is the industry-technology construction of Miller & Blair (2009, Ch. 5) and of BEA's handbook (Horowitz & Planting 2009). In the detail tables, two commodities are produced by no industry: *noncomparable imports* and *used and secondhand goods*. For these *q* = 0, so their columns of *D* are set to zero, and purchases of them count as inputs from outside the network, like imports and value added. (At the summary level every commodity has positive output.) Total rows and columns ("Total Intermediate", "Total Industry Output", and so on) are excluded.

BEA's published industry-by-industry tables also make a small adjustment for scrap, which we omit.

## 3. Selecting 400 industries

The detail tables have 402 industries. We keep the **400 with the largest 2017 gross output**, so that the network has exactly the order of the strongly regular benchmark *H*(2,20). The two dropped industries are the smallest:

| Code | Industry | 2017 gross output |
|---|---|---|
| 333316 | Photographic and photocopying equipment manufacturing | $1.6 bn |
| 335110 | Electric lamp bulb and part manufacturing | $1.5 bn |

Restricting *A* to the 400 treats these two as outside suppliers. No remaining industry gets more than **3.6%** of its inputs from them (`build_summary.json`, `dropped_input_share_max`).

The 400 are numbered in **alphabetical order of BEA code**. Each has a sector (15 BEA sectors) and a summary code (71 summary industries) from the crosswalk. Manufacturing accounts for 230 of the 400, because the detail level splits manufacturing most finely.

Two of the 400 are accounting constructs rather than producers: **"Customs duties"** (4200ID) and **"Private households"** (814000). Neither buys intermediate inputs, so both have output multiplier exactly 1 and are isolated in every threshold graph in the sweep. They are kept, because dropping them would break the match with *H*(2,20), and the text names them.

## 4. Cleaning: negative entries set to zero

A handful of Use entries are negative. In the detail tables these are federal general government (nondefense) purchases of farm commodities, such as oilseeds, grains and dairy. They record sales of farm products *from* government inventories (Commodity Credit Corporation), which BEA enters as negative purchases.

| Table | Negative entries | Total |
|---|---|---|
| Detail 2012 | 8 | −$1,161 million |
| Detail 2017 | 7 | −$697 million |
| Summary 2012–2023 | 5–7 per year | — |

They are set to **zero**, so that *A* is entrywise non-negative. Non-negativity is what the Perron–Frobenius theorem and the walk interpretation need. For scale, total gross output in 2017 was $34.5 trillion. Before the change, the most negative entry of *A* was −0.0006. Setting them to zero changes no reported number at the precision quoted (checked by running the analysis both ways).

## 5. A data feature, not an error: one column sum above 1

Column *j* of *A* sums to the share of industry *j*'s output spent on intermediate inputs, normally below 1. **State and local government passenger transit** (S00201) has a column sum of **1.59** in 2017: its intermediate purchases exceed the value of its output, so its value added in the table is negative (−$10.1 bn). This is how BEA records a heavily subsidised government enterprise, not a data error, and it is kept. Its consequence for the mathematics is discussed in [`METHODS.md`](METHODS.md) §2: the column-sum bound on ρ(*A*) cannot be used, so ρ(*A*) is computed directly.

## 6. Processed files (committed)

| File | Contents |
|---|---|
| `industries.csv` | index (1–400), code, name, sector code, sector, summary code, gross output 2017 and 2012 ($ million) |
| `dropped_industries.csv` | the two industries not used |
| `A_detail_2017.csv`, `A_detail_2012.csv` | 400 × 400 direct requirements, rows and columns in the vertex ordering |
| `A_summary_2012.csv` … `A_summary_2023.csv` | 71 × 71 direct requirements |
| `summary_industries.csv` | the 71 summary codes and their sectors |
| `build_summary.json` | counts, column-sum checks, negative entries, the pinned commit |

`tests/test_data.py` checks shapes, ordering, non-negativity and ρ(*A*) < 1 for every table.

## 7. Licences

- **BEA tables:** works of the US federal government, not subject to copyright in the US.
- **`useeior`:** MIT licence, © 2021 U.S. EPA. Only data files from it are used, not code.
- **This repository's code:** MIT (see `LICENSE`).

## 8. Citing the data in the article

`refs-supply-chain.bib` cites the BEA tables (`bea2025io`, with the commit in `howpublished`) and the peer-reviewed paper describing `useeior` (`li2022useeior`). Under the advisor's rule, books and papers support the *mathematics*. The BEA tables are cited as the source of the *numbers*, which is standard practice.
