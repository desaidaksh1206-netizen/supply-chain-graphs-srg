"""Download, verify and assemble the US input-output data.

Source
------
The U.S. Bureau of Economic Analysis (BEA) publishes Make and Use tables of the
US economy. We take them as redistributed, unchanged, by the U.S. Environmental
Protection Agency's open-source package ``useeior`` (Li et al. 2022), which
stores each table as an R data file. Every file is fetched from one pinned git
commit and checked against a SHA-256 digest.

Tables used (all "before redefinitions", producer prices, 2017 industry codes):

* Detail level (402 industries): 2017 and 2012.
* Summary level (71 industries): every year 2012--2023.

Construction (Miller & Blair 2009, Ch. 5; Horowitz & Planting 2009)
-------------------------------------------------------------------
U  commodity-by-industry intermediate use      ($ of commodity c bought by industry j)
V  industry-by-commodity make                  ($ of commodity c made by industry i)
g  industry output = row sums of V;  q  commodity output = column sums of V
B = U diag(g)^-1   commodity inputs per $ of industry output
D = V diag(q)^-1   market shares: fraction of commodity c made by industry i
A = D B            industry-by-industry direct requirements:
                   a_ij = $ of industry i's output used per $ of industry j's output

Commodities that no industry produces (in the detail tables: noncomparable
imports, and used and secondhand goods) have q = 0; their column of D is set to
zero, so they count as inputs from outside the network, like imports and value
added.
"""

from __future__ import annotations

import hashlib
import json
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"

USEEIOR_COMMIT = "1565a70b7da2e50ff7ab5fc344fc0ba95f307a8a"  # useeior v1.7.2, 2025-07-28
BASE = f"https://raw.githubusercontent.com/USEPA/useeior/{USEEIOR_COMMIT}/data/"

DETAIL_YEARS = (2012, 2017)
SUMMARY_YEARS = tuple(range(2012, 2024))
MAIN_YEAR = 2017
N_INDUSTRIES = 400  # order of the Hamming graph H(2,20)
NEGATIVE_USES: list[tuple[int, float]] = []  # (count, $ million) per table built


def _files() -> list[str]:
    f = ["MasterCrosswalk2017", "Detail_IndustryCodeName_2017", "Detail_CommodityCodeName_2017",
         "Detail_GrossOutput_IO_17sch", "Summary_IndustryCodeName_2017",
         "Sector_IndustryCodeName_2017"]
    for y in DETAIL_YEARS:
        f += [f"Detail_Use_{y}_PRO_BeforeRedef_17sch", f"Detail_Make_{y}_BeforeRedef_17sch"]
    for y in SUMMARY_YEARS:
        f += [f"Summary_Use_{y}_PRO_BeforeRedef_17sch", f"Summary_Make_{y}_BeforeRedef_17sch"]
    return [x + ".rda" for x in f]


PINNED = ROOT / "data" / "SHA256SUMS.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(force: bool = False) -> dict[str, str]:
    """Download every table into data/raw and verify it against data/SHA256SUMS.json.

    If the checksum file does not exist yet (first run by the author) it is created.
    """
    RAW.mkdir(parents=True, exist_ok=True)
    pinned = json.loads(PINNED.read_text()) if PINNED.exists() else {}
    digests = {}
    for name in _files():
        path = RAW / name
        if force or not path.exists():
            print(f"downloading {name}")
            with urllib.request.urlopen(BASE + name, timeout=120) as r:
                path.write_bytes(r.read())
        d = sha256(path)
        if name in pinned and pinned[name] != d:
            raise RuntimeError(f"{name}: checksum {d} does not match pinned {pinned[name]}")
        digests[name] = d
    if not pinned:
        PINNED.write_text(json.dumps(digests, indent=2) + "\n")
    return digests


def _read(name: str) -> pd.DataFrame:
    import pyreadr  # only needed when building from the raw R files

    r = pyreadr.read_r(str(RAW / f"{name}.rda"))
    return r[list(r)[0]]


def direct_requirements(use: pd.DataFrame, make: pd.DataFrame,
                        industries: list[str], commodities: list[str]) -> tuple[pd.DataFrame, pd.Series]:
    """A = D B (industry-by-industry) and industry output g, from BEA tables."""
    U = use.loc[commodities, industries].astype(float).fillna(0.0).to_numpy()
    # A handful of entries are negative: federal government sales of farm
    # products from inventory (Commodity Credit Corporation), which BEA records
    # as negative purchases. They are set to zero so that A is non-negative.
    NEGATIVE_USES.append((int((U < 0).sum()), float(-U[U < 0].sum())))
    U = np.clip(U, 0.0, None)
    V = make.loc[industries, commodities].astype(float).fillna(0.0).to_numpy()
    g = V.sum(axis=1)
    q = V.sum(axis=0)
    if (g <= 0).any():
        raise ValueError("an industry has no output")
    B = U / g[None, :]
    D = np.divide(V, q[None, :], out=np.zeros_like(V), where=q[None, :] > 0)
    A = D @ B
    return (pd.DataFrame(A, index=industries, columns=industries),
            pd.Series(g, index=industries))


def _codes(level: str) -> tuple[list[str], list[str]]:
    if level == "Detail":
        ind = _read("Detail_IndustryCodeName_2017").iloc[:, 0].tolist()
        com = _read("Detail_CommodityCodeName_2017").iloc[:, 0].tolist()
    else:
        ind = _read("Summary_IndustryCodeName_2017").iloc[:, 0].tolist()
        com = None
    return ind, com


def build(n_industries: int = N_INDUSTRIES) -> dict:
    """Assemble data/processed from data/raw. Returns a summary dict."""
    PROCESSED.mkdir(parents=True, exist_ok=True)
    cw = _read("MasterCrosswalk2017").drop_duplicates("BEA_2017_Detail_Code")
    cw = cw.set_index("BEA_2017_Detail_Code")
    sector_names = _read("Sector_IndustryCodeName_2017").set_index(
        "BEA_2017_Sector_Industry_Code").iloc[:, 0]
    names = _read("Detail_IndustryCodeName_2017").set_index(
        "BEA_2017_Detail_Industry_Code").iloc[:, 0]

    ind, com = _codes("Detail")
    out = {}
    detail = {}
    NEGATIVE_USES.clear()
    for y in DETAIL_YEARS:
        A, g = direct_requirements(_read(f"Detail_Use_{y}_PRO_BeforeRedef_17sch"),
                                   _read(f"Detail_Make_{y}_BeforeRedef_17sch"), ind, com)
        detail[y] = (A, g)
        out[f"negative_uses_detail_{y}"] = {"count": NEGATIVE_USES[-1][0],
                                            "total_musd": NEGATIVE_USES[-1][1]}

    # The 400 largest industries by 2017 gross output, in alphabetical order of code.
    A17, g17 = detail[MAIN_YEAR]
    keep_by_size = g17.sort_values(ascending=False).index[:n_industries]
    dropped = sorted(set(ind) - set(keep_by_size))
    chosen = sorted(keep_by_size)

    industries = pd.DataFrame({
        "index": range(1, n_industries + 1),
        "code": chosen,
        "name": [names[c] for c in chosen],
        "sector_code": [cw.loc[c, "BEA_2017_Sector_Code"] for c in chosen],
        "sector": [sector_names[cw.loc[c, "BEA_2017_Sector_Code"]] for c in chosen],
        "summary_code": [cw.loc[c, "BEA_2017_Summary_Code"] for c in chosen],
        "gross_output_2017_musd": [g17[c] for c in chosen],
        "gross_output_2012_musd": [detail[2012][1][c] for c in chosen],
    })
    industries.to_csv(PROCESSED / "industries.csv", index=False)
    pd.DataFrame({"code": dropped, "name": [names[c] for c in dropped],
                  "gross_output_2017_musd": [g17[c] for c in dropped]}).to_csv(
        PROCESSED / "dropped_industries.csv", index=False)

    for y, (A, g) in detail.items():
        # Restricting to the chosen industries treats the dropped ones as outside
        # suppliers, exactly like imports; their share of inputs is reported below.
        A.loc[chosen, chosen].to_csv(PROCESSED / f"A_detail_{y}.csv", float_format="%.10g")
        colsum_full = A.sum(axis=0)
        out[f"detail_{y}"] = out.get(f"detail_{y}", {})
        out[f"detail_{y}"] = {
            "n_industries_source": len(ind),
            "max_column_sum": float(colsum_full.max()),
            "mean_column_sum": float(colsum_full.mean()),
            "dropped_input_share_max": float(A.loc[dropped, chosen].sum(axis=0).max()),
            "total_output_busd": float(g.sum() / 1000),
        }

    summary_ind, _ = _codes("Summary")
    sector_of_summary = cw.drop_duplicates("BEA_2017_Summary_Code").set_index(
        "BEA_2017_Summary_Code")["BEA_2017_Sector_Code"]
    for y in SUMMARY_YEARS:
        use = _read(f"Summary_Use_{y}_PRO_BeforeRedef_17sch")
        make = _read(f"Summary_Make_{y}_BeforeRedef_17sch")
        s_ind = [c for c in summary_ind if c in use.columns and c in make.index]
        s_com = [c for c in make.columns if c in use.index and not str(c).startswith("Total")]
        A, g = direct_requirements(use, make, s_ind, s_com)
        A.to_csv(PROCESSED / f"A_summary_{y}.csv", float_format="%.10g")
        out[f"summary_{y}"] = {"n_industries": len(s_ind),
                               "max_column_sum": float(A.sum(axis=0).max()),
                               "negative_uses_count": NEGATIVE_USES[-1][0]}
    pd.DataFrame({"code": s_ind,
                  "sector_code": [sector_of_summary.get(c, "G") for c in s_ind]}).to_csv(
        PROCESSED / "summary_industries.csv", index=False)

    out.update({"n_industries": n_industries, "dropped": dropped,
                "detail_years": list(DETAIL_YEARS), "summary_years": list(SUMMARY_YEARS),
                "useeior_commit": USEEIOR_COMMIT})
    (PROCESSED / "build_summary.json").write_text(json.dumps(out, indent=2) + "\n")
    return out


def load(year: int = MAIN_YEAR) -> tuple[np.ndarray, pd.DataFrame]:
    """Return (A as an N x N array in the vertex ordering, industries table)."""
    industries = pd.read_csv(PROCESSED / "industries.csv")
    A = pd.read_csv(PROCESSED / f"A_detail_{year}.csv", index_col=0)
    assert list(A.index) == list(industries["code"]) == list(A.columns), "ordering mismatch"
    return A.to_numpy(), industries


def load_summary(year: int) -> np.ndarray:
    return pd.read_csv(PROCESSED / f"A_summary_{year}.csv", index_col=0).to_numpy()
