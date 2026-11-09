"""Step 3: every computation, figure and number in the paper.

Writes
  results/sweep_2017.csv, results/sweep_2012.csv   one row per threshold (detail level)
  results/years_summary.csv                       one row per year, 2012-2023 (summary level)
  results/walk_shares.csv                         s_k for k = 1..12
  results/linkages.csv                            output multipliers and forward linkages
  results/summary.json                            every scalar quoted in the paper
  paper/generated/numbers.tex                     the same scalars as LaTeX macros
  paper/figures/*.pdf, *.png

  python scripts/run_analysis.py           full run (about 10 minutes)
  python scripts/run_analysis.py --reuse   reuse the saved sweeps; rebuild the rest
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from supplygraphs import analysis as an  # noqa: E402
from supplygraphs import data  # noqa: E402
from supplygraphs import figures as fg  # noqa: E402
from supplygraphs import graphs as gr  # noqa: E402
from supplygraphs import leontief as lf  # noqa: E402

RESULTS = ROOT / "results"
FIGS = ROOT / "paper" / "figures"
GEN = ROOT / "paper" / "generated"
KMAX = 12


def clean(d: dict) -> dict:
    return {k: (None if (isinstance(v, float) and not np.isfinite(v)) else
                (float(v) if isinstance(v, (np.floating, np.integer)) else v))
            for k, v in d.items()}


def main():
    reuse = "--reuse" in sys.argv
    RESULTS.mkdir(exist_ok=True)
    GEN.mkdir(parents=True, exist_ok=True)
    build = json.loads((data.PROCESSED / "build_summary.json").read_text())

    A, ind = data.load(2017)
    A12, _ = data.load(2012)
    n = A.shape[0]
    W, W12 = gr.link_weights(A), gr.link_weights(A12)

    # ------------------------------------------------------------ supply chains as walks
    rho, right, left = lf.perron(A)
    shares = lf.walk_shares(A, KMAX)
    apl = lf.average_chain_length(A)
    L = lf.leontief_inverse(A)
    mult = L.sum(axis=0)
    forward = L.sum(axis=1)
    pd.DataFrame({"k": np.arange(1, KMAX + 1), "share": shares}).to_csv(
        RESULTS / "walk_shares.csv", index=False)
    link = pd.DataFrame({"code": ind["code"], "name": ind["name"], "sector": ind["sector"],
                         "output_multiplier": mult, "forward_linkage": forward,
                         "perron_left": left, "weighted_out_degree": A.sum(axis=1) - np.diag(A),
                         "weighted_in_degree": A.sum(axis=0) - np.diag(A)})
    link.to_csv(RESULTS / "linkages.csv", index=False)
    neumann_err = {k: float(np.abs(lf.neumann_partial_sum(A, k) - L).max()) for k in (5, 10, 20)}
    colsum = A.sum(axis=0)
    rho12 = lf.spectral_radius(A12)
    apl12 = lf.average_chain_length(A12)
    s12 = lf.walk_shares(A12, 3)

    # ------------------------------------------------------------ graphs
    theta_star = an.theta_for_density(W, an.TARGET_DENSITY)
    star = an.at_density(W, an.TARGET_DENSITY, np.random.default_rng(an.SEED + 10))
    a_star = gr.threshold_adjacency(W, theta_star)
    sec = an.group_profile(a_star, ind["sector"])
    sub = an.group_profile(a_star, ind["summary_code"])
    er_struct = gr.structure_profile(gr.erdos_renyi(n, star["density"], np.random.default_rng(an.SEED + 1)))
    deg = a_star.sum(axis=1)
    hubs = [{"name": ind["name"][i], "sector": ind["sector"][i], "degree": int(deg[i])}
            for i in np.argsort(-deg)[:10]]
    isolated = list(ind["name"][deg == 0])
    hub_idx = set(np.argsort(-deg)[:10])
    a2 = gr.square(a_star)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n) if a_star[i, j]]
    hub_edges = [e for e in edges if e[0] in hub_idx or e[1] in hub_idx]
    frac_edges_hub = len(hub_edges) / len(edges)
    star12 = an.at_density(W12, an.TARGET_DENSITY, np.random.default_rng(an.SEED + 12))

    if reuse:
        df = pd.read_csv(RESULTS / "sweep_2017.csv")
        df12 = pd.read_csv(RESULTS / "sweep_2012.csv")
        years = pd.read_csv(RESULTS / "years_summary.csv")
    else:
        df = an.sweep(W)
        df.to_csv(RESULTS / "sweep_2017.csv", index=False)
        df12 = an.sweep(W12, n_er=100, n_dp=10)
        df12.to_csv(RESULTS / "sweep_2012.csv", index=False)
        rows = []
        rng = np.random.default_rng(an.SEED + 100)
        for y in data.SUMMARY_YEARS:
            As = data.load_summary(y)
            r = an.at_density(gr.link_weights(As), an.TARGET_DENSITY, rng)
            r.update({"year": y, "n": As.shape[0], "rho": lf.spectral_radius(As),
                      "apl": lf.average_chain_length(As), "s1": float(lf.walk_shares(As, 1)[0])})
            rows.append(r)
        years = pd.DataFrame(rows)
        years.to_csv(RESULTS / "years_summary.csv", index=False)

    nulls = df.dropna(subset=["er_delta"])
    th = nulls["theta"].to_numpy()
    rng_ = {
        "gap_above_dp_band": an.ranges(th, (nulls["gap"] > nulls["dp_gap_hi"]).to_numpy()),
        "gap_above_er_band": an.ranges(th, (nulls["gap"] > nulls["er_gap_hi"]).to_numpy()),
        "delta_above_er_band": an.ranges(th, (nulls["delta"] > nulls["er_delta_hi"]).to_numpy()),
        "delta_below_dp_band": an.ranges(th, (nulls["delta"] < nulls["dp_delta_lo"]).to_numpy()),
    }
    yrs = {
        "gap_above_dp_all": bool((years["gap"] > years["dp_gap_hi"]).all()),
        "gap_above_er_all": bool((years["gap"] > years["er_gap_hi"]).all()),
        "delta_between_all": bool(((years["delta"] > years["er_delta_hi"])
                                   & (years["delta"] < years["dp_delta_lo"])).all()),
        "n_gap_above_dp": int((years["gap"] > years["dp_gap_hi"]).sum()),
        "n_delta_between": int(((years["delta"] > years["er_delta_hi"])
                                & (years["delta"] < years["dp_delta_lo"])).sum()),
        "rho_min": float(years["rho"].min()), "rho_max": float(years["rho"].max()),
        "rho_min_year": int(years.loc[years["rho"].idxmin(), "year"]),
        "rho_max_year": int(years.loc[years["rho"].idxmax(), "year"]),
        "apl_min": float(years["apl"].min()), "apl_max": float(years["apl"].max()),
        "gap_min": float(years["gap"].min()), "gap_max": float(years["gap"].max()),
        "delta_min": float(years["delta"].min()), "delta_max": float(years["delta"].max()),
        "n_summary": int(years["n"].iloc[0]),
        "years": [int(years["year"].min()), int(years["year"].max())],
    }

    # ------------------------------------------------------------ figures
    fg.fig_adjacency(a_star, ind["sector"], FIGS, theta_star)
    fg.fig_walks(shares, rho, apl, FIGS)
    fg.fig_delta_sweep(df, theta_star, FIGS)
    fg.fig_terms(df, FIGS)
    fg.fig_common_neighbours(a_star, np.random.default_rng(an.SEED + 2), FIGS, theta_star)
    fg.fig_spectra(a_star, np.random.default_rng(an.SEED + 3), FIGS)
    fg.fig_years(years, FIGS)
    fg.fig_detail_years([("2017", df), ("2012", df12)], FIGS)
    fg.fig_requirements_heatmap(A, ind["sector"], FIGS)

    top_mult = link.sort_values("output_multiplier", ascending=False).head(5)
    top_fwd = link.sort_values("forward_linkage", ascending=False).head(5)
    top_perron = link.sort_values("perron_left", ascending=False).head(5)
    summary = {
        "data": build,
        "n": n,
        "leontief": {
            "rho": rho, "rho_2012": rho12, "apl": apl, "apl_2012": apl12,
            "shares": [float(x) for x in shares], "shares_2012": [float(x) for x in s12],
            "share_beyond_3": float(1 - shares[:3].sum()),
            "max_column_sum": float(colsum.max()),
            "max_column_sum_industry": ind["name"][int(np.argmax(colsum))],
            "n_column_sum_above_1": int((colsum > 1).sum()),
            "mean_column_sum": float(colsum.mean()),
            "neumann_max_error": neumann_err,
            "multiplier_min": float(mult.min()), "multiplier_max": float(mult.max()),
            "multiplier_mean": float(mult.mean()),
            "top_multipliers": top_mult[["name", "output_multiplier"]].to_dict(orient="records"),
            "top_forward": top_fwd[["name", "forward_linkage"]].to_dict(orient="records"),
            "top_perron": top_perron[["name", "perron_left"]].to_dict(orient="records"),
        },
        "theta_star": theta_star,
        "at_theta_star": clean(star),
        "at_theta_star_2012": clean(star12),
        "er_sample_at_theta_star": clean(er_struct),
        "sector_at_theta_star": sec,
        "summary_group_at_theta_star": sub,
        "hubs": hubs, "isolated": isolated,
        "frac_edges_touching_top10": frac_edges_hub,
        "ranges": rng_,
        "years": yrs,
        "seed": an.SEED, "n_er": an.N_ER, "n_dp": an.N_DP, "n_dp_star": 50,
        "thetas": [float(df["theta"].min()), float(df["theta"].max()), int(len(df))],
    }
    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    write_macros(summary)
    print(json.dumps({k: summary[k] for k in ("theta_star", "ranges", "years")}, indent=2))


def write_macros(s):
    a, a12 = s["at_theta_star"], s["at_theta_star_2012"]
    lo = s["leontief"]
    sec, sub = s["sector_at_theta_star"], s["summary_group_at_theta_star"]
    er = s["er_sample_at_theta_star"]
    y = s["years"]

    def pct(x, d=0):
        return f"{100 * x:.{d}f}"

    def th(x):
        return f"{x:.4f}" if x < 0.01 else f"{x:.3f}"

    def rg(key):
        r = s["ranges"][key]
        return ", ".join(f"$[{th(lo_)},{th(hi_)}]$" for lo_, hi_ in r) or "no threshold"

    m = {
        "NInd": s["n"], "NIndSource": s["data"]["detail_2017"]["n_industries_source"],
        "NSummary": y["n_summary"], "YearFirst": y["years"][0], "YearLast": y["years"][1],
        "DroppedShare": pct(s["data"]["detail_2017"]["dropped_input_share_max"], 1),
        "Rho": f"{lo['rho']:.3f}", "RhoTwelve": f"{lo['rho_2012']:.3f}",
        "APL": f"{lo['apl']:.2f}", "APLTwelve": f"{lo['apl_2012']:.2f}",
        "ShareOne": pct(lo["shares"][0], 1), "ShareTwo": pct(lo["shares"][1], 1),
        "ShareThree": pct(lo["shares"][2], 1), "ShareBeyondThree": pct(lo["share_beyond_3"], 1),
        "ShareOneTwelve": pct(lo["shares_2012"][0], 1),
        "RatioMax": f"{max(np.array(lo['shares'][1:9]) / np.array(lo['shares'][:8])):.3f}",
        "RatioMin": f"{min(np.array(lo['shares'][1:9]) / np.array(lo['shares'][:8])):.3f}",
        "MaxColSum": f"{lo['max_column_sum']:.2f}", "MeanColSum": f"{lo['mean_column_sum']:.2f}",
        "NColAboveOne": lo["n_column_sum_above_1"],
        "MultMin": f"{lo['multiplier_min']:.2f}", "MultMax": f"{lo['multiplier_max']:.2f}",
        "MultMean": f"{lo['multiplier_mean']:.2f}",
        "NeumannTwenty": f"{lo['neumann_max_error'][20]:.1e}".replace("e-0", r"\times 10^{-").replace("e-", r"\times 10^{-") + "}",
        "ThetaStar": th(s["theta_star"]), "ThetaStarCents": f"{100 * s['theta_star']:.2f}",
        "ThetaStarTwelve": th(a12["theta"]),
        "EdgesStar": int(a["m"]), "DensityStar": f"{a['density']:.3f}",
        "MeanDegStar": f"{a['mean_degree']:.1f}", "SdDegStar": f"{a['sd_degree']:.1f}",
        "MaxDegStar": int(a["max_degree"]), "ComponentsStar": int(a["n_components"]),
        "IsolatedStar": int(a["n_isolated"]), "CliqueStar": int(a["clique_number"]),
        "DistinctStar": int(a["n_distinct_eigenvalues"]),
        "LambdaOneAStar": f"{a['lambda1']:.1f}", "SpecGapStar": f"{a['spectral_gap']:.1f}",
        "ERDistinct": er["n_distinct_eigenvalues"], "ERLambdaOne": f"{er['lambda1']:.1f}",
        "ERSpecGap": f"{er['spectral_gap']:.1f}", "ERClique": int(er["clique_number"]),
        "ERSdDeg": f"{er['sd_degree']:.1f}",
        "ERCommon": f"{(s['n'] - 2) * a['density'] ** 2:.2f}",
        "LambdaFitStar": f"{a['lambda_fit']:.2f}", "MuFitStar": f"{a['mu_fit']:.2f}",
        "GapStar": f"{a['gap']:.2f}", "SdLambdaStar": f"{a['sd_lambda']:.2f}",
        "SdMuStar": f"{a['sd_mu']:.2f}", "DeltaStar": f"{a['delta']:.3f}",
        "FracDegreeStar": pct(a["frac_degree"]), "FracEdgeStar": pct(a["frac_edge"]),
        "FracNonedgeStar": pct(a["frac_nonedge"]),
        "ERDeltaStar": f"{a['er_delta']:.3f}", "ERDeltaLoStar": f"{a['er_delta_lo']:.3f}",
        "ERDeltaHiStar": f"{a['er_delta_hi']:.3f}",
        "ERGapStar": f"{a['er_gap']:.2f}", "ERGapLoStar": f"{a['er_gap_lo']:.2f}",
        "ERGapHiStar": f"{a['er_gap_hi']:.2f}",
        "DPDeltaStar": f"{a['dp_delta']:.3f}", "DPDeltaLoStar": f"{a['dp_delta_lo']:.3f}",
        "DPDeltaHiStar": f"{a['dp_delta_hi']:.3f}",
        "DPGapStar": f"{a['dp_gap']:.2f}", "DPGapLoStar": f"{a['dp_gap_lo']:.2f}",
        "DPGapHiStar": f"{a['dp_gap_hi']:.2f}",
        "GapTwelve": f"{a12['gap']:.2f}", "DeltaTwelve": f"{a12['delta']:.3f}",
        "DPGapHiTwelve": f"{a12['dp_gap_hi']:.2f}", "ERDeltaHiTwelve": f"{a12['er_delta_hi']:.3f}",
        "DPDeltaLoTwelve": f"{a12['dp_delta_lo']:.3f}",
        "RangeGapAboveDP": rg("gap_above_dp_band"), "RangeGapAboveER": rg("gap_above_er_band"),
        "RangeDeltaAboveER": rg("delta_above_er_band"), "RangeDeltaBelowDP": rg("delta_below_dp_band"),
        "FracWithinSector": pct(sec["frac_edges_within"]),
        "FracWithinSectorExp": pct(sec["frac_pairs_within"]),
        "FracWithinSub": pct(sub["frac_edges_within"], 1),
        "FracWithinSubExp": pct(sub["frac_pairs_within"], 1),
        "LambdaWithinSub": f"{sub['lambda_within']:.1f}", "LambdaAcrossSub": f"{sub['lambda_across']:.1f}",
        "MuWithinSub": f"{sub['mu_within']:.1f}", "MuAcrossSub": f"{sub['mu_across']:.1f}",
        "FracHubEdges": pct(s["frac_edges_touching_top10"]),
        "HubOne": s["hubs"][0]["name"].replace("&", r"\&"), "HubOneDeg": s["hubs"][0]["degree"],
        "HubTwo": s["hubs"][1]["name"].replace("&", r"\&"), "HubTwoDeg": s["hubs"][1]["degree"],
        "HubThree": s["hubs"][2]["name"].replace("&", r"\&"), "HubThreeDeg": s["hubs"][2]["degree"],
        "HubFour": s["hubs"][3]["name"].replace("&", r"\&"), "HubFourDeg": s["hubs"][3]["degree"],
        "NER": s["n_er"], "NDP": s["n_dp"], "NDPStar": s["n_dp_star"], "Seed": s["seed"],
        "NThetas": s["thetas"][2],
        "YearsGapAboveDP": y["n_gap_above_dp"], "YearsDeltaBetween": y["n_delta_between"],
        "NYears": y["years"][1] - y["years"][0] + 1,
        "YearsRhoMin": f"{y['rho_min']:.3f}", "YearsRhoMax": f"{y['rho_max']:.3f}",
        "YearsRhoMinYear": y["rho_min_year"], "YearsRhoMaxYear": y["rho_max_year"],
        "YearsAPLMin": f"{y['apl_min']:.2f}", "YearsAPLMax": f"{y['apl_max']:.2f}",
        "YearsGapMin": f"{y['gap_min']:.2f}", "YearsGapMax": f"{y['gap_max']:.2f}",
        "YearsDeltaMin": f"{y['delta_min']:.3f}", "YearsDeltaMax": f"{y['delta_max']:.3f}",
    }
    lines = ["% Generated by scripts/run_analysis.py -- do not edit by hand.",
             f"% Seed {s['seed']}. Every number quoted in sec-supply-chain.tex comes from here."]
    for k, v in m.items():
        lines.append(f"\\newcommand{{\\{k}}}{{{v}}}")
    (GEN / "numbers.tex").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
