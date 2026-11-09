"""Every figure, as PDF (for LaTeX) and PNG (for the README).

Colour roles are fixed across figures: blue = the industry network,
orange = Erdos--Renyi benchmark, aqua = degree-preserving benchmark,
near-black = the strongly regular benchmark H(2,20). In the common-neighbour
figure blue/orange instead mean adjacent / non-adjacent pairs, as its legend says.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

from . import graphs as gr  # noqa: E402

BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"

SHORT = {
    "Agriculture, forestry, fishing, and hunting": "Agriculture",
    "Arts, entertainment, recreation, accommodation, and food services": "Arts, food, lodging",
    "Educational services, health care, and social assistance": "Education & health",
    "Finance, insurance, real estate, rental, and leasing": "Finance & real estate",
    "Professional and business services": "Professional services",
    "Other services, except government": "Other services",
    "Transportation and warehousing": "Transportation",
}

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 9,
    "axes.edgecolor": INK2,
    "axes.labelcolor": INK,
    "axes.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "legend.frameon": False,
    "legend.fontsize": 8,
    "lines.linewidth": 2,
    "savefig.bbox": "tight",
    "savefig.dpi": 200,
    "mathtext.fontset": "cm",
    "pdf.fonttype": 42,
})


def _save(fig, outdir: Path, name: str):
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{name}.pdf")
    fig.savefig(outdir / f"{name}.png")
    plt.close(fig)


def _sector_order(sectors):
    s = np.asarray(sectors)
    order = sorted(range(len(s)), key=lambda i: (s[i], i))
    sorted_s = s[order]
    bounds = [i for i in range(1, len(s)) if sorted_s[i] != sorted_s[i - 1]]
    return np.array(order), bounds


def fig_adjacency(a: np.ndarray, sectors, outdir: Path, theta: float):
    order, bounds = _sector_order(sectors)
    n = len(order)
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 3.3))
    for ax, idx, title in [(axes[0], np.arange(n), "Ordered by BEA industry code"),
                           (axes[1], order, "Grouped by BEA sector")]:
        ax.imshow(a[np.ix_(idx, idx)], cmap="Greys", interpolation="nearest", vmin=0, vmax=1)
        ax.set_title(title, fontsize=9, color=INK)
        ticks = [0, n // 2 - 1, n - 1]
        ax.set_xticks(ticks, [str(t + 1) for t in ticks])
        ax.set_yticks(ticks, [str(t + 1) for t in ticks])
        ax.grid(False)
    for b in bounds:
        axes[1].axhline(b - 0.5, color=ORANGE, lw=0.3, alpha=0.8)
        axes[1].axvline(b - 0.5, color=ORANGE, lw=0.3, alpha=0.8)
    fig.suptitle(rf"Adjacency matrix of $G_\theta$ at $\theta^*={theta:.4f}$ (black = edge)",
                 fontsize=9, color=INK, y=0.98)
    _save(fig, outdir, "fig-adjacency")


def fig_walks(shares: np.ndarray, rho: float, apl: float, outdir: Path):
    k = np.arange(1, len(shares) + 1)
    fig, ax = plt.subplots(figsize=(6.0, 2.9))
    ax.bar(k, shares, color=BLUE, width=0.7, label=r"$s_k$: share of requirements via chains of length $k$")
    ax.plot(k, shares[0] * rho ** (k - 1), color=INK, lw=1.2, ls="--", marker="o", ms=3,
            label=rf"$s_1\rho^{{k-1}}$, $\rho(A)={rho:.3f}$")
    for x, y in zip(k[:4], shares[:4]):
        ax.text(x + 0.4, y + 0.004, f"{100 * y:.1f}%", ha="left", fontsize=7, color=INK2)
    ax.set_xticks(k)
    ax.set_xlabel("supply-chain length $k$ (steps)")
    ax.set_ylabel("share of total requirements")
    ax.set_title(f"Average supply-chain length {apl:.2f} steps", fontsize=9, loc="left", color=INK)
    ax.legend(loc="upper right")
    _save(fig, outdir, "fig-walks")


def _band(ax, x, df, prefix, color, label):
    ax.fill_between(x, df[f"{prefix}_lo"], df[f"{prefix}_hi"], color=color, alpha=0.22, lw=0)
    ax.plot(x, df[prefix], color=color, lw=1.2, label=label)


def fig_delta_sweep(df: pd.DataFrame, theta_star: float, outdir: Path):
    d = df.dropna(subset=["er_delta"])
    x = d["theta"]
    fig, axes = plt.subplots(2, 1, figsize=(6.0, 5.2), sharex=True)
    ax = axes[0]
    _band(ax, x, d, "er_delta", ORANGE, "Erdős–Rényi, same density")
    _band(ax, x, d, "dp_delta", AQUA, "Same degrees, rewired")
    ax.plot(x, d["delta"], color=BLUE, label=r"Industry network $G_\theta$")
    ax.axhline(0, color=INK, lw=1.5)
    ax.text(x.max(), 0.025, r"$\delta=0$: every strongly regular graph", color=INK,
            fontsize=8, ha="right")
    ax.set_ylabel(r"$\delta(G)$")
    ax.set_ylim(-0.03, 1.0)
    ax.legend(loc="upper left", fontsize=7.5)
    ax.set_title("(a) Distance from strong regularity", fontsize=9, loc="left", color=INK)
    ax = axes[1]
    _band(ax, x, d, "er_gap", ORANGE, None)
    _band(ax, x, d, "dp_gap", AQUA, None)
    ax.plot(x, d["gap"], color=BLUE)
    ax.axhline(0, color=GRID, lw=0.8)
    ax.plot([theta_star], [16], marker="s", ms=6, color=INK, ls="none")
    ax.annotate(r"$H(2,20)$: $\lambda-\mu=16$", xy=(theta_star, 16),
                xytext=(theta_star * 2.2, 24), fontsize=8, color=INK,
                arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
    ax.set_ylabel(r"$\bar\lambda-\bar\mu$")
    ax.set_xlabel(r"threshold $\theta$ (dollars of input per dollar of output, log scale)")
    ax.set_title("(b) Gap between mean common neighbours on edges and on non-edges",
                 fontsize=9, loc="left", color=INK)
    for a_ in axes:
        a_.set_xscale("log")
        a_.axvline(theta_star, color=INK2, lw=0.8, ls=":")
    axes[0].text(theta_star * 1.08, 0.12, r"$\theta^*$", color=INK2, fontsize=9)
    _save(fig, outdir, "fig-delta-sweep")


def fig_terms(df: pd.DataFrame, outdir: Path):
    d = df.dropna(subset=["frac_degree"])
    fig, ax = plt.subplots(figsize=(6.0, 2.8))
    ax.stackplot(d["theta"], d["frac_degree"], d["frac_edge"], d["frac_nonedge"],
                 colors=[BLUE, ORANGE, AQUA], edgecolor="white", linewidth=0.8,
                 labels=[r"degrees $n\,\mathrm{Var}(\mathcal{D})$",
                         r"edges $2m\,\mathrm{Var}(\mathcal{L})$",
                         r"non-edges $2\bar m\,\mathrm{Var}(\mathcal{M})$"])
    ax.set_xscale("log")
    ax.set_ylim(0, 1)
    ax.set_xlim(d["theta"].min(), d["theta"].max())
    ax.set_ylabel("share of squared residual")
    ax.set_xlabel(r"threshold $\theta$")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3)
    _save(fig, outdir, "fig-residual-terms")


def fig_common_neighbours(a: np.ndarray, rng: np.random.Generator, outdir: Path, theta: float):
    from scipy.stats import binom

    n = a.shape[0]
    p = a.sum() / (n * (n - 1))
    er = gr.erdos_renyi(n, p, rng)
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.9), sharey=True)
    for ax, g, title in [(axes[0], a, rf"Industry network $G_{{\theta^*}}$"),
                         (axes[1], er, rf"Erdős–Rényi sample, $p={p:.3f}$")]:
        _, lam, mu = gr.common_neighbour_lists(g)
        top = int(np.quantile(np.concatenate([lam, mu]), 0.995)) + 2
        bins = np.arange(-0.5, top + 0.5, 1)
        ax.hist(mu, bins=bins, density=True, histtype="stepfilled", color=ORANGE, alpha=0.35, lw=0)
        ax.hist(mu, bins=bins, density=True, histtype="step", color=ORANGE, lw=1.2,
                label=r"non-adjacent pairs $\mathcal{M}$")
        ax.hist(lam, bins=bins, density=True, histtype="step", color=BLUE, lw=1.5,
                label=r"adjacent pairs $\mathcal{L}$")
        k = np.arange(0, top)
        ax.plot(k, binom.pmf(k, n - 2, p * p), color=MUTED, lw=1.0, ls="--",
                label=r"$\mathrm{Bin}(n-2,p^2)$")
        ax.set_title(title, fontsize=9, color=INK)
        ax.set_xlabel("common neighbours $c(u,v)$")
        ax.set_xlim(-0.5, top)
    ymax = axes[0].get_ylim()[1] * 1.25
    for x, lbl, y in [(18, r"$\lambda=18$", 0.62), (2, r"$\mu=2$", 0.80)]:
        axes[0].axvline(x, color=INK, lw=1.2)
        axes[0].text(x + 1.0, ymax * y, lbl, fontsize=8, color=INK,
                     bbox=dict(facecolor="white", edgecolor="none", pad=0.5))
    axes[0].text(axes[0].get_xlim()[1] * 0.98, ymax * 0.93, r"black lines: $H(2,20)$",
                 fontsize=8, color=INK, ha="right")
    axes[0].set_ylabel("proportion of pairs")
    axes[0].set_ylim(0, ymax * 1.25)
    axes[1].legend(loc="upper right", fontsize=7)
    _save(fig, outdir, "fig-common-neighbours")


def fig_spectra(a: np.ndarray, rng: np.random.Generator, outdir: Path):
    n = a.shape[0]
    p = a.sum() / (n * (n - 1))
    graphs = [(gr.hamming_graph(2, 20), r"$H(2,20)$: strongly regular, three distinct eigenvalues"),
              (gr.erdos_renyi(n, p, rng), "Erdős–Rényi graph of the same density"),
              (a, r"Industry network $G_{\theta^*}$")]
    lo = min(gr.spectrum(g).min() for g, _ in graphs) - 1
    hi = max(gr.spectrum(g).max() for g, _ in graphs) + 1
    fig, axes = plt.subplots(3, 1, figsize=(6.0, 5.0), sharex=True, gridspec_kw={"hspace": 0.45})
    bins = np.arange(np.floor(lo) - 0.5, np.ceil(hi) + 1, 1.0)
    for ax, (g, title) in zip(axes, graphs):
        ax.hist(gr.spectrum(g), bins=bins, color=BLUE)
        ax.set_title(title, fontsize=9, loc="left", color=INK)
        ax.set_ylabel("count")
    for x, mult in [(38, 1), (18, 38), (-2, 361)]:
        axes[0].text(x, min(mult, 330) + 15, f"×{mult}", ha="center", fontsize=7, color=INK2)
    axes[0].set_ylim(0, 400)
    axes[-1].set_xlabel("eigenvalue of the adjacency matrix")
    _save(fig, outdir, "fig-spectra")


def fig_years(years: pd.DataFrame, outdir: Path):
    x = years["year"]
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.6))
    ax = axes[0]
    ax.plot(x, years["rho"], color=BLUE, marker="o", ms=3)
    ax.set_title(r"(a) $\rho(A)$", fontsize=9, loc="left", color=INK)
    ax = axes[1]
    _band(ax, x, years, "er_gap", ORANGE, "Erdős–Rényi")
    _band(ax, x, years, "dp_gap", AQUA, "same degrees")
    ax.plot(x, years["gap"], color=BLUE, marker="o", ms=3, label="industry network (71 industries)")
    ax.set_title(r"(b) $\bar\lambda-\bar\mu$", fontsize=9, loc="left", color=INK)
    ax = axes[2]
    _band(ax, x, years, "er_delta", ORANGE, None)
    _band(ax, x, years, "dp_delta", AQUA, None)
    ax.plot(x, years["delta"], color=BLUE, marker="o", ms=3)
    ax.set_title(r"(c) $\delta$", fontsize=9, loc="left", color=INK)
    for a_ in axes:
        a_.set_xticks([2012, 2015, 2018, 2021], ["2012", "2015", "2018", "2021"])
        a_.axvspan(2019.6, 2020.4, color=GRID, alpha=0.6, lw=0)
    handles, labels = axes[1].get_legend_handles_labels()
    order = [2, 0, 1]
    fig.legend([handles[i] for i in order], [labels[i] for i in order], loc="lower center",
               ncol=3, bbox_to_anchor=(0.5, -0.06), fontsize=8)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    _save(fig, outdir, "fig-years")


def fig_detail_years(curves: list[tuple[str, pd.DataFrame]], outdir: Path):
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.8), sharex=True)
    styles = [dict(color=BLUE, ls="-"), dict(color=ORANGE, ls="--")]
    for (label, df), st in zip(curves, styles):
        axes[0].plot(df["theta"], df["delta"], label=label, lw=1.5, **st)
        axes[1].plot(df["theta"], df["gap"], lw=1.5, **st)
    axes[0].set_ylabel(r"$\delta(G_\theta)$")
    axes[1].set_ylabel(r"$\bar\lambda-\bar\mu$")
    for ax in axes:
        ax.set_xscale("log")
        ax.set_xlabel(r"threshold $\theta$")
    axes[0].legend(loc="upper left")
    _save(fig, outdir, "fig-detail-years")


def fig_requirements_heatmap(A: np.ndarray, sectors, outdir: Path):
    order, bounds = _sector_order(sectors)
    s = np.asarray(sectors)[order]
    M = A[np.ix_(order, order)].copy()
    M[M <= 0] = np.nan
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    im = ax.imshow(M, cmap="Blues", norm=LogNorm(vmin=1e-5, vmax=0.3), interpolation="nearest")
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    starts = [0] + bounds + [len(s)]
    for a_, b_ in zip(starts[:-1], starts[1:]):
        if b_ - a_ >= 9:
            ax.text(-4, (a_ + b_) / 2, SHORT.get(s[a_], s[a_]), ha="right", va="center",
                fontsize=6, color=INK2)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04,
                 label=r"$a_{ij}$ (\$ of row industry per \$ of column industry)")
    ax.set_title("Direct requirements $A$, 2017, grouped by sector", fontsize=9, color=INK)
    _save(fig, outdir, "fig-requirements-heatmap")
