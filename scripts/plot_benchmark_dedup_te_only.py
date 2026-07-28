#!/usr/bin/env python3
"""
Generate evaluation plots for deduplication benchmark results.

Reads TSV/CSV tables from benchmark_dedup_te_only/ and writes PNG figures
to benchmark_dedup_te_only/plots/, then appends a figures section to REPORT.md.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib_venn import venn3

# Colorblind-friendly palette
COLORS = {
    "none": "#4C78A8",
    "fastq": "#F58518",
    "bam": "#54A24B",
}
METHOD_ORDER = ["none", "fastq", "bam"]
METHOD_LABELS = {
    "none": "none (baseline)",
    "fastq": "fastq (clumpify)",
    "bam": "bam (multimapper)",
}


def read_table(path: Path) -> pd.DataFrame:
    """Auto-detect comma vs tab separator; strip BOM/CR."""
    text = path.read_text(encoding="utf-8-sig")
    sep = "," if text.splitlines()[0].count(",") > text.splitlines()[0].count("\t") else "\t"
    df = pd.read_csv(path, sep=sep)
    df.columns = [c.strip().lstrip("\ufeff") for c in df.columns]
    if "sample" in df.columns:
        df["sample"] = df["sample"].astype(str).str.strip()
    return df


def style_ax(ax, title: str, ylabel: str = "", xlabel: str = ""):
    ax.set_title(title, fontsize=12, pad=10)
    if ylabel:
        ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle=":", alpha=0.4)


def save(fig, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {path}")


def plot_fragments_grouped(reads: pd.DataFrame, out: Path):
    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(reads))
    width = 0.25
    for i, m in enumerate(METHOD_ORDER):
        ax.bar(
            x + (i - 1) * width,
            reads[f"{m}_fragments"] / 1e6,
            width,
            label=METHOD_LABELS[m],
            color=COLORS[m],
        )
    ax.set_xticks(x)
    ax.set_xticklabels(reads["sample"], rotation=35, ha="right")
    style_ax(ax, "Final usable fragments by sample (Allo BAM read1)", "Fragments (millions)")
    ax.legend(frameon=False)
    save(fig, out / "01_fragments_by_sample.png")


def plot_pct_removed(reads: pd.DataFrame, out: Path):
    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(reads))
    width = 0.35
    ax.bar(x - width / 2, reads["fastq_pct_removed_vs_baseline"], width, label=METHOD_LABELS["fastq"], color=COLORS["fastq"])
    ax.bar(x + width / 2, reads["bam_pct_removed_vs_baseline"], width, label=METHOD_LABELS["bam"], color=COLORS["bam"])
    ax.set_xticks(x)
    ax.set_xticklabels(reads["sample"], rotation=35, ha="right")
    style_ax(ax, "Percent fragments removed vs none baseline", "% removed")
    ax.legend(frameon=False)
    save(fig, out / "02_pct_reads_removed.png")


def plot_mean_removed_summary(reads: pd.DataFrame, out: Path):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    means_frag = [reads[f"{m}_fragments"].mean() / 1e6 for m in METHOD_ORDER]
    means_rm = [
        0.0,
        reads["fastq_pct_removed_vs_baseline"].mean(),
        reads["bam_pct_removed_vs_baseline"].mean(),
    ]
    axes[0].bar(METHOD_ORDER, means_frag, color=[COLORS[m] for m in METHOD_ORDER])
    style_ax(axes[0], "Mean fragments", "Fragments (millions)")
    axes[0].set_xticks(range(len(METHOD_ORDER)))
    axes[0].set_xticklabels([METHOD_LABELS[m] for m in METHOD_ORDER], rotation=20, ha="right")

    axes[1].bar(["fastq", "bam"], means_rm[1:], color=[COLORS["fastq"], COLORS["bam"]])
    style_ax(axes[1], "Mean % removed vs none", "% removed")
    save(fig, out / "03_mean_fragments_and_removal.png")


def plot_peak_counts(peaks: pd.DataFrame, out: Path):
    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(peaks))
    width = 0.25
    for i, m in enumerate(METHOD_ORDER):
        ax.bar(
            x + (i - 1) * width,
            peaks[f"{m}_n_peaks"] / 1000,
            width,
            label=METHOD_LABELS[m],
            color=COLORS[m],
        )
    ax.set_xticks(x)
    ax.set_xticklabels(peaks["sample"], rotation=35, ha="right")
    style_ax(ax, "MACS3 peak counts by sample", "Peaks (thousands)")
    ax.legend(frameon=False)
    save(fig, out / "04_peak_counts_by_sample.png")


def plot_peak_delta(peaks: pd.DataFrame, out: Path):
    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(peaks))
    width = 0.35
    ax.bar(x - width / 2, peaks["fastq_delta_peaks_vs_baseline"], width, label="fastq − none", color=COLORS["fastq"])
    ax.bar(x + width / 2, peaks["bam_delta_peaks_vs_baseline"], width, label="bam − none", color=COLORS["bam"])
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(peaks["sample"], rotation=35, ha="right")
    style_ax(ax, "Peak count change vs none baseline", "Δ peaks")
    ax.legend(frameon=False)
    save(fig, out / "05_peak_delta_vs_baseline.png")


def plot_jaccard_bars(ov: pd.DataFrame, out: Path):
    # Support both column naming schemes
    pairs = []
    for a, b, key in [
        ("bam", "fastq", "bam_fastq"),
        ("bam", "none", "bam_none"),
        ("fastq", "none", "fastq_none"),
    ]:
        col = f"jaccard_peakset_{a}_{b}" if f"jaccard_peakset_{a}_{b}" in ov.columns else f"jaccard_{a}_{b}"
        pairs.append((f"{a}↔{b}", col))

    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(ov))
    width = 0.25
    for i, (label, col) in enumerate(pairs):
        ax.bar(x + (i - 1) * width, ov[col], width, label=label)
    ax.set_xticks(x)
    ax.set_xticklabels(ov["sample"], rotation=35, ha="right")
    ax.set_ylim(0, 1.05)
    style_ax(ax, "Peak-set Jaccard index by sample", "Jaccard")
    ax.legend(frameon=False)
    save(fig, out / "06_jaccard_by_sample.png")


def plot_jaccard_heatmap(ov: pd.DataFrame, out: Path):
    def jcol(a, b):
        c1 = f"jaccard_peakset_{a}_{b}"
        c2 = f"jaccard_{a}_{b}"
        return c1 if c1 in ov.columns else c2

    labels = METHOD_ORDER
    mat = np.eye(3)
    means = {
        ("bam", "fastq"): ov[jcol("bam", "fastq")].mean(),
        ("bam", "none"): ov[jcol("bam", "none")].mean(),
        ("fastq", "none"): ov[jcol("fastq", "none")].mean(),
    }
    idx = {m: i for i, m in enumerate(labels)}
    for (a, b), v in means.items():
        mat[idx[a], idx[b]] = v
        mat[idx[b], idx[a]] = v

    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    im = ax.imshow(mat, cmap="Blues", vmin=0.6, vmax=1.0)
    ax.set_xticks(range(3))
    ax.set_yticks(range(3))
    ax.set_xticklabels([METHOD_LABELS[m] for m in labels], rotation=25, ha="right")
    ax.set_yticklabels([METHOD_LABELS[m] for m in labels])
    for i in range(3):
        for j in range(3):
            ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", color="black")
    ax.set_title("Mean peak-set Jaccard across samples")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    save(fig, out / "07_jaccard_heatmap.png")


def plot_shared_all_three(ov: pd.DataFrame, out: Path):
    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(ov))
    width = 0.25
    series = [
        ("from bam", "pct_shared_all_three_from_bam", "shared_all_three_from_bam", "bam_total", COLORS["bam"]),
        ("from fastq", "pct_shared_all_three_from_fastq", "shared_all_three_from_fastq", "fastq_total", COLORS["fastq"]),
        ("from none", "pct_shared_all_three_from_none", "shared_all_three_from_none", "none_total", COLORS["none"]),
    ]
    for i, (label, pct_col, cnt_col, tot_col, color) in enumerate(series):
        if pct_col in ov.columns:
            vals = ov[pct_col]
        else:
            vals = 100.0 * ov[cnt_col] / ov[tot_col]
        ax.bar(x + (i - 1) * width, vals, width, label=label, color=color)
    ax.set_xticks(x)
    ax.set_xticklabels(ov["sample"], rotation=35, ha="right")
    ax.set_ylim(0, 105)
    style_ax(ax, "Peaks overlapping both other methods", "% of that method's peaks")
    ax.legend(frameon=False)
    save(fig, out / "08_pct_shared_all_three.png")


def venn_subsets_from_row(row: pd.Series) -> tuple:
    """
    Approximate 3-set Venn region sizes from directional overlaps.
    Uses min() of reciprocal overlaps for pairwise/triple intersections.
    """
    nA, nB, nC = int(row["bam_total"]), int(row["fastq_total"]), int(row["none_total"])
    AB = min(int(row["bam_overlap_fastq"]), int(row["fastq_overlap_bam"]))
    AC = min(int(row["bam_overlap_none"]), int(row["none_overlap_bam"]))
    BC = min(int(row["fastq_overlap_none"]), int(row["none_overlap_fastq"]))
    ABC = min(
        int(row["shared_all_three_from_bam"]),
        int(row["shared_all_three_from_fastq"]),
        int(row["shared_all_three_from_none"]),
    )
    only_AB = max(AB - ABC, 0)
    only_AC = max(AC - ABC, 0)
    only_BC = max(BC - ABC, 0)
    only_A = max(nA - only_AB - only_AC - ABC, 0)
    only_B = max(nB - only_AB - only_BC - ABC, 0)
    only_C = max(nC - only_AC - only_BC - ABC, 0)
    # matplotlib_venn order: (Ab, aB, AB, abC, AbC, aBC, ABC)
    return (only_A, only_B, only_AB, only_C, only_AC, only_BC, ABC)


def plot_venn_mean(ov: pd.DataFrame, out: Path):
    # Average region sizes across samples
    regions = np.array([venn_subsets_from_row(row) for _, row in ov.iterrows()], dtype=float)
    mean_regions = tuple(regions.mean(axis=0))
    fig, ax = plt.subplots(figsize=(7, 6))
    v = venn3(
        subsets=mean_regions,
        set_labels=("bam", "fastq", "none"),
        ax=ax,
        set_colors=(COLORS["bam"], COLORS["fastq"], COLORS["none"]),
        alpha=0.55,
    )
    for text in v.set_labels or []:
        if text:
            text.set_fontsize(11)
    for text in v.subset_labels or []:
        if text:
            text.set_fontsize(8)
            # show integers
            try:
                text.set_text(f"{float(text.get_text()):.0f}")
            except ValueError:
                pass
    ax.set_title("Approximate peak-set Venn (mean across samples)\nintersection sizes from reciprocal overlaps")
    save(fig, out / "09_venn_mean.png")


def plot_venn_per_sample(ov: pd.DataFrame, out: Path):
    n = len(ov)
    ncols = 3
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 4 * nrows))
    axes = np.array(axes).reshape(-1)
    for i, (_, row) in enumerate(ov.iterrows()):
        ax = axes[i]
        subsets = venn_subsets_from_row(row)
        v = venn3(
            subsets=subsets,
            set_labels=("bam", "fastq", "none"),
            ax=ax,
            set_colors=(COLORS["bam"], COLORS["fastq"], COLORS["none"]),
            alpha=0.55,
        )
        for text in v.subset_labels or []:
            if text:
                text.set_fontsize(7)
                try:
                    text.set_text(f"{float(text.get_text()):.0f}")
                except ValueError:
                    pass
        for text in v.set_labels or []:
            if text:
                text.set_fontsize(9)
        ax.set_title(row["sample"], fontsize=10)
    for j in range(i + 1, len(axes)):
        axes[j].axis("off")
    fig.suptitle("Approximate peak-set Venn diagrams per sample", y=1.01, fontsize=13)
    save(fig, out / "10_venn_per_sample.png")


def plot_scatter_reads_vs_peaks(reads: pd.DataFrame, peaks: pd.DataFrame, out: Path):
    merged = reads.merge(peaks, on="sample")
    fig, ax = plt.subplots(figsize=(7, 5.5))
    for m in METHOD_ORDER:
        ax.scatter(
            merged[f"{m}_fragments"] / 1e6,
            merged[f"{m}_n_peaks"] / 1000,
            s=70,
            color=COLORS[m],
            label=METHOD_LABELS[m],
            alpha=0.85,
            edgecolors="white",
            linewidths=0.5,
        )
        for _, row in merged.iterrows():
            ax.annotate(
                row["sample"].replace("_REP", " R"),
                (row[f"{m}_fragments"] / 1e6, row[f"{m}_n_peaks"] / 1000),
                fontsize=6,
                alpha=0.55,
                xytext=(3, 3),
                textcoords="offset points",
            )
    style_ax(ax, "Library depth vs peak yield", "Peaks (thousands)", "Fragments (millions)")
    ax.legend(frameon=False)
    save(fig, out / "11_scatter_fragments_vs_peaks.png")


def plot_concordance_summary(ov: pd.DataFrame, out: Path):
    """Grouped bars: mean pairwise % of peaks retained in the other method."""
    metrics = {
        "bam→fastq": (ov["bam_overlap_fastq"] / ov["bam_total"] * 100).mean(),
        "bam→none": (ov["bam_overlap_none"] / ov["bam_total"] * 100).mean(),
        "fastq→bam": (ov["fastq_overlap_bam"] / ov["fastq_total"] * 100).mean(),
        "fastq→none": (ov["fastq_overlap_none"] / ov["fastq_total"] * 100).mean(),
        "none→bam": (ov["none_overlap_bam"] / ov["none_total"] * 100).mean(),
        "none→fastq": (ov["none_overlap_fastq"] / ov["none_total"] * 100).mean(),
    }
    fig, ax = plt.subplots(figsize=(9, 4.5))
    keys = list(metrics.keys())
    vals = list(metrics.values())
    bars = ax.bar(keys, vals, color="#72B7B2")
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5, f"{val:.1f}%", ha="center", fontsize=8)
    ax.set_ylim(0, 110)
    style_ax(ax, "Mean directional peak recovery (A→B = % of A peaks overlapping B)", "% of source peaks")
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels(keys, rotation=25, ha="right")
    save(fig, out / "12_directional_overlap_means.png")


def update_report(report: Path, plot_dir: Path, rel_prefix: str = "plots"):
    plots = sorted(plot_dir.glob("*.png"))
    if not plots:
        return
    section = [
        "",
        "## Figures",
        "",
        "Generated by `scripts/plot_benchmark_dedup_te_only.py`.",
        "",
    ]
    captions = {
        "01_fragments_by_sample.png": "Final Allo BAM fragments (read1) per sample.",
        "02_pct_reads_removed.png": "Percent fragments removed vs none.",
        "03_mean_fragments_and_removal.png": "Mean depth and mean % removal summary.",
        "04_peak_counts_by_sample.png": "MACS3 peak counts per sample.",
        "05_peak_delta_vs_baseline.png": "Peak count change vs none.",
        "06_jaccard_by_sample.png": "Peak-set Jaccard by sample.",
        "07_jaccard_heatmap.png": "Mean pairwise Jaccard heatmap.",
        "08_pct_shared_all_three.png": "Percent of peaks overlapping both other methods.",
        "09_venn_mean.png": "Approximate mean Venn of peak sets.",
        "10_venn_per_sample.png": "Approximate per-sample Venn diagrams.",
        "11_scatter_fragments_vs_peaks.png": "Fragments vs peak yield.",
        "12_directional_overlap_means.png": "Mean directional peak recovery.",
    }
    for p in plots:
        cap = captions.get(p.name, p.stem.replace("_", " "))
        section.append(f"### {cap}")
        section.append("")
        section.append(f"![{cap}]({rel_prefix}/{p.name})")
        section.append("")

    text = report.read_text()
    marker = "## Figures"
    if marker in text:
        text = text.split(marker)[0].rstrip() + "\n"
    report.write_text(text + "\n".join(section))
    print(f"Updated {report}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--benchmark-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "benchmark_dedup_te_only",
    )
    args = ap.parse_args()
    bdir = args.benchmark_dir
    plot_dir = bdir / "plots"

    reads = read_table(bdir / "read_counts.tsv")
    peaks = read_table(bdir / "peak_counts.tsv")
    ov = read_table(bdir / "peak_overlaps.tsv")

    # Stable sample order
    sample_order = list(reads["sample"])
    peaks = peaks.set_index("sample").loc[sample_order].reset_index()
    ov = ov.set_index("sample").loc[sample_order].reset_index()

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 12,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )

    plot_fragments_grouped(reads, plot_dir)
    plot_pct_removed(reads, plot_dir)
    plot_mean_removed_summary(reads, plot_dir)
    plot_peak_counts(peaks, plot_dir)
    plot_peak_delta(peaks, plot_dir)
    plot_jaccard_bars(ov, plot_dir)
    plot_jaccard_heatmap(ov, plot_dir)
    plot_shared_all_three(ov, plot_dir)
    plot_venn_mean(ov, plot_dir)
    plot_venn_per_sample(ov, plot_dir)
    plot_scatter_reads_vs_peaks(reads, peaks, plot_dir)
    plot_concordance_summary(ov, plot_dir)
    update_report(bdir / "REPORT.md", plot_dir)


if __name__ == "__main__":
    main()
