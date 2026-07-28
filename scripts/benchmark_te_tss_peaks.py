#!/usr/bin/env python3
"""
TE-proximal peak analysis for dedup benchmarking.

For each method/sample, read HOMER annotatePeaks.txt (te_gtf_only runs), keep peaks
with Distance to TSS in [-2000, +2000], and:
  - count filtered peaks and unique TE loci (Nearest PromoterID)
  - bar-plot counts
  - Venn of TE locus sets across methods
  - write unique TE CSVs per method (and per-sample tables)

Requires: matplotlib, matplotlib-venn, pandas
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib_venn import venn3

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

# Map results dirs -> method keys (from params dedup_level)
DEFAULT_RUNS = {
    "bam": "results_test3_te_only",
    "fastq": "results_test4_te_only",
    "none": "results_test5_te_only",
}


def normalize_colnames(fieldnames: List[str]) -> Dict[str, str]:
    """Map logical names -> actual header strings."""
    mapping = {}
    for col in fieldnames:
        c = col.strip()
        if c.startswith("PeakID"):
            mapping["PeakID"] = col
        else:
            mapping[c] = col
    return mapping


def parse_distance(val: str) -> Optional[int]:
    if val is None:
        return None
    s = str(val).strip()
    if s == "" or s.upper() == "NA" or s == ".":
        return None
    try:
        return int(float(s))
    except ValueError:
        return None


def load_te_proximal_peaks(
    path: Path,
    tss_min: int = -2000,
    tss_max: int = 2000,
) -> Tuple[pd.DataFrame, Set[str]]:
    """Return filtered peak table and unique TE IDs."""
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        colmap = normalize_colnames(reader.fieldnames or [])
        required = ["Distance to TSS", "Nearest PromoterID"]
        for req in required:
            if req not in colmap:
                raise KeyError(f"{path}: missing column {req!r}; have {list(colmap)}")

        rows = []
        tes: Set[str] = set()
        for row in reader:
            dist = parse_distance(row[colmap["Distance to TSS"]])
            if dist is None or dist < tss_min or dist > tss_max:
                continue
            te = (row.get(colmap["Nearest PromoterID"]) or "").strip()
            if not te or te.upper() == "NA":
                entrez_col = colmap.get("Entrez ID")
                te = (row.get(entrez_col) or "").strip() if entrez_col else ""
            if not te or te.upper() == "NA":
                continue
            peak_col = colmap.get("PeakID")
            peak_id = row.get(peak_col, "") if peak_col else ""
            def g(name: str, default: str = "") -> str:
                key = colmap.get(name)
                return row.get(key, default) if key else default

            rows.append(
                {
                    "peak_id": peak_id,
                    "chr": g("Chr"),
                    "start": g("Start"),
                    "end": g("End"),
                    "strand": g("Strand"),
                    "peak_score": g("Peak Score"),
                    "annotation": g("Annotation"),
                    "distance_to_tss": dist,
                    "te_id": te,
                    "entrez_id": g("Entrez ID"),
                    "gene_name": g("Gene Name"),
                }
            )
            tes.add(te)
    return pd.DataFrame(rows), tes


def style_ax(ax, title: str, ylabel: str = ""):
    ax.set_title(title, fontsize=12, pad=10)
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle=":", alpha=0.4)


def savefig(fig, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    pct = path
    plt.close(fig)
    print(f"Wrote {pct}")


def venn_from_sets(sets: Dict[str, Set[str]]) -> tuple:
    a, b, c = (sets["bam"], sets["fastq"], sets["none"])
    only_a = a - b - c
    only_b = b - a - c
    only_c = c - a - b
    ab = (a & b) - c
    ac = (a & c) - b
    bc = (b & c) - a
    abc = a & b & c
    return (len(only_a), len(only_b), len(ab), len(only_c), len(ac), len(bc), len(abc))


def plot_venn(sets: Dict[str, Set[str]], title: str, path: Path):
    fig, ax = plt.subplots(figsize=(7, 6))
    subsets = venn_from_sets(sets)
    v = venn3(
        subsets=subsets,
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
            text.set_fontsize(9)
    ax.set_title(title)
    # Add set sizes in subtitle-ish text
    note = (
        f"|bam|={len(sets['bam'])}  |fastq|={len(sets['fastq'])}  |none|={len(sets['none'])}  "
        f"|union|={len(sets['bam']|sets['fastq']|sets['none'])}"
    )
    ax.text(0.5, -0.05, note, transform=ax.transAxes, ha="center", fontsize=9)
    savefig(fig, path)


def write_unique_te_csv(path: Path, te_ids: Iterable[str], meta: Optional[dict] = None):
    rows = sorted(set(te_ids))
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["te_id"] + (list(meta.keys()) if meta else []))
        for te in rows:
            w.writerow([te] + (list(meta.values()) if meta else []))
    print(f"Wrote {path} ({len(rows)} TEs)")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project-dir", type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Default: <project>/benchmark_dedup_te_only/te_tss",
    )
    ap.add_argument("--tss-min", type=int, default=-2000)
    ap.add_argument("--tss-max", type=int, default=2000)
    args = ap.parse_args()

    project = args.project_dir
    outdir = args.outdir or (project / "benchmark_dedup_te_only" / "te_tss")
    plot_dir = outdir / "plots"
    csv_dir = outdir / "unique_tes"
    outdir.mkdir(parents=True, exist_ok=True)
    csv_dir.mkdir(parents=True, exist_ok=True)

    # Discover samples from baseline
    base_peak_dir = project / DEFAULT_RUNS["none"] / "bowtie2/merged_library/macs3/narrow_peak"
    samples = sorted(
        p.name.replace("_peaks.annotatePeaks.txt", "")
        for p in base_peak_dir.glob("*_peaks.annotatePeaks.txt")
    )
    if not samples:
        raise SystemExit(f"No annotatePeaks.txt under {base_peak_dir}")

    count_rows = []
    # method -> sample -> set(te)
    te_by_method_sample: Dict[str, Dict[str, Set[str]]] = {m: {} for m in METHOD_ORDER}
    # method -> union TE set
    te_by_method: Dict[str, Set[str]] = {m: set() for m in METHOD_ORDER}
    # keep filtered peak tables optionally light: only write summary + unique lists

    for sample in samples:
        row = {"sample": sample}
        for method, run in DEFAULT_RUNS.items():
            ann = (
                project
                / run
                / "bowtie2/merged_library/macs3/narrow_peak"
                / f"{sample}_peaks.annotatePeaks.txt"
            )
            if not ann.exists():
                row[f"{method}_n_peaks_tss"] = None
                row[f"{method}_n_unique_tes"] = None
                te_by_method_sample[method][sample] = set()
                continue
            df, tes = load_te_proximal_peaks(ann, args.tss_min, args.tss_max)
            te_by_method_sample[method][sample] = tes
            te_by_method[method] |= tes
            row[f"{method}_n_peaks_tss"] = len(df)
            row[f"{method}_n_unique_tes"] = len(tes)

            # per-sample unique TE csv
            write_unique_te_csv(
                csv_dir / f"{method}__{sample}__unique_tes.csv",
                tes,
                meta={"method": method, "sample": sample},
            )

            # also save filtered peak table (compact)
            peaks_out = outdir / "filtered_peaks" / f"{method}__{sample}__tss{args.tss_min}_{args.tss_max}.csv"
            peaks_out.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(peaks_out, index=False)

        # pairwise TE Jaccard / recovery for this sample
        for a, b in [("bam", "fastq"), ("bam", "none"), ("fastq", "none")]:
            sa, sb = te_by_method_sample[a][sample], te_by_method_sample[b][sample]
            inter = sa & sb
            union = sa | sb
            row[f"jaccard_te_{a}_{b}"] = round(len(inter) / len(union), 4) if union else 0.0
            row[f"{a}_tes_in_{b}"] = len(sa & sb)
            row[f"{a}_pct_tes_in_{b}"] = round(100 * len(sa & sb) / len(sa), 2) if sa else None
            row[f"{b}_pct_tes_in_{a}"] = round(100 * len(sa & sb) / len(sb), 2) if sb else None
        shared3 = (
            te_by_method_sample["bam"][sample]
            & te_by_method_sample["fastq"][sample]
            & te_by_method_sample["none"][sample]
        )
        row["n_tes_shared_all_three"] = len(shared3)
        union3 = (
            te_by_method_sample["bam"][sample]
            | te_by_method_sample["fastq"][sample]
            | te_by_method_sample["none"][sample]
        )
        row["n_tes_union"] = len(union3)
        row["pct_tes_shared_all_three_of_union"] = (
            round(100 * len(shared3) / len(union3), 2) if union3 else None
        )
        count_rows.append(row)

    counts = pd.DataFrame(count_rows)
    counts.to_csv(outdir / "te_tss_counts.tsv", sep="\t", index=False)

    # Method-level unique TE CSVs (union across samples)
    for method in METHOD_ORDER:
        write_unique_te_csv(
            csv_dir / f"{method}__ALL_SAMPLES__unique_tes.csv",
            te_by_method[method],
            meta={"method": method, "sample": "ALL"},
        )

    # Exclusive / shared TE lists (pooled across samples)
    a, b, c = te_by_method["bam"], te_by_method["fastq"], te_by_method["none"]
    partitions = {
        "only_bam": a - b - c,
        "only_fastq": b - a - c,
        "only_none": c - a - b,
        "bam_fastq_not_none": (a & b) - c,
        "bam_none_not_fastq": (a & c) - b,
        "fastq_none_not_bam": (b & c) - a,
        "all_three": a & b & c,
        "union": a | b | c,
    }
    for name, s in partitions.items():
        write_unique_te_csv(csv_dir / f"pooled__{name}.csv", s)

    # Summary JSON
    def mean_col(col):
        return round(float(counts[col].mean()), 2)

    summary = {
        "tss_window": [args.tss_min, args.tss_max],
        "te_id_column": "Nearest PromoterID",
        "n_samples": len(samples),
        "samples": samples,
        "mean_n_peaks_tss": {m: mean_col(f"{m}_n_peaks_tss") for m in METHOD_ORDER},
        "mean_n_unique_tes": {m: mean_col(f"{m}_n_unique_tes") for m in METHOD_ORDER},
        "pooled_n_unique_tes": {m: len(te_by_method[m]) for m in METHOD_ORDER},
        "pooled_venn": {
            "only_bam": len(partitions["only_bam"]),
            "only_fastq": len(partitions["only_fastq"]),
            "only_none": len(partitions["only_none"]),
            "bam_fastq_not_none": len(partitions["bam_fastq_not_none"]),
            "bam_none_not_fastq": len(partitions["bam_none_not_fastq"]),
            "fastq_none_not_bam": len(partitions["fastq_none_not_bam"]),
            "all_three": len(partitions["all_three"]),
            "union": len(partitions["union"]),
        },
        "mean_jaccard_te": {
            "bam_vs_fastq": mean_col("jaccard_te_bam_fastq"),
            "bam_vs_none": mean_col("jaccard_te_bam_none"),
            "fastq_vs_none": mean_col("jaccard_te_fastq_none"),
        },
        "mean_pct_shared_all_three_of_union": mean_col("pct_tes_shared_all_three_of_union"),
    }

    # Difference assessment
    mean_tes = summary["mean_n_unique_tes"]
    base = mean_tes["none"]
    pct_change = {
        m: round(100 * (mean_tes[m] - base) / base, 2) if base else None
        for m in ("bam", "fastq")
    }
    pooled = summary["pooled_venn"]
    exclusive_frac = {
        "only_bam_of_union": round(100 * pooled["only_bam"] / pooled["union"], 2),
        "only_fastq_of_union": round(100 * pooled["only_fastq"] / pooled["union"], 2),
        "only_none_of_union": round(100 * pooled["only_none"] / pooled["union"], 2),
        "all_three_of_union": round(100 * pooled["all_three"] / pooled["union"], 2),
    }
    bam_vs_fastq_j = summary["mean_jaccard_te"]["bam_vs_fastq"]
    bam_fastq_big = bool(
        bam_vs_fastq_j < 0.9
        or exclusive_frac["only_bam_of_union"] >= 2
        or exclusive_frac["only_fastq_of_union"] >= 2
    )
    max_abs_change = max(abs(pct_change["bam"] or 0), abs(pct_change["fastq"] or 0))
    dedup_vs_none_big = bool(max_abs_change >= 10 or exclusive_frac["only_none_of_union"] >= 10)
    summary["pct_change_unique_tes_vs_none"] = pct_change
    summary["exclusive_fractions_of_pooled_union"] = exclusive_frac
    summary["bam_vs_fastq_big_difference"] = bam_fastq_big
    summary["dedup_vs_none_big_difference"] = dedup_vs_none_big
    summary["big_difference"] = bam_fastq_big or dedup_vs_none_big
    summary["interpretation"] = (
        f"Mean unique TE loci near TSS vs none: bam {pct_change['bam']}%, fastq {pct_change['fastq']}%. "
        f"Mean TE Jaccard bam↔fastq={bam_vs_fastq_j}, "
        f"bam↔none={summary['mean_jaccard_te']['bam_vs_none']}, "
        f"fastq↔none={summary['mean_jaccard_te']['fastq_vs_none']}. "
        f"Pooled union: {exclusive_frac['all_three_of_union']}% in all three; "
        f"none-only {exclusive_frac['only_none_of_union']}%; "
        f"bam-only {exclusive_frac['only_bam_of_union']}%; "
        f"fastq-only {exclusive_frac['only_fastq_of_union']}%. "
        + (
            "bam vs fastq: small difference. "
            if not bam_fastq_big
            else "bam vs fastq: notable difference. "
        )
        + (
            "Dedup vs none: moderate difference (none recovers extra TE loci)."
            if dedup_vs_none_big
            else "Dedup vs none: small difference."
        )
    )

    with (outdir / "summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2)

    # -------- Plots --------
    plt.rcParams.update({"font.size": 10, "figure.facecolor": "white", "axes.facecolor": "white"})

    # Bar: filtered peak counts
    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(counts))
    width = 0.25
    for i, m in enumerate(METHOD_ORDER):
        ax.bar(x + (i - 1) * width, counts[f"{m}_n_peaks_tss"], width, label=METHOD_LABELS[m], color=COLORS[m])
    ax.set_xticks(x)
    ax.set_xticklabels(counts["sample"], rotation=35, ha="right")
    style_ax(ax, f"Peaks with Distance to TSS in [{args.tss_min}, {args.tss_max}]", "Peak count")
    ax.legend(frameon=False)
    savefig(fig, plot_dir / "13_te_tss_peak_counts.png")

    # Bar: unique TE counts
    fig, ax = plt.subplots(figsize=(10, 5))
    for i, m in enumerate(METHOD_ORDER):
        ax.bar(x + (i - 1) * width, counts[f"{m}_n_unique_tes"], width, label=METHOD_LABELS[m], color=COLORS[m])
    ax.set_xticks(x)
    ax.set_xticklabels(counts["sample"], rotation=35, ha="right")
    style_ax(ax, f"Unique TE loci with a peak within TSS [{args.tss_min}, {args.tss_max}]", "Unique TE count")
    ax.legend(frameon=False)
    savefig(fig, plot_dir / "14_te_tss_unique_te_counts.png")

    # Mean summary bars
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].bar(METHOD_ORDER, [summary["mean_n_peaks_tss"][m] for m in METHOD_ORDER], color=[COLORS[m] for m in METHOD_ORDER])
    style_ax(axes[0], "Mean TSS-proximal peaks", "Peaks")
    axes[0].set_xticks(range(3))
    axes[0].set_xticklabels([METHOD_LABELS[m] for m in METHOD_ORDER], rotation=20, ha="right")
    axes[1].bar(METHOD_ORDER, [summary["mean_n_unique_tes"][m] for m in METHOD_ORDER], color=[COLORS[m] for m in METHOD_ORDER])
    style_ax(axes[1], "Mean unique TE loci", "TEs")
    axes[1].set_xticks(range(3))
    axes[1].set_xticklabels([METHOD_LABELS[m] for m in METHOD_ORDER], rotation=20, ha="right")
    savefig(fig, plot_dir / "15_te_tss_means.png")

    # Jaccard bars
    fig, ax = plt.subplots(figsize=(10, 5))
    pairs = [
        ("bam↔fastq", "jaccard_te_bam_fastq"),
        ("bam↔none", "jaccard_te_bam_none"),
        ("fastq↔none", "jaccard_te_fastq_none"),
    ]
    for i, (lab, col) in enumerate(pairs):
        ax.bar(x + (i - 1) * width, counts[col], width, label=lab)
    ax.set_xticks(x)
    ax.set_xticklabels(counts["sample"], rotation=35, ha="right")
    ax.set_ylim(0, 1.05)
    style_ax(ax, "TE-locus set Jaccard (TSS-proximal peaks)", "Jaccard")
    ax.legend(frameon=False)
    savefig(fig, plot_dir / "16_te_tss_jaccard_by_sample.png")

    # Pooled Venn
    plot_venn(
        te_by_method,
        f"Unique TE loci near TSS [{args.tss_min}, {args.tss_max}] (pooled across samples)",
        plot_dir / "17_te_tss_venn_pooled.png",
    )

    # Per-sample Venns
    n = len(samples)
    ncols = 3
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 4 * nrows))
    axes = np.array(axes).reshape(-1)
    for i, sample in enumerate(samples):
        ax = axes[i]
        sets = {m: te_by_method_sample[m][sample] for m in METHOD_ORDER}
        subsets = venn_from_sets(sets)
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
        for text in v.set_labels or []:
            if text:
                text.set_fontsize(8)
        ax.set_title(sample, fontsize=10)
    for j in range(i + 1, len(axes)):
        axes[j].axis("off")
    fig.suptitle(f"TE-locus Venns per sample (TSS [{args.tss_min}, {args.tss_max}])", y=1.01)
    savefig(fig, plot_dir / "18_te_tss_venn_per_sample.png")

    # Markdown report (no tabulate dependency)
    simple = [
        f"# TE-proximal peak benchmark (Distance to TSS ∈ [{args.tss_min}, {args.tss_max}])\n",
        "TE ID = HOMER `Nearest PromoterID` from `*_peaks.annotatePeaks.txt`.\n",
        "## Method mapping\n",
        "| Method | Results dir |",
        "|---|---|",
        "| bam | `results_test3_te_only` |",
        "| fastq | `results_test4_te_only` |",
        "| none | `results_test5_te_only` |",
        "",
        "## Counts per sample\n",
        "| sample | none peaks | fastq peaks | bam peaks | none TEs | fastq TEs | bam TEs | Jaccard bam-fastq | shared all3 / union % |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for _, r in counts.iterrows():
        simple.append(
            f"| {r['sample']} | {r['none_n_peaks_tss']} | {r['fastq_n_peaks_tss']} | {r['bam_n_peaks_tss']} | "
            f"{r['none_n_unique_tes']} | {r['fastq_n_unique_tes']} | {r['bam_n_unique_tes']} | "
            f"{r['jaccard_te_bam_fastq']} | {r['pct_tes_shared_all_three_of_union']} |"
        )
    simple += [
        "",
        "## Summary\n",
        f"- Mean TSS-proximal peaks: none={summary['mean_n_peaks_tss']['none']}, "
        f"fastq={summary['mean_n_peaks_tss']['fastq']}, bam={summary['mean_n_peaks_tss']['bam']}",
        f"- Mean unique TE loci: none={summary['mean_n_unique_tes']['none']}, "
        f"fastq={summary['mean_n_unique_tes']['fastq']}, bam={summary['mean_n_unique_tes']['bam']}",
        f"- Mean TE Jaccard: {summary['mean_jaccard_te']}",
        f"- Pooled Venn: {summary['pooled_venn']}",
        f"- Exclusive fractions of pooled union: {exclusive_frac}",
        "",
        "## Does dedup make a big difference?\n",
        f"- **bam vs fastq:** {'Yes' if bam_fastq_big else 'No (small difference)'}",
        f"- **dedup vs none:** {'Yes (moderate)' if dedup_vs_none_big else 'No (small difference)'}",
        "",
        summary["interpretation"],
        "",
        "## Outputs\n",
        f"- `{outdir / 'te_tss_counts.tsv'}`",
        f"- `{outdir / 'summary.json'}`",
        f"- Unique TE CSVs under `{csv_dir}/` "
        "(per method×sample, method unions, and pooled Venn partitions)",
        f"- Filtered peak CSVs under `{outdir / 'filtered_peaks'}/`",
        "",
        "## Figures\n",
    ]
    for p in sorted(plot_dir.glob("*.png")):
        simple.append(f"### {p.stem}\n")
        simple.append(f"![{p.stem}](plots/{p.name})\n")

    (outdir / "REPORT_TE_TSS.md").write_text("\n".join(simple) + "\n")
    print(f"Wrote {outdir / 'REPORT_TE_TSS.md'}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
