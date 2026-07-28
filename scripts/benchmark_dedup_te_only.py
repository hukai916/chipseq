#!/usr/bin/env python3
"""
Benchmark deduplication strategies on te_gtf_only ChIP-seq runs.

Compares:
  - Final usable read counts (post BAM-filter + Allo) vs no-dedup baseline
  - MACS3 narrowPeak counts per sample
  - Peak set overlaps (any base overlap via bedtools) between methods

Method labels are taken from each run's pipeline_info/params_*.json (dedup_level).
Requires: bedtools on PATH.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence


FLAGSTAT_TOTAL_RE = re.compile(r"^(\d+) \+ \d+ in total")
FLAGSTAT_PRIMARY_RE = re.compile(r"^(\d+) \+ \d+ primary$")
FLAGSTAT_READ1_RE = re.compile(r"^(\d+) \+ \d+ read1")


@dataclass(frozen=True)
class MethodRun:
    key: str
    label: str
    dedup_level: str
    root: Path


def latest_params(root: Path) -> dict:
    files = sorted(root.glob("pipeline_info/params_*.json"))
    if not files:
        raise FileNotFoundError(f"No params_*.json under {root}/pipeline_info")
    with files[-1].open() as fh:
        return json.load(fh)


def parse_flagstat(path: Path) -> dict:
    text = path.read_text()

    def _first(pat: re.Pattern) -> Optional[int]:
        for line in text.splitlines():
            m = pat.match(line)
            if m:
                return int(m.group(1))
        return None

    total = _first(FLAGSTAT_TOTAL_RE)
    primary = _first(FLAGSTAT_PRIMARY_RE)
    read1 = _first(FLAGSTAT_READ1_RE)
    if read1 is not None:
        fragments = read1
    elif primary is not None:
        fragments = primary // 2
    elif total is not None:
        fragments = total // 2
    else:
        fragments = None
    return {
        "total_reads": total,
        "primary_reads": primary,
        "read1": read1,
        "fragments": fragments,
    }


def find_final_flagstat(merged_stats: Path, sample: str) -> Optional[Path]:
    for name in (
        f"{sample}.mLb.clN.allo.sorted.bam.flagstat",
        f"{sample}.mLb.clN.sorted.bam.flagstat",
    ):
        p = merged_stats / name
        if p.exists():
            return p
    return None


def count_lines(path: Path) -> int:
    n = 0
    with path.open() as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                n += 1
    return n


def narrowpeak_to_bed6(src: Path, dst: Path) -> int:
    """Write coordinate-sorted BED6 from narrowPeak; return peak count."""
    unsorted = dst.with_suffix(dst.suffix + ".unsorted")
    n = 0
    with src.open() as fin, unsorted.open("w") as fout:
        for i, line in enumerate(fin):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            chrom, start, end = parts[0], parts[1], parts[2]
            name = parts[3] if len(parts) > 3 else f"peak_{i}"
            score = parts[4] if len(parts) > 4 else "0"
            strand = parts[5] if len(parts) > 5 else "."
            fout.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")
            n += 1
    with dst.open("w") as out:
        subprocess.run(
            ["bedtools", "sort", "-i", str(unsorted)],
            check=True,
            stdout=out,
        )
    unsorted.unlink(missing_ok=True)
    return n


def bedtools_count_u(a: Path, b: Path) -> int:
    """Number of intervals in A that overlap B (any base)."""
    proc = subprocess.run(
        ["bedtools", "intersect", "-u", "-a", str(a), "-b", str(b)],
        check=True,
        capture_output=True,
        text=True,
    )
    return sum(1 for line in proc.stdout.splitlines() if line.strip())


def bedtools_count_u_both(a: Path, b: Path, c: Path) -> int:
    """Number of intervals in A that overlap both B and C."""
    step1 = subprocess.run(
        ["bedtools", "intersect", "-u", "-a", str(a), "-b", str(b)],
        check=True,
        capture_output=True,
        text=True,
    )
    proc = subprocess.run(
        ["bedtools", "intersect", "-u", "-a", "stdin", "-b", str(c)],
        input=step1.stdout,
        check=True,
        capture_output=True,
        text=True,
    )
    return sum(1 for line in proc.stdout.splitlines() if line.strip())


def bedtools_jaccard_bp(a: Path, b: Path) -> float:
    """Genomic base-pair Jaccard from bedtools jaccard."""
    proc = subprocess.run(
        ["bedtools", "jaccard", "-a", str(a), "-b", str(b)],
        check=True,
        capture_output=True,
        text=True,
    )
    lines = [ln for ln in proc.stdout.splitlines() if ln.strip()]
    if len(lines) < 2:
        return 0.0
    # header: intersection union jaccard n_intersections
    cols = lines[1].split("\t")
    return float(cols[2])


def peak_set_jaccard(n_a: int, n_b: int, a_in_b: int, b_in_a: int) -> float:
    shared = min(a_in_b, b_in_a)
    union = n_a + n_b - shared
    return round(shared / union, 4) if union else 0.0


def discover_peak_samples(peak_dir: Path) -> List[str]:
    return [
        p.name.replace("_peaks.narrowPeak", "")
        for p in sorted(peak_dir.glob("*_peaks.narrowPeak"))
    ]


def write_csv(path: Path, rows: List[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore", delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def pct(num: Optional[float], den: Optional[float]) -> Optional[float]:
    if num is None or den is None or den == 0:
        return None
    return round(100.0 * num / den, 2)


def mean(vals: Iterable[Optional[float]]) -> Optional[float]:
    xs = [v for v in vals if v is not None]
    return round(sum(xs) / len(xs), 2) if xs else None


def md_table(rows: Sequence[dict], cols: Sequence[str]) -> str:
    out = [
        "| " + " | ".join(cols) + " |",
        "|" + "|".join(["---"] * len(cols)) + "|",
    ]
    for r in rows:
        out.append(
            "| "
            + " | ".join("" if r.get(c) is None else str(r.get(c)) for c in cols)
            + " |"
        )
    return "\n".join(out)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project-dir", type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument("--outdir", type=Path, default=None)
    ap.add_argument(
        "--runs",
        nargs="+",
        default=["results_test3_te_only", "results_test4_te_only", "results_test5_te_only"],
    )
    ap.add_argument("--baseline-dedup", default="none")
    args = ap.parse_args()

    if not shutil.which("bedtools"):
        raise SystemExit("bedtools not found on PATH (try: module load bedtools)")

    project = args.project_dir
    outdir = args.outdir or (project / "benchmark_dedup_te_only")
    outdir.mkdir(parents=True, exist_ok=True)

    methods: List[MethodRun] = []
    for name in args.runs:
        root = project / name
        params = latest_params(root)
        level = params.get("dedup_level", "unknown")
        mark = params.get("markduplicates", "")
        if level == "bam":
            label = f"bam ({mark})" if mark else "bam"
        elif level == "fastq":
            label = "fastq (clumpify)"
        elif level == "none":
            label = "none (baseline)"
        else:
            label = str(level)
        methods.append(MethodRun(key=level, label=label, dedup_level=level, root=root))

    key_counts: Dict[str, int] = defaultdict(int)
    unique_methods: List[MethodRun] = []
    for m in methods:
        key_counts[m.key] += 1
        key = m.key if key_counts[m.key] == 1 else f"{m.key}_{key_counts[m.key]}"
        unique_methods.append(MethodRun(key=key, label=m.label, dedup_level=m.dedup_level, root=m.root))
    methods = unique_methods

    baseline = next((m for m in methods if m.dedup_level == args.baseline_dedup), None)
    if baseline is None:
        raise SystemExit(f"Baseline dedup_level={args.baseline_dedup!r} not found")

    peak_dir_b = baseline.root / "bowtie2/merged_library/macs3/narrow_peak"
    samples = discover_peak_samples(peak_dir_b)
    if not samples:
        raise SystemExit(f"No narrowPeak files in {peak_dir_b}")

    # -------- Read counts --------
    read_rows = []
    for sample in samples:
        row: dict = {"sample": sample}
        counts = {}
        for m in methods:
            fs = find_final_flagstat(m.root / "bowtie2/merged_library/samtools_stats", sample)
            if fs is None:
                row[f"{m.key}_fragments"] = None
                row[f"{m.key}_total_reads"] = None
                counts[m.key] = None
                continue
            parsed = parse_flagstat(fs)
            row[f"{m.key}_fragments"] = parsed["fragments"]
            row[f"{m.key}_total_reads"] = parsed["total_reads"]
            counts[m.key] = parsed["fragments"]

        base_n = counts.get(baseline.key)
        row["baseline_fragments"] = base_n
        for m in methods:
            n = counts.get(m.key)
            if m.key == baseline.key:
                row[f"{m.key}_reads_removed_vs_baseline"] = 0
                row[f"{m.key}_pct_removed_vs_baseline"] = 0.0
            else:
                removed = (base_n - n) if (base_n is not None and n is not None) else None
                row[f"{m.key}_reads_removed_vs_baseline"] = removed
                row[f"{m.key}_pct_removed_vs_baseline"] = pct(removed, base_n)
        read_rows.append(row)
    write_csv(outdir / "read_counts.tsv", read_rows)

    # -------- Peak counts + overlaps (bedtools) --------
    method_keys = [m.key for m in methods]
    if len(method_keys) != 3:
        raise SystemExit(f"Expected 3 methods, got {method_keys}")
    a, b, c = method_keys

    peak_count_rows = []
    overlap_rows = []

    with tempfile.TemporaryDirectory(prefix="dedup_bench_") as tmp:
        tmpdir = Path(tmp)
        beds: Dict[str, Dict[str, Path]] = {m.key: {} for m in methods}
        counts: Dict[str, Dict[str, int]] = {m.key: {} for m in methods}

        for sample in samples:
            row = {"sample": sample}
            for m in methods:
                np_path = (
                    m.root
                    / "bowtie2/merged_library/macs3/narrow_peak"
                    / f"{sample}_peaks.narrowPeak"
                )
                if not np_path.exists():
                    row[f"{m.key}_n_peaks"] = None
                    continue
                bed = tmpdir / f"{m.key}__{sample}.bed"
                n = narrowpeak_to_bed6(np_path, bed)
                beds[m.key][sample] = bed
                counts[m.key][sample] = n
                row[f"{m.key}_n_peaks"] = n

            base_p = row.get(f"{baseline.key}_n_peaks")
            for m in methods:
                n = row.get(f"{m.key}_n_peaks")
                if base_p is not None and n is not None:
                    row[f"{m.key}_delta_peaks_vs_baseline"] = n - base_p
                    row[f"{m.key}_pct_peaks_vs_baseline"] = pct(n, base_p)
                else:
                    row[f"{m.key}_delta_peaks_vs_baseline"] = None
                    row[f"{m.key}_pct_peaks_vs_baseline"] = None
            peak_count_rows.append(row)

            if not all(sample in beds[k] for k in method_keys):
                continue

            ba, bb, bc = beds[a][sample], beds[b][sample], beds[c][sample]
            na, nb, nc = counts[a][sample], counts[b][sample], counts[c][sample]

            a_in_b = bedtools_count_u(ba, bb)
            a_in_c = bedtools_count_u(ba, bc)
            b_in_a = bedtools_count_u(bb, ba)
            b_in_c = bedtools_count_u(bb, bc)
            c_in_a = bedtools_count_u(bc, ba)
            c_in_b = bedtools_count_u(bc, bb)

            a_both = bedtools_count_u_both(ba, bb, bc)
            b_both = bedtools_count_u_both(bb, ba, bc)
            c_both = bedtools_count_u_both(bc, ba, bb)

            overlap_rows.append(
                {
                    "sample": sample,
                    f"{a}_total": na,
                    f"{b}_total": nb,
                    f"{c}_total": nc,
                    f"{a}_overlap_{b}": a_in_b,
                    f"{a}_overlap_{c}": a_in_c,
                    f"{b}_overlap_{a}": b_in_a,
                    f"{b}_overlap_{c}": b_in_c,
                    f"{c}_overlap_{a}": c_in_a,
                    f"{c}_overlap_{b}": c_in_b,
                    f"{a}_pct_overlap_{b}": pct(a_in_b, na),
                    f"{a}_pct_overlap_{c}": pct(a_in_c, na),
                    f"{b}_pct_overlap_{a}": pct(b_in_a, nb),
                    f"{b}_pct_overlap_{c}": pct(b_in_c, nb),
                    f"{c}_pct_overlap_{a}": pct(c_in_a, nc),
                    f"{c}_pct_overlap_{b}": pct(c_in_b, nc),
                    f"shared_all_three_from_{a}": a_both,
                    f"shared_all_three_from_{b}": b_both,
                    f"shared_all_three_from_{c}": c_both,
                    f"pct_shared_all_three_from_{a}": pct(a_both, na),
                    f"pct_shared_all_three_from_{b}": pct(b_both, nb),
                    f"pct_shared_all_three_from_{c}": pct(c_both, nc),
                    f"jaccard_peakset_{a}_{b}": peak_set_jaccard(na, nb, a_in_b, b_in_a),
                    f"jaccard_peakset_{a}_{c}": peak_set_jaccard(na, nc, a_in_c, c_in_a),
                    f"jaccard_peakset_{b}_{c}": peak_set_jaccard(nb, nc, b_in_c, c_in_b),
                    f"jaccard_bp_{a}_{b}": round(bedtools_jaccard_bp(ba, bb), 4),
                    f"jaccard_bp_{a}_{c}": round(bedtools_jaccard_bp(ba, bc), 4),
                    f"jaccard_bp_{b}_{c}": round(bedtools_jaccard_bp(bb, bc), 4),
                }
            )

    write_csv(outdir / "peak_counts.tsv", peak_count_rows)
    write_csv(outdir / "peak_overlaps.tsv", overlap_rows)

    summary = {
        "methods": [
            {
                "key": m.key,
                "label": m.label,
                "dedup_level": m.dedup_level,
                "results_dir": m.root.name,
                "is_baseline": m.key == baseline.key,
            }
            for m in methods
        ],
        "n_samples": len(samples),
        "samples": samples,
        "note": (
            "fragments = read1 from *.mLb.clN.allo.sorted.bam.flagstat. "
            "Peak overlap = bedtools intersect -u (any base). "
            "shared_all_three_from_X = peaks in X overlapping both other methods."
        ),
        "mean_fragments": {m.key: mean(r.get(f"{m.key}_fragments") for r in read_rows) for m in methods},
        "mean_reads_removed_vs_baseline": {
            m.key: mean(r.get(f"{m.key}_reads_removed_vs_baseline") for r in read_rows) for m in methods
        },
        "mean_pct_removed_vs_baseline": {
            m.key: mean(r.get(f"{m.key}_pct_removed_vs_baseline") for r in read_rows) for m in methods
        },
        "mean_n_peaks": {m.key: mean(r.get(f"{m.key}_n_peaks") for r in peak_count_rows) for m in methods},
        "mean_jaccard_peakset": {
            f"{a}_vs_{b}": mean(r.get(f"jaccard_peakset_{a}_{b}") for r in overlap_rows),
            f"{a}_vs_{c}": mean(r.get(f"jaccard_peakset_{a}_{c}") for r in overlap_rows),
            f"{b}_vs_{c}": mean(r.get(f"jaccard_peakset_{b}_{c}") for r in overlap_rows),
        },
        "mean_jaccard_bp": {
            f"{a}_vs_{b}": mean(r.get(f"jaccard_bp_{a}_{b}") for r in overlap_rows),
            f"{a}_vs_{c}": mean(r.get(f"jaccard_bp_{a}_{c}") for r in overlap_rows),
            f"{b}_vs_{c}": mean(r.get(f"jaccard_bp_{b}_{c}") for r in overlap_rows),
        },
        "mean_shared_all_three": {
            f"from_{a}": mean(r.get(f"shared_all_three_from_{a}") for r in overlap_rows),
            f"from_{b}": mean(r.get(f"shared_all_three_from_{b}") for r in overlap_rows),
            f"from_{c}": mean(r.get(f"shared_all_three_from_{c}") for r in overlap_rows),
        },
        "mean_pct_shared_all_three": {
            f"from_{a}": mean(r.get(f"pct_shared_all_three_from_{a}") for r in overlap_rows),
            f"from_{b}": mean(r.get(f"pct_shared_all_three_from_{b}") for r in overlap_rows),
            f"from_{c}": mean(r.get(f"pct_shared_all_three_from_{c}") for r in overlap_rows),
        },
    }
    with (outdir / "summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2)

    # -------- Markdown report --------
    lines = [
        "# Deduplication benchmark (te_gtf_only runs)\n",
        "## Methods\n",
        "| Results dir | dedup_level | Label | Role |",
        "|---|---|---|---|",
    ]
    for m in methods:
        role = "baseline" if m.key == baseline.key else "comparison"
        lines.append(f"| `{m.root.name}` | `{m.dedup_level}` | {m.label} | {role} |")
    lines += [
        "",
        "> **Note:** Labels follow `dedup_level` in each run's `pipeline_info/params_*.json` "
        "(not directory number). In these runs: "
        "`results_test3_te_only` = **bam**, `results_test4_te_only` = **fastq**, "
        "`results_test5_te_only` = **none**.",
        "",
        "## 1. Reads removed vs baseline (`none`)\n",
        "Metric: `read1` from final Allo BAM flagstat "
        "(`*.mLb.clN.allo.sorted.bam.flagstat`). "
        "Removed = baseline − method.\n",
        "### Fragments (read1)\n",
        md_table(read_rows, ["sample"] + [f"{m.key}_fragments" for m in methods]),
        "\n### Reads removed vs baseline\n",
        md_table(
            read_rows,
            ["sample"]
            + [f"{m.key}_reads_removed_vs_baseline" for m in methods if m.key != baseline.key],
        ),
        "\n### Percent removed vs baseline\n",
        md_table(
            read_rows,
            ["sample"]
            + [f"{m.key}_pct_removed_vs_baseline" for m in methods if m.key != baseline.key],
        ),
        "\n### Means\n",
        "| Method | Mean fragments | Mean removed | Mean % removed |",
        "|---|---:|---:|---:|",
    ]
    for m in methods:
        lines.append(
            f"| {m.label} | {summary['mean_fragments'][m.key]} | "
            f"{summary['mean_reads_removed_vs_baseline'][m.key]} | "
            f"{summary['mean_pct_removed_vs_baseline'][m.key]} |"
        )

    lines += [
        "\n## 2. MACS3 peak counts (te_gtf_only annotation runs)\n",
        md_table(peak_count_rows, ["sample"] + [f"{m.key}_n_peaks" for m in methods]),
        "\n### Mean peak counts\n",
        "| Method | Mean n_peaks | Mean Δ vs baseline |",
        "|---|---:|---:|",
    ]
    for m in methods:
        mean_delta = mean(r.get(f"{m.key}_delta_peaks_vs_baseline") for r in peak_count_rows)
        lines.append(f"| {m.label} | {summary['mean_n_peaks'][m.key]} | {mean_delta} |")

    lines += [
        "\n## 3. Peak overlaps (any base, bedtools)\n",
        "`X_overlap_Y` = peaks in X overlapping ≥1 bp of a peak in Y. "
        "`shared_all_three_from_X` = peaks in X overlapping both other methods.\n",
        md_table(
            overlap_rows,
            [
                "sample",
                f"{a}_total",
                f"{b}_total",
                f"{c}_total",
                f"{a}_overlap_{b}",
                f"{a}_overlap_{c}",
                f"{b}_overlap_{a}",
                f"{b}_overlap_{c}",
                f"{c}_overlap_{a}",
                f"{c}_overlap_{b}",
                f"shared_all_three_from_{a}",
                f"shared_all_three_from_{b}",
                f"shared_all_three_from_{c}",
                f"pct_shared_all_three_from_{a}",
                f"pct_shared_all_three_from_{b}",
                f"pct_shared_all_three_from_{c}",
                f"jaccard_peakset_{a}_{b}",
                f"jaccard_peakset_{a}_{c}",
                f"jaccard_peakset_{b}_{c}",
                f"jaccard_bp_{a}_{b}",
                f"jaccard_bp_{a}_{c}",
                f"jaccard_bp_{b}_{c}",
            ],
        ),
        "\n### Mean peak-set Jaccard\n",
        "| Pair | Mean Jaccard |",
        "|---|---:|",
    ]
    for k, v in summary["mean_jaccard_peakset"].items():
        lines.append(f"| {k} | {v} |")
    lines += [
        "\n### Mean genomic (bp) Jaccard\n",
        "| Pair | Mean Jaccard |",
        "|---|---:|",
    ]
    for k, v in summary["mean_jaccard_bp"].items():
        lines.append(f"| {k} | {v} |")
    lines += [
        "\n### Mean peaks shared with both other methods\n",
        "| Perspective | Mean count | Mean % of that method |",
        "|---|---:|---:|",
        f"| from_{a} | {summary['mean_shared_all_three'][f'from_{a}']} | {summary['mean_pct_shared_all_three'][f'from_{a}']} |",
        f"| from_{b} | {summary['mean_shared_all_three'][f'from_{b}']} | {summary['mean_pct_shared_all_three'][f'from_{b}']} |",
        f"| from_{c} | {summary['mean_shared_all_three'][f'from_{c}']} | {summary['mean_pct_shared_all_three'][f'from_{c}']} |",
        "\n## Output files\n",
        f"- `{outdir / 'read_counts.tsv'}`",
        f"- `{outdir / 'peak_counts.tsv'}`",
        f"- `{outdir / 'peak_overlaps.tsv'}`",
        f"- `{outdir / 'summary.json'}`",
        f"- `{outdir / 'REPORT.md'}`",
    ]

    report_path = outdir / "REPORT.md"
    report_path.write_text("\n".join(lines) + "\n")
    print(f"Wrote report: {report_path}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
