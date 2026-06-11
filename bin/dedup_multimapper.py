#!/usr/bin/env python3

# Optimized deduplication for BAM files with multi-mapped reads.
# Compared to 03a: uses a signature inverted index and exact frozenset cache
# to avoid scanning the full fingerprint_db for every read group.

__version__ = "1.0.0"

import argparse
import sys
from collections import defaultdict
from pathlib import Path
import pysam
from dedup_bam_checks import prepare_name_sorted_bam

def get_fragment_signature(read):
    """Generates a coordinate signature for a single alignment record."""
    chrom = read.reference_name
    if chrom is None:
        return None

    if read.is_paired:
        frag_start = min(read.reference_start, read.next_reference_start)
        frag_len = abs(read.template_length)
        return f"{chrom}:{frag_start}:{frag_len}"
    else:
        return f"{chrom}:{read.reference_start}:{read.reference_end}"

class FingerprintIndex:
    """Indexes unique read-group signatures for fast duplicate lookup."""

    def __init__(self):
        self.fingerprint_db = {}   # read_id -> frozenset of signatures
        self.sig_index = defaultdict(set)  # signature -> read_ids containing it
        self.exact_index = {}      # frozenset -> representative read_id

    def find_duplicate(self, current_signatures, representative_id, min_overlap_pct, log_file=None):
        """
        Return True if current_signatures overlaps any indexed group at or above min_overlap_pct. Logic matches 03a.
        """
        if not current_signatures:
            return False

        total_current = len(current_signatures)

        matched_id = self.exact_index.get(current_signatures)
        if matched_id is not None:
            if log_file:
                log_file.write(
                    f"{representative_id}\t{matched_id}\t{total_current}\t"
                    f"{total_current}\t100.00\n"
                )
            return True

        counts = defaultdict(int)
        for sig in current_signatures:
            for uid in self.sig_index.get(sig, ()):
                counts[uid] += 1

        for uid, shared_coords in counts.items():
            overlap_pct = (shared_coords / total_current) * 100.0
            if overlap_pct >= min_overlap_pct:
                if log_file:
                    log_file.write(
                        f"{representative_id}\t{uid}\t{total_current}\t"
                        f"{shared_coords}\t{overlap_pct:.2f}\n"
                    )
                return True

        return False

    def register_unique(self, representative_id, current_signatures):
        """Add a non-duplicate read group to all indexes."""
        self.fingerprint_db[representative_id] = current_signatures
        self.exact_index[current_signatures] = representative_id
        for sig in current_signatures:
            self.sig_index[sig].add(representative_id)


def collect_signatures(alignments):
    """Build a frozenset of fragment signatures for a read group."""
    signatures = set()
    for read in alignments:
        sig = get_fragment_signature(read)
        if sig:
            signatures.add(sig)
    return frozenset(signatures)


def process_read_group(alignments, index, min_overlap_pct, log_file=None):
    """Determine duplicate status and update index for unique groups."""
    current_signatures = collect_signatures(alignments)
    if not current_signatures:
        return False

    representative_id = alignments[0].query_name
    if index.find_duplicate(current_signatures, representative_id, min_overlap_pct, log_file):
        return True

    index.register_unique(representative_id, current_signatures)
    return False


def write_read_group(group, outfile, mode, is_dup):
    if is_dup and mode == "drop":
        return
    for read in group:
        if is_dup and mode == "mark":
            read.flag |= 1024  # 0x400 duplicate flag
        outfile.write(read)


def coordinate_sort_bam(input_bam, output_bam, threads=1):
    """Coordinate-sort a BAM with pysam.sort (samtools sort)."""
    pysam.set_verbosity(0)
    print(
        f"Coordinate-sorting output: {input_bam} -> {output_bam}",
        file=sys.stderr,
    )
    pysam.sort("-o", str(output_bam), "-@", str(threads), str(input_bam))


def main():
    parser = argparse.ArgumentParser(
        description="Deduplicate multi-mapped ChIP-seq alignments (optimized index)."
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s Version: {__version__}",
    )
    parser.add_argument("-i", "--input", required=True, help="Input name-sorted BAM")
    parser.add_argument("-o", "--output", required=True, help="Output BAM")
    parser.add_argument(
        "-p", "--overlap-pct", type=float, default=100.0,
        help="Minimum matching %% (Default: 100.0)",
    )
    parser.add_argument(
        "-m", "--mode", choices=["drop", "mark"], default="mark",
        help="Action for duplicates: 'drop' deletes rows, 'mark' sets 0x400 flag (Default: mark)",
    )
    parser.add_argument("-l", "--log", default=None, help="Optional duplicate log path")
    parser.add_argument(
        "--sort-threads", type=int, default=1,
        help="Threads for BAM sorting (name-sort input, coordinate-sort output) (Default: 1)",
    )
    parser.add_argument(
        "--coord-sort",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Coordinate-sort output BAM after dedup (default: on). Use --no-coord-sort for name order.",
    )
    args = parser.parse_args()

    pysam.set_verbosity(0)  # silence htslib missing-.bai notice; index not used for streaming

    input_bam, name_sort_tmp = prepare_name_sorted_bam(args.input, threads=args.sort_threads)

    output_path = Path(args.output)
    coord_sort_tmp = None
    if args.coord_sort:
        dedup_bam = output_path.with_name(f"{output_path.stem}.dedup.tmp{output_path.suffix}")
        coord_sort_tmp = dedup_bam
    else:
        dedup_bam = output_path

    index = FingerprintIndex()
    total_groups = 0
    duplicate_groups = 0

    log_fh = open(args.log, "w") if args.log else None
    if log_fh:
        log_fh.write(
            "Duplicate_Read_ID\tMatched_Unique_Read_ID\t"
            "Total_Coordinates\tShared_Coordinates\tOverlap_Percentage\n"
        )

    try:
        with pysam.AlignmentFile(input_bam, "rb") as infile, \
             pysam.AlignmentFile(str(dedup_bam), "wb", template=infile) as outfile:

            current_query = None
            current_group = []

            for read in infile:
                if read.query_name != current_query:
                    if current_group:
                        total_groups += 1
                        is_dup = process_read_group(
                            current_group, index, args.overlap_pct, log_fh
                        )
                        if is_dup:
                            duplicate_groups += 1
                        write_read_group(current_group, outfile, args.mode, is_dup)

                    current_group = [read]
                    current_query = read.query_name
                else:
                    current_group.append(read)

            if current_group:
                total_groups += 1
                is_dup = process_read_group(
                    current_group, index, args.overlap_pct, log_fh
                )
                if is_dup:
                    duplicate_groups += 1
                write_read_group(current_group, outfile, args.mode, is_dup)

        if args.coord_sort:
            coordinate_sort_bam(coord_sort_tmp, output_path, threads=args.sort_threads)

    finally:
        pct = (duplicate_groups / total_groups) * 100 if total_groups else 0.0
        summary_lines = [
            f"Deduplication Complete [{args.mode.upper()} MODE].",
            f"Total Unique Molecules Processed: {total_groups}",
            f"Duplicate Molecules Handled:     {duplicate_groups} ({pct:.2f}%)",
        ]
        if args.coord_sort:
            summary_lines.append(f"Output coordinate-sorted: {output_path}")
        for line in summary_lines:
            print(line)
        if log_fh:
            log_fh.write("\n" + "\n".join(summary_lines) + "\n")
            log_fh.close()
        if name_sort_tmp:
            Path(name_sort_tmp).unlink(missing_ok=True)
        if coord_sort_tmp:
            Path(coord_sort_tmp).unlink(missing_ok=True)

if __name__ == "__main__":
    main()
