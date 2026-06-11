"""Shared BAM validation helpers for dedup scripts."""

import sys
from pathlib import Path
import pysam

NAME_SORT_ORDER = "queryname"


def is_name_sorted(bam_path):
    """Return True if the BAM is name-sorted."""
    pysam.set_verbosity(0)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        so = bam.header.get("HD", {}).get("SO")

    if so == NAME_SORT_ORDER:
        return True
    if so in ("coordinate", "unsorted"):
        return False
    if so is not None:
        return False

    return _has_non_decreasing_query_names(bam_path)


def prepare_name_sorted_bam(bam_path, threads=1):
    """
    Return (bam_to_read, temp_path).

    If already name-sorted, returns (bam_path, None).
    Otherwise name-sorts to a sibling temp file via pysam.sort and returns
    (sorted_path, sorted_path) so the caller can delete temp_path when done.
    """
    if is_name_sorted(bam_path):
        return bam_path, None

    path = Path(bam_path)
    sorted_path = path.with_name(f"{path.stem}.namesort.tmp{path.suffix}")

    print(
        f"Input BAM is not name-sorted; sorting with pysam.sort: "
        f"{bam_path} -> {sorted_path}",
        file=sys.stderr,
    )
    pysam.set_verbosity(0)
    pysam.sort("-n", "-@", str(threads), "-o", str(sorted_path), str(bam_path))
    return str(sorted_path), str(sorted_path)


def ensure_name_sorted(bam_path):
    """
    Exit with an error unless the BAM is name-sorted.

    Checks @HD SO:queryname first. If the SO tag is missing, verifies that
    query names are in non-decreasing order (all alignments per read grouped).
    """
    if is_name_sorted(bam_path):
        return

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        so = bam.header.get("HD", {}).get("SO")

    if so is None:
        _exit_unsorted_query_names(bam_path)
    _exit_not_name_sorted(bam_path, so)


def _has_non_decreasing_query_names(bam_path):
    prev = None
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            qname = read.query_name
            if prev is not None and qname < prev:
                return False
            prev = qname
    return True


def _exit_not_name_sorted(bam_path, so):
    label = so if so is not None else "not set"
    print(
        f"Error: input BAM must be name-sorted (@HD SO:queryname); found SO:{label!r}.",
        file=sys.stderr,
    )
    print(f"Fix: samtools sort -n -o name_sorted.bam {bam_path}", file=sys.stderr)
    sys.exit(1)


def _exit_unsorted_query_names(bam_path):
    prev = None
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            qname = read.query_name
            if prev is not None and qname < prev:
                print(
                    "Error: BAM header lacks SO:queryname and read names are not "
                    f"name-sorted (saw {qname!r} after {prev!r}).",
                    file=sys.stderr,
                )
                print(f"Fix: samtools sort -n -o name_sorted.bam {bam_path}", file=sys.stderr)
                sys.exit(1)
            prev = qname
