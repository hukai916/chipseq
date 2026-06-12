//
// Allocate multi-mapped reads with Allo on BAM_FILTER_BAMTOOLS output.
// Filter steps (bamtools filter, PE orphan removal) are not repeated here.
// Allo rewrites alignments, so coordinate-sort, index, and samtools QC are rerun.
// Emits the same channels as BAM_FILTER_BAMTOOLS for drop-in use downstream.
//

include { BAM_SORT_STATS_SAMTOOLS } from '../../nf-core/bam_sort_stats_samtools'
include { ALLO_REDUCE             } from '../../../modules/local/allo_reduce'

workflow BAM_ALLO_REDUCE {

    take:
    ch_bam   // channel: [ val(meta), path(bam) ] BAM_FILTER_BAMTOOLS output
    ch_fasta // channel: [ val(meta), path(fasta) ]

    main:

    //
    // Collate + Allo run inside ALLO_REDUCE
    //
    ALLO_REDUCE ( ch_bam )

    //
    // Re-sort, index, and QC (Allo SAM output is not coordinate-sorted)
    //
    BAM_SORT_STATS_SAMTOOLS (
        ALLO_REDUCE.out.bam,
        ch_fasta
    )

    emit:
    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), path(bam) ]
    bai      = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), path(bai) ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), path(stats) ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), path(idxstats) ]
}

//
// Wrapper for standalone runs (run_test.nf -entry RUN_TEST).
// Qualifies task names as RUN_BAM_ALLO_REDUCE:BAM_ALLO_REDUCE:PROCESS for modules.config.
//
workflow RUN_BAM_ALLO_REDUCE {

    take:
    ch_bam   // channel: [ val(meta), path(bam) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]

    main:
    BAM_ALLO_REDUCE(ch_bam, ch_fasta)
}
