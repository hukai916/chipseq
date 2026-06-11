//
// Allocate multi-mapped reads with Allo, then coordinate-sort, index, and QC.
// Emits the same channels as BAM_FILTER_BAMTOOLS for drop-in use downstream.
//

include { SAMTOOLS_SORT                            } from '../../../modules/nf-core/samtools/sort'
include { BAM_SORT_STATS_SAMTOOLS                    } from '../../nf-core/bam_sort_stats_samtools'
include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_PE } from '../../nf-core/bam_sort_stats_samtools'
include { ALLO_REDUCE                              } from '../../../modules/local/allo_reduce'
include { BAM_REMOVE_ORPHANS                       } from '../../../modules/local/bam_remove_orphans'

workflow BAM_ALLO_REDUCE {

    take:
    ch_bam   // channel: [ val(meta), path(bam) ] coordinate-sorted filtered BAM
    ch_fasta // channel: [ val(meta), path(fasta) ]

    main:

    ch_bam
        .branch {
            meta, bam ->
                single_end: meta.single_end
                    return [ meta, bam ]
                paired_end: !meta.single_end
                    return [ meta, bam ]
        }
        .set { ch_bam_branched }

    //
    // Allo requires name-sorted (queryname) input (SE and PE)
    //
    SAMTOOLS_SORT (
        ch_bam_branched.single_end.mix(ch_bam_branched.paired_end),
        ch_fasta,
        ''
    )

    ALLO_REDUCE ( SAMTOOLS_SORT.out.bam )

    ALLO_REDUCE
        .out
        .bam
        .branch {
            meta, bam ->
                single_end: meta.single_end
                    return [ meta, bam ]
                paired_end: !meta.single_end
                    return [ meta, bam ]
        }
        .set { ch_allo_bam }

    //
    // SE: coordinate-sort, index, and samtools QC
    //
    BAM_SORT_STATS_SAMTOOLS (
        ch_allo_bam.single_end,
        ch_fasta
    )

    //
    // PE: remove orphans, then coordinate-sort, index, and samtools QC
    //
    BAM_REMOVE_ORPHANS (
        ch_allo_bam.paired_end
    )

    BAM_SORT_STATS_PE (
        BAM_REMOVE_ORPHANS.out.bam,
        ch_fasta
    )

    emit:
    name_bam = SAMTOOLS_SORT.out.bam                                                              // channel: [ val(meta), path(bam) ] name-sorted input to Allo
    bam      = BAM_SORT_STATS_PE.out.bam.mix(BAM_SORT_STATS_SAMTOOLS.out.bam)                     // channel: [ val(meta), path(bam) ]
    bai      = BAM_SORT_STATS_PE.out.bai.mix(BAM_SORT_STATS_SAMTOOLS.out.bai)                     // channel: [ val(meta), path(bai) ]
    stats    = BAM_SORT_STATS_PE.out.stats.mix(BAM_SORT_STATS_SAMTOOLS.out.stats)                 // channel: [ val(meta), path(stats) ]
    flagstat = BAM_SORT_STATS_PE.out.flagstat.mix(BAM_SORT_STATS_SAMTOOLS.out.flagstat)           // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_SORT_STATS_PE.out.idxstats.mix(BAM_SORT_STATS_SAMTOOLS.out.idxstats)           // channel: [ val(meta), path(idxstats) ]
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
