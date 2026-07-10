//
// Index a (already sorted) BAM file and run samtools stats, flagstat and idxstats.
// Used when BAM-level deduplication is skipped (e.g. FASTQ-level dedup with Clumpify),
// so that downstream steps receive the same channel interface as the
// BAM_MARKDUPLICATES_* subworkflows.
//

include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../../nf-core/bam_stats_samtools/main'

workflow BAM_INDEX_STATS_SAMTOOLS {

    take:
    ch_bam   // channel: [ val(meta), path(bam) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]

    main:
    SAMTOOLS_INDEX ( ch_bam )

    ch_bam_bai = ch_bam.join(SAMTOOLS_INDEX.out.bai, by: [0])

    BAM_STATS_SAMTOOLS ( ch_bam_bai, ch_fasta )

    emit:
    bam      = ch_bam                            // channel: [ val(meta), path(bam) ]
    bai      = SAMTOOLS_INDEX.out.bai            // channel: [ val(meta), path(bai) ]
    metrics  = channel.empty()                   // channel: [ val(meta), path(metrics) ] (no dedup metrics)

    stats    = BAM_STATS_SAMTOOLS.out.stats      // channel: [ val(meta), path(stats) ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat   // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats   // channel: [ val(meta), path(idxstats) ]
}
