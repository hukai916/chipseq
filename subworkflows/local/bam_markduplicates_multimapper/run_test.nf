#!/usr/bin/env nextflow
/*
 * Standalone test entry for BAM_MARKDUPLICATES_MULTIMAPPER.
 *
 * Run from subworkflows/local/bam_markduplicates_multimapper/ directory:
    // cd /path/to/02_chipseq/subworkflows/local/bam_markduplicates_multimapper/
    // nextflow run run_test.nf -c test2.config -profile mamba,lsf,singularity -resume
*/

include { BAM_MARKDUPLICATES_MULTIMAPPER } from './main'

workflow {
    if (!params.test_bam) {
        error('Set --test_bam to a merged coordinate-sorted BAM (e.g. PICARD_MERGESAMFILES output)')
    }

    ch_reads = Channel.of([
        [ id: params.test_sample_id, single_end: false ],
        file(params.test_bam, checkIfExists: true),
    ])

    ch_fasta = Channel.of([
        [ id: 'genome' ],
        file(params.fasta, checkIfExists: true),
    ])

    ch_fai = Channel.of([
        [ id: 'genome' ],
        file(params.fai, checkIfExists: true),
    ])

    BAM_MARKDUPLICATES_MULTIMAPPER(ch_reads, ch_fasta, ch_fai)
}
