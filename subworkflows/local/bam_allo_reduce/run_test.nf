#!/usr/bin/env nextflow
/*
 * Standalone test entry for BAM_ALLO_REDUCE.
 *
 * Input should be a coordinate-sorted BAM with multimapping reads retained
 * (typically BAM_FILTER_BAMTOOLS output: *.mLb.clN.sorted.bam with keep_multi_map true).
 *
 * Uses RUN_BAM_ALLO_REDUCE wrapper so task names match modules.config withName selectors.
 *
 * From this directory:
 *   nextflow run run_test.nf -entry RUN_TEST -c test2.config -profile mamba,lsf,singularity -resume
 *
 * From pipeline root (02_chipseq/):
 *   nextflow run subworkflows/local/bam_allo_reduce/run_test.nf \
 *     -entry RUN_TEST \
 *     -c subworkflows/local/bam_allo_reduce/test2.config \
 *     -profile mamba,lsf,singularity \
 *     -resume
 */

include { RUN_BAM_ALLO_REDUCE } from './main'

workflow RUN_TEST {

    if (!params.test_bam) {
        error('Set --test_bam to a coordinate-sorted filtered BAM (e.g. *.mLb.clN.sorted.bam from BAM_FILTER_BAMTOOLS)')
    }

    ch_bam = Channel.of([
        [
            id: params.test_sample_id,
            single_end: params.test_single_end,
            is_control: params.test_is_control,
        ],
        file(params.test_bam, checkIfExists: true),
    ])

    ch_fasta = Channel.of([
        [ id: 'genome' ],
        file(params.fasta, checkIfExists: true),
    ])

    RUN_BAM_ALLO_REDUCE(ch_bam, ch_fasta)
}
