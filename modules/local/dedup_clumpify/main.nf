process DEDUP_CLUMPIFY {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data' :
        'community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9' }"

    input:
    tuple val(meta), path(reads1, stageAs: 'input1/*'), path(reads2, stageAs: 'input2/*')

    output:
    tuple val(meta), path('*.clumped.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')             , emit: log
    tuple val("${task.process}"), val('bbmap'), eval('bbversion.sh | grep -v "Duplicate cpuset"'), emit: versions_bbmap, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // reads1 / reads2 may each contain multiple files (one per technical replicate /
    // resequencing run of the same library). They are concatenated so that Clumpify
    // deduplicates across the whole pooled (merged) library at once.
    if (meta.single_end) {
        """
        cat ${reads1} > ${prefix}.merged.fastq.gz

        clumpify.sh \\
            in=${prefix}.merged.fastq.gz \\
            out=${prefix}.clumped.fastq.gz \\
            threads=${task.cpus} \\
            ${args} \\
            &> ${prefix}.clumpify.log

        rm -f ${prefix}.merged.fastq.gz
        """
    } else {
        """
        cat ${reads1} > ${prefix}_R1.merged.fastq.gz
        cat ${reads2} > ${prefix}_R2.merged.fastq.gz

        clumpify.sh \\
            in1=${prefix}_R1.merged.fastq.gz \\
            in2=${prefix}_R2.merged.fastq.gz \\
            out1=${prefix}_1.clumped.fastq.gz \\
            out2=${prefix}_2.clumped.fastq.gz \\
            threads=${task.cpus} \\
            ${args} \\
            &> ${prefix}.clumpify.log

        rm -f ${prefix}_R1.merged.fastq.gz ${prefix}_R2.merged.fastq.gz
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_command = meta.single_end ?
        "echo '' | gzip > ${prefix}.clumped.fastq.gz" :
        "echo '' | gzip > ${prefix}_1.clumped.fastq.gz ; echo '' | gzip > ${prefix}_2.clumped.fastq.gz"
    """
    touch ${prefix}.clumpify.log
    ${output_command}
    """
}
