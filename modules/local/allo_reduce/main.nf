process ALLO_REDUCE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity'
        ? 'docker://hukai916/allo:0.2'
        : 'hukai916/allo:0.2' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val("${task.process}"), val('allo'), eval("python -c \"import importlib.metadata as m; print(m.version('bio-allo'))\" 2>/dev/null || echo unknown"), topic: versions, emit: versions_allo

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.mLb.allo"
    def is_control = meta.is_control != null
        ? meta.is_control
        : (meta.control instanceof Boolean ? meta.control : !meta.control)
    def args = task.ext.args ?: ''
    if (is_control) {
        args = "${args} --random".trim()
    }
    def seq_mode = meta.single_end ? 'se' : 'pe'
    """
    samtools collate \\
        -@ ${task.cpus} \\
        -o collated.bam \\
        ${bam}

    allo \\
        collated.bam \\
        -seq ${seq_mode} \\
        -o ${prefix} \\
        -p ${task.cpus} \\
        ${args}

    if [ -f ${prefix}.bam ]; then
        :
    elif [ -f ${prefix} ]; then
        mv ${prefix} ${prefix}.bam
    else
        echo "ERROR: Allo did not produce expected output (tried ${prefix}.bam and ${prefix})" >&2
        exit 1
    fi

    rm -f collated.bam
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.mLb.allo"
    """
    touch ${prefix}.bam
    """
}
