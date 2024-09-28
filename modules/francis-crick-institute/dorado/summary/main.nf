process DORADO_SUMMARY {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/thecrick/pipetech_dorado:0.8.0-linux-x64'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    dorado summary \\
        $args \\
        $bam \\
        > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(dorado --version 2>&1)
    END_VERSIONS
    """
}
