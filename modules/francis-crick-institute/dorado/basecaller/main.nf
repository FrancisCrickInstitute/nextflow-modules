process DORADO_BASECALLER {
    tag "$meta.id"
    label 'process_high'

    container 'docker.io/thecrick/pipetech_dorado:0.8.0-linux-x64'

    input:
    tuple val(meta), path("pod5s/*")
    path(bam)
    val(model)
    val(bc_kit)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def bc_kit_arg = bc_kit ? "--kit-name ${bc_kit}"  : ''
    def resume_bam = bam ? "--resume-from resume_pod5.bam" : ''

    """
    export LC_ALL=C
    RANDOM_ID=\$(head /dev/urandom | tr -dc 'a-zA-Z0-9' | head -c 8 || true)

    if [ -L "pod5.bam" ]; then
        mv pod5.bam resume_pod5.bam
    fi

    dorado basecaller \\
        $model \\
        pod5s/ \\
        $bc_kit_arg \\
        $resume_bam \\
        $args \\
        > ${prefix}_\${RANDOM_ID}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(dorado --version 2>&1)
    END_VERSIONS
    """
}
