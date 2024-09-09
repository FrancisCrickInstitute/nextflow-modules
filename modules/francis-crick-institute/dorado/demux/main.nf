process DORADO_DEMUX {
    label 'process_high'

    container 'docker.io/thecrick/pipetech_dorado:0.7.3-linux-x64'

    input:
    tuple val(meta), path(bam)
    val emit_bam

    output:
    tuple val(meta), path("*.bam")     , emit: bam , optional: true
    tuple val(meta), path("*.fastq.gz"), emit: fastq , optional: true
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def emit = emit_bam ? '' : '--emit-fastq'
    def bc_map_bash = meta.collect { "${it.barcode}=${it.id}" }.join(' ')

    """
    dorado demux \\
        --output-dir ./ \\
        --no-classify \\
        $emit \\
        $args \\
        $bam

    if ls ./*.fastq 1> /dev/null 2>&1; then
        for file in ./*.fastq; do
            gzip \"\$file\"
        done
    fi

    for f in ./*.bam ./*.fastq.gz; do
        if [ -f "\$f" ]; then
            base_name=\$(basename "\$f" | sed 's/\\.bam//' | sed 's/\\.fastq\\.gz//')
            for pair in ${bc_map_bash}; do
                barcode=\$(echo \$pair | cut -d'=' -f1)
                sample_id=\$(echo \$pair | cut -d'=' -f2)

                if [[ "\$f" == *.fastq.gz ]]; then
                    ext="fastq.gz"
                else
                    ext="bam"
                fi

                if [ "\$base_name" == "\$barcode" ]; then
                    mv "\$f" "\$sample_id.\$ext"
                    break
                fi
            done
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(dorado --version 2>&1)
    END_VERSIONS
    """
}
