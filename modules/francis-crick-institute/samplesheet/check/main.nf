process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'

    container "docker.io/python:3.11"

    input:
    path samplesheet

    output:
    path '*.csv'        , emit: csv
    when:
    task.ext.when == null || task.ext.when

    shell:
    sample       = samplesheet
    process_name = task.process
    output       = task.ext.output ?: 'samplesheet.valid.csv'
    template 'samplesheet_check.py'
}
