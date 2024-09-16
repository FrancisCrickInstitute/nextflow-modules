/*
Subworkflow for basecalling and demultiplexing nanopore data
*/

include { SAMPLESHEET_CHECK                   } from '../../../../modules/francis-crick-institute/samplesheet/check/main'
include { DORADO_BASECALLER                   } from '../../../../modules/francis-crick-institute/dorado/basecaller/main'
include { DORADO_DEMUX                        } from '../../../../modules/francis-crick-institute/dorado/demux/main'
include { SAMTOOLS_MERGE as MERGE_BASECALLING } from '../../../../modules/nf-core/samtools/merge/main'

// Valid BC Kits
def valid_bc_kits = [
    "EXP-NBD103",
    "EXP-NBD104",
    "EXP-NBD114",
    "EXP-NBD196",
    "EXP-PBC001",
    "EXP-PBC096",
    "SQK-16S024",
    "SQK-16S114-24",
    "SQK-LWB001",
    "SQK-MLK111-96-XL",
    "SQK-MLK114-96-XL",
    "SQK-NBD111-24",
    "SQK-NBD111-96",
    "SQK-NBD114-24",
    "SQK-NBD114-96",
    "SQK-PBK004",
    "SQK-PCB109",
    "SQK-PCB110",
    "SQK-PCB111-24",
    "SQK-PCB114-24",
    "SQK-RAB201",
    "SQK-RAB204",
    "SQK-RBK001",
    "SQK-RBK004",
    "SQK-RBK110-96",
    "SQK-RBK111-24",
    "SQK-RBK111-96",
    "SQK-RBK114-24",
    "SQK-RBK114-96",
    "SQK-RLB001",
    "SQK-RPB004",
    "SQK-RPB114-24",
    "TWIST-ALL",
    "VSK-PTC001",
    "VSK-VMK001",
    "VSK-VMK004",
    "VSK-VPS001"
]

def find_bc_kit(run_dir, valid_bc_kits) {
    // Find the summary file if it exists
    def run_dir_path = new File(run_dir)
    def summary_file_name = run_dir_path.list().find { it.contains('final_summary') && it.endsWith('.txt') }
    if (!summary_file_name) {
        return null
    }

    // Load summary file content
    def summary_file = new File(run_dir_path, summary_file_name)
    def summary_content = summary_file.readLines().join('\n')

    // Find the first matching barcode kit in the summary file
    def bc_kit = valid_bc_kits.find { summary_content.contains(it) }
    return bc_kit
}

workflow ONT_DEMULTIPLEX {

    take:
    val_model         // Specify a short code like HAC or SUP
    val_bc_kit        // The barcode kit to demultiplex against
    val_check_barcode // Specifies if the run directory should be searched for a valid barcode kit
    val_bc_parse_pos  // The parse position to substring the barcode from in the meta data
    val_append_bc     // Appends bc-kit to sample matching
    val_batch_num     // Number of files in a dorado batch
    val_run_dir       // The string path of the nanopore run directory
    val_bam           // The string path of a BAM file to use instead of basecalling
    val_resume_bam    // The string path of a BAM file to resume basecalling from
    val_emit_bam      // Defines whether demux outputs a bam or fastq file
    val_samplesheet   // The string path of the samplesheet to parse for metadata

    main:

    // Init
    ch_versions = Channel.empty()

    // Check bc-kit is valid if supplied by user and assign if it is
    def bc_kit = null
    if(val_bc_kit != null && !(val_bc_kit in valid_bc_kits)) {
        exit 1, "Invalid barcode kit specified: ${val_bc_kit}"
    }
    else if (val_bc_kit != null) {
        bc_kit = val_bc_kit
    }

    // Try to resolve bc_kit
    if(val_run_dir != null && val_bc_kit == null) {
        bc_kit = find_bc_kit(val_run_dir, valid_bc_kits)
        if(bc_kit) { 
            log.info("Barcode Kit found from summary file: ${bc_kit}")
        }
    }

    //
    // MODULE: Check samplesheet
    //
    ch_samplesheet = Channel.from(file(val_samplesheet, checkIfExists: true))
    SAMPLESHEET_CHECK (
        ch_samplesheet
    )

    //
    // CHANNEL: Parse samplesheet into metadata
    //
    ch_sample_meta = SAMPLESHEET_CHECK.out.csv
        .splitCsv (header:true, sep:",")
        .map {
            it.group = it.group.replaceAll(" ", "_").toLowerCase()
            it.user = it.user.replaceAll(" ", "_").toLowerCase()
            if(val_bc_parse_pos != null) {
                it.barcode = "barcode" + it.barcode.substring(val_bc_parse_pos, val_bc_parse_pos + 2)
            }
            if(bc_kit != null && val_append_bc == true) {
                it.barcode = bc_kit + "_" + it.barcode
            }
            it
        }

    //
    // CHANNEL: Add all pod5 files
    //
    ch_pod5_files = Channel.empty()
    if(val_run_dir != null) {
        ch_pod5_files         = Channel.fromPath("${val_run_dir}/pod5/*.pod5")
        ch_pod5_files_pass    = Channel.fromPath("${val_run_dir}/pod5_pass/*.pod5")
        ch_pod5_files_fail    = Channel.fromPath("${val_run_dir}/pod5_fail/*.pod5")
        ch_pod5_files_skipped = Channel.fromPath("${val_run_dir}/pod5_skipped/*.pod5")
        ch_pod5_files         = ch_pod5_files_pass.mix(ch_pod5_files_fail).mix(ch_pod5_files_skipped).mix(ch_pod5_files)

        //
        // CHANNEL: Collate pod5 files into batches
        //
        ch_pod5_files = ch_pod5_files
            .collate(val_batch_num)
            .map{ [[ id: it[0].simpleName.substring(0, 26) ], it ] }
    }

    //
    // CHANNEL: Add bam file to a channel if it exists
    //
    ch_bam = Channel.empty()
    if (val_bam != null) {
        ch_bam = Channel.from(file(val_bam, checkIfExists: true))
            .map{ [ [ id: it.simpleName ], it ] }
    }

    //
    // CHANNEL: Add resume bam file to a channel if it exists
    //
    ch_resume_bam = Channel.empty()
    if (val_resume_bam != null) {
        ch_resume_bam = Channel.from(file(val_resume_bam, checkIfExists: true))
            .map{ [ [ id: it.simpleName ], it ] }
    }

    // Only run if a bam file wasnt supplied
    if(val_bam == null) {
        //
        // MODULE: Generate a bam file using pod5 files and any supplied bam to resume from
        //
        DORADO_BASECALLER (
            ch_pod5_files,
            val_resume_bam ? ch_resume_bam.map{it[1]} : [],
            val_model,
            bc_kit ?: []
        )
        ch_versions = ch_versions.mix(DORADO_BASECALLER.out.versions)
        ch_bam      = DORADO_BASECALLER.out.bam

        //
        // CHANNEL: Create basecalling merge channels
        //
        ch_bc_merge = ch_bam
            .collect{ it[1] }
            .branch {
                tomerge: it.size() > 1
                    return [[ id: it[0].simpleName.substring(0, 26) ], it ]
                pass: true
                    return [[ id: it[0].simpleName.substring(0, 26) ], it ]
            }

        //
        // MODULE: Merged basecalled bams if required
        //
        MERGE_BASECALLING (
            ch_bc_merge.tomerge,
            [[],[]],
            [[],[]]
        )
        ch_versions = ch_versions.mix(MERGE_BASECALLING.out.versions)
        ch_bam      = MERGE_BASECALLING.out.bam.mix(ch_bc_merge.pass)
    }

    //
    // CHANNEL: Add sample meta to main data channel
    //
    ch_demux_input = ch_sample_meta
        .collect()
        .map{[it]}
        .combine(ch_bam.collect{it[1]})

    //
    // MODULE: Generate demultiplexed bam or fastq files
    //
    DORADO_DEMUX (
        ch_demux_input,
        val_emit_bam
    )
    ch_versions    = ch_versions.mix(DORADO_DEMUX.out.versions)
    ch_demux_bam   = DORADO_DEMUX.out.bam
    ch_demux_fastq = DORADO_DEMUX.out.fastq

    //
    // CHANNEL: Merge metadata to the demultiplexed fastq file
    //
    ch_demux_fastq_meta = ch_sample_meta
        .map { [it.id, it] }
        .join( ch_demux_fastq.map{it[1]}.flatten().map{ [ it.simpleName, it ] } )
        .map { [ it[1], it[2] ] }

    //
    // CHANNEL: Merge metadata to the demultiplexed bam file
    //
    ch_demux_bam_meta = ch_sample_meta
        .map { [it.id, it] }
        .join( ch_demux_bam.map{it[1]}.flatten().map{ [ it.simpleName, it ] } )
        .map { [ it[1], it[2] ] }

    //
    // CHANNEL: Check for no matches
    //
    meta_list = ch_sample_meta.toSortedList()
    if(val_emit_bam) {
        ch_demux_bam_meta.count().map{
            if(it == 0) {
                log.error(meta_list.toString())
                exit 1, "Sample output channel empty"
            }
        }
    }
    else {
        ch_demux_fastq_meta.count().map{
            if(it == 0) {
                log.error(meta_list.toString())
                exit 1, "Sample output channel empty"
            }
        }
    }

    emit:
    pod5        = ch_pod5_files       // channel: [ path(pod5) ]
    bam         = ch_bam              // channel: [ val(meta), path(bam) ]
    demux_fastq = ch_demux_fastq_meta // channel: [ val(meta), path(fastq) ]
    demux_bam   = ch_demux_bam_meta   // channel: [ val(meta), path(bam) ]
    versions    = ch_versions         // channel: path(versions.yml)
}
