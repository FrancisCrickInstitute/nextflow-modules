include { WISECONDORX_CONVERT } from '../../../modules/nf-core/wisecondorx/convert/main'
include { WISECONDORX_PREDICT } from '../../../modules/nf-core/wisecondorx/predict/main'

workflow ONT_DEMULTIPLEX {

    take:
    val_dorado_auto_model // Specify a short code like HAC or SUP
    val_dorado_model      // Specify an exact model which will override the previous

    // ch_bam          // channel: [ val(meta), path(bam), path(bai) ]
    // ch_fasta        // channel: [ val(meta2), path(fasta) ]
    // ch_fai          // channel: [ val(meta3), path(fai) ]
    // ch_ref          // channel: [ val(meta4), path(reference) ]
    // ch_blacklist    // channel: [ val(meta5), path(blacklist) ]

    main:

    ch_versions = Channel.empty()


    emit:
    // aberrations_bed = WISECONDORX_PREDICT.out.aberrations_bed   // channel: [ val(meta), path(bed) ]
    // bins_bed        = WISECONDORX_PREDICT.out.bins_bed          // channel: [ val(meta), path(bed) ]
    // segments_bed    = WISECONDORX_PREDICT.out.segments_bed      // channel: [ val(meta), path(bed) ]
    // chr_statistics  = WISECONDORX_PREDICT.out.chr_statistics    // channel: [ val(meta), path(txt) ]
    // chr_plots       = WISECONDORX_PREDICT.out.chr_plots         // channel: [ val(meta), [ path(png), path(png), ... ] ]
    // genome_plot     = WISECONDORX_PREDICT.out.genome_plot       // channel: [ val(meta), path(png) ]

    // versions        = ch_versions                               // channel: path(versions.yml)
}


// Select dorado model
// Dorado models can be selected by auto selection (e.g. hac selects the latest compatible hac model) or by direct model selection.
// In the case of direct selection, we need to add the container path to the model.
dorado_model = params.dorado_auto_model
if(params.dorado_model) {
    dorado_model = "/home/" + params.dorado_model
}

// Check bc-kit
def dorado_bc_kit = null
if(params.dorado_bc_kit && !(params.dorado_bc_kit in params.bc_kits)) {
    exit 1, "Invalid barcode kit specified: ${params.dorado_bc_kit}"
}
// Extract barcode kit value from the summary file
if (!params.dorado_bc_kit) {
    // Find the summary file
    def summaryFileDir = new File(params.run_dir)
    def summaryFileName = summaryFileDir.list().find { it.contains('final_summary') && it.endsWith('.txt') }
    // Check if a summary file was found
    if (!summaryFileName) {
        exit 1, "No final summary file found in ${summaryFileDir}."
    }
    // Create a File object from the filename
    def summaryFile = new File(summaryFileDir, summaryFileName)
    // Read the entire content of the summary file as a single string
    def summaryContent = summaryFile.readLines().join('\n')

    // Find the first matching barcode kit in the summary file
    def extrapolatedBcKit = params.bc_kits.find { bc_kit ->
        summaryContent.contains(bc_kit)
    }

    // Set `params.dorado_bc_kit` to the found kit or null if no match
    dorado_bc_kit = extrapolatedBcKit ?: null
}

// Extract run_id
def runid = file(params.run_dir).name


  //
    // CHANNEL: Adding all pod5 files
    //
    ch_pod5_files         = Channel.fromPath("${params.run_dir}/pod5/*.pod5")
    ch_pod5_files_pass    = Channel.fromPath("${params.run_dir}/pod5_pass/*.pod5")
    ch_pod5_files_fail    = Channel.fromPath("${params.run_dir}/pod5_fail/*.pod5")
    ch_pod5_files_skipped = Channel.fromPath("${params.run_dir}/pod5_skipped/*.pod5")
    ch_pod5_files         = ch_pod5_files_pass.mix(ch_pod5_files_fail).mix(ch_pod5_files_skipped).mix(ch_pod5_files)
    ch_collected_pod5     = ch_pod5_files.collect().ifEmpty([])

    //
    // CHANNEL: Put all pod5 generated files and their corresponding sample IDs into a single channel 
    //
    ch_pod5_files = ch_pod5_files
        .collate(params.dorado_batch_num)
        .map{ [[ id: it[0].simpleName.substring(0, 26) ], it ] }

    //
    // CHANNEL: Adding bam files to a channel if it exists
    //
    ch_bam = Channel.empty()
    if (params.bam) {
        ch_bam = Channel.from(file(params.bam, checkIfExists: true))
            .map{ [ [ id: it.simpleName ], it ] }
    }

    //
    // CHANNEL: Load samplesheet
    //
    if(params.samplesheet) {
        ch_samplesheet = Channel.from(file(params.samplesheet))
    }

    //
    // SUBWORKFLOW: check input samplesheet and add relevant info to metadata
    //
    SAMPLESHEET_PARSE (
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(SAMPLESHEET_PARSE.out.versions)
    ch_meta     = SAMPLESHEET_PARSE.out.meta

    //
    // CHANNEL: extract run ID name and assign to metadata
    //
    ch_meta = ch_meta.map{
        it.run_id = runid 
        it.id = it.sample_id
        it.remove("sample_id")
        it
    }

    //
    // CHANNEL: Collect barcode names
    //
    ch_barcodes = ch_meta.map{ it.barcode }.toSortedList()

    if (params.run_basecaller) {
        //
        // MODULE: Generate a bam file using pod5 files and any supplied bam to resume from
        //
        DORADO_BASECALLER (
            ch_pod5_files,
            params.bam ? ch_bam.map{it[1]} : [],
            dorado_model,
            dorado_bc_kit ?: []
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

    if (params.run_demux) {
        //
        // MODULE: Generate demultiplexed bam or fastq files
        //
        DORADO_DEMUX (
            ch_bam
        )
        ch_versions    = ch_versions.mix(DORADO_DEMUX.out.versions)
        ch_demux_bam   = DORADO_DEMUX.out.bam
        ch_demux_fastq = DORADO_DEMUX.out.fastq

        //
        // CHANNEL: Merge metadata to the demultiplexed fastq file
        //
        ch_demux_fastq = ch_meta
            .map { [it.barcode, it] }
            .join( ch_demux_fastq.map{it[1]}.flatten().map{ [ it.simpleName, it ] } )
            .map { [ it[1], it[2] ] }

        //
        // CHANNEL: Merge metadata to the demultiplexed bam file
        //
        ch_demux_bam = ch_meta
            .map { [it.barcode, it] }
            .join( ch_demux_bam.map{it[1]}.flatten().map{ [ it.simpleName, it ] } )
            .map { [ it[1], it[2] ] }
    }