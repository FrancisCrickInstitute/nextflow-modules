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
