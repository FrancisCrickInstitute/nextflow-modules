//
// Dump pipeline parameters in a json file
//
def dump_parameters(workflow, params) {
    def timestamp = new java.util.Date().format( 'yyyy-MM-dd_HH-mm-ss')
    def filename = "params_${timestamp}.json"
    def temp_pf = new File(workflow.launchDir.toString(), ".${filename}")
    def jsonStr =  groovy.json.JsonOutput.toJson(params)
    temp_pf.text =  groovy.json.JsonOutput.prettyPrint(jsonStr)

    nextflow.extension.FilesEx.copyTo(temp_pf.toPath(), "${params.outdir}/pipeline_info/params_${timestamp}.json")
    temp_pf.delete()
}

//
// Dump channel meta into CSV
//
def dump_meta(meta, path) {
    def csvFile = new File(path)
    csvFile.parentFile.mkdirs()
    csvFile.withWriter { writer ->
        def headers = meta[0].keySet()
        writer.writeLine(headers.join(','))

        meta.each { map ->
            def row = headers.collect { map[it] }
            writer.writeLine(row.join(','))
        }
    }
}

//
// Generate a workflow complete file
//
def workflow_complete_summary(workflow, path) {
    def outputFile = new File(path)
    outputFile.parentFile.mkdirs()
    outputFile.withWriter { writer ->
        writer.writeLine("complete\tduration")
        writer.writeLine(workflow.complete.toString() + "\t" + workflow.duration.toString())
    }
}

//
// Generate a random ID
//
def gen_id(length) {
    def chars = (('A'..'Z') + ('a'..'z') + ('0'..'9')).join()
    def random = new java.security.SecureRandom()
    def sb = new StringBuilder(length)
    (1..length).each {
        sb.append(chars[random.nextInt(chars.length())])
    }
    return sb.toString()
}
