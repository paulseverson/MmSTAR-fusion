#!/usr/bin/env nextflow


/*
================================================================================
                         SET UP CONFIGURATION VARIABLES
================================================================================
*/

// Show help message
if (params.help) exit 0, helpMessage()

running_tools = []
visualization_tools = []
reference = [
    star_fusion: false
]

// Check if genome exists in the config file
if (params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

if (!Channel.fromPath(params.genomes_base, checkIfExists: true)) {exit 1, "Directory ${params.genomes_base} doesn't exist."}

params.star_fusion_ref = params.genome ? params.genomes[params.genome].star_fusion_ref ?: null : null

if (!params.star_index && (!params.fasta && !params.gtf)) exit 1, "Either specify STAR-INDEX or Fasta and GTF!"

if (params.star_fusion) {
    running_tools.add("STAR-Fusion")
    reference.star_fusion = Channel.value(file(params.star_fusion_ref)).ifEmpty{exit 1, "Star-Fusion reference directory not found!"}
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)

/*
 * Create a channel for input read files
 */
if(params.readPaths) {
    if(params.single_end) {
        Channel.from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty{exit 1, "params.readPaths was empty - no input files supplied" }
            .into{ch_read_files_fastqc; read_files_multiqc; read_files_star_fusion}
    } else {
        Channel.from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty{exit 1, "params.readPaths was empty - no input files supplied" }
            .into{ch_read_files_fastqc; read_files_multiqc; read_files_star_fusion}
    }
} else {
    Channel.fromFilePairs( params.reads, size: params.single_end ? 1 : 2 )
        .ifEmpty{exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --single_end on the command line." }
        .into{ch_read_files_fastqc; read_files_multiqc; read_files_star_fusion}
}

/*
================================================================================
                                PRINTING SUMMARY
================================================================================
*/

// Header log info
log.info nfcoreHeader()
def summary = [:]
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Fasta Ref']    = params.fasta
summary['GTF Ref']      = params.gtf
summary['STAR Index']   = params.star_index ? params.star_index : 'Not specified, building'
summary['Fusion tools']        = running_tools.size() == 0 ? 'None' : running_tools.join(", ")
summary['Data Type']        = params.single_end ? 'Single-End' : 'Paired-End'
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']   = params.outdir
summary['Launch dir']   = workflow.launchDir
summary['Working dir']  = workflow.workDir
summary['Script dir']   = workflow.projectDir
summary['User']         = workflow.userName
if(workflow.profile == 'awsbatch') {
   summary['AWS Region']    = params.awsregion
   summary['AWS Queue']     = params.awsqueue
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-rnafusion-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/rnafusion Workflow Summary'
    section_href: 'https://github.com/nf-core/rnafusion'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * STAR-Fusion
 */
process star_fusion {
    tag "${sample}"
    label 'process_high'

    publishDir "${params.outdir}/tools/Star-Fusion/${sample}", mode: 'copy'

    input:
        set val(sample), file(reads) from read_files_star_fusion
        file(reference) from reference.star_fusion
        file(star_index) from ch_star_index

    output:
        set val(sample), file("${sample}_star-fusion.tsv") optional true into star_fusion_fusions
        file("*.{tsv,txt}") into star_fusion_output

    when: params.star_fusion || (params.star_fusion && params.debug)

    script:
    def avail_mem = task.memory ? "--limitBAMsortRAM ${task.memory.toBytes() - 100000000}" : ''
    option = params.single_end ? "--left_fq ${reads[0]}" : "--left_fq ${reads[0]} --right_fq ${reads[1]}"
    def extra_params = params.star_fusion_opt ? params.star_fusion_opt : ''
    """
    STAR-Fusion \\
        --genome_lib_dir ${reference} \\
        --left_fq ${reads}[0] \\
        --right_fq ${reads}[1] \\
        --CPU ${task.cpus} \\
        -O . ${extra_params}

    mv star-fusion.fusion_predictions.tsv ${sample}_star-fusion.tsv
    mv star-fusion.fusion_predictions.abridged.tsv ${sample}_abridged.tsv
    """
}

star_fusion_fusions = star_fusion_fusions.dump(tag:'star_fusion_fusions')


/*************************************************************
 * Quality check & software verions
 ************************************************************/

/*
 * FastQC
 */
process fastqc {
    tag "$name"
    label 'process_medium'
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: { filename ->
                      filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
                }

    input:
    set val(name), file(reads) from ch_read_files_fastqc

    output:
    file "*_fastqc.{zip,html}" into ch_fastqc_results

    when: !params.debug

    script:
    """
    fastqc --quiet --threads $task.cpus $reads
    """
}

/*
 * MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file (multiqc_config) from ch_multiqc_config
    file (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
    file ('fastqc/*') from ch_fastqc_results.collect().ifEmpty([])
    file (fusions_mq) from summary_fusions_mq.collect().ifEmpty([])
    file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")

    output:
    file "*multiqc_report.html" into ch_multiqc_report
    file "*_data"
    file "multiqc_plots"

    when: !params.debug

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc -f $rtitle $rfilename $custom_config_file .
    """
}

