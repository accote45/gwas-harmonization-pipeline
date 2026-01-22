#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.json_file = null
params.output_dir = null
params.dual_build = true  // Set false for single build output only
params.help = false

def helpMessage() {
    log.info"""
    ========================================
    GWAS Harmonization Pipeline
    ========================================
    
    Usage:
      nextflow run harmonize.nf --json_file <json> --output_dir <dir>
    
    Required:
      --json_file      Path to JSON configuration file
      --output_dir     Output directory for QC'd data
    
    Optional:
      --dual_build     Generate dual-build merged output [default: true]
    
    Example (from qc_scripts directory):
      nextflow run harmonize.nf \\
        --json_file json_configs/cad_nonbiobank.json \\
        --output_dir /sc/arion/projects/paul_oreilly/data/GWASs/NonBiobanks/qced_data
    ========================================
    """.stripIndent()
}

if (params.help) { helpMessage(); exit 0 }
if (!params.json_file || !params.output_dir) {
    log.error "ERROR: --json_file and --output_dir are required"
    helpMessage()
    exit 1
}

process HARMONIZE {
    tag "${trait}_${build_type}"
    publishDir "${params.output_dir}/${trait}/logs/${build_type}", mode: 'copy', pattern: "logs/*"
    
    input:
    tuple val(trait), val(source), val(original_build), val(build_type), path(json_file)
    
    output:
    tuple val(trait), val(source), path("${trait}_${build_name}.txt.gz"), emit: harmonized
    path "logs/*", optional: true
    
    script:
    build_name = (build_type == 'original') ? 
        (original_build == '37' ? 'hg19' : 'hg38') :
        (original_build == '37' ? 'hg38' : 'hg19')
    
    """
    Rscript ${projectDir}/harmonize_gwas.R \\
        ${json_file} \\
        ${trait}_${build_name}.txt.gz \\
        ${build_type}
    """
}

process MERGE_BUILDS {
    tag "${trait}"
    publishDir "${params.output_dir}/${trait}", mode: 'copy'
    
    input:
    tuple val(trait), val(source), path(hg19_file), path(hg38_file)
    
    output:
    path "${trait}_munged.txt.gz"
    
    """
    Rscript ${projectDir}/merge_builds.R \\
        ${trait} \\
        ${hg19_file} \\
        ${hg38_file} \\
        ${trait}_munged.txt.gz
    """
}

workflow {
    // Parse JSON directly
    def json = new groovy.json.JsonSlurper().parse(file(params.json_file))
    def config = json instanceof List ? json[0] : json
    def trait = config.trait
    def source = config.source ?: 'Unknown'
    def build = config.build ?: '37'
    
    // Create input channel
    input_ch = Channel.of(
        [trait, source, build, 'original', file(params.json_file)],
        [trait, source, build, 'lifted', file(params.json_file)]
    )
    
    // Run harmonization
    HARMONIZE(input_ch)
    
    // Optionally merge builds
    if (params.dual_build) {
        merge_input = HARMONIZE.out.harmonized
            .groupTuple(by: [0, 1])
            .map { t, s, files -> 
                // Sort files by name to ensure hg19 comes before hg38
                def sorted = files.sort { it.name }
                def hg19 = sorted.find { it.name.contains('hg19') }
                def hg38 = sorted.find { it.name.contains('hg38') }
                tuple(t, s, hg19, hg38)
            }
        
        MERGE_BUILDS(merge_input)
    }
}