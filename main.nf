#!/usr/bin/env nextflow
/*
========================================================================================
    Cyclomics/CycloSeq Informed pipeline
========================================================================================
    Github : https://github.com/cyclomics/cycloseq
    Website: https://cyclomics.com
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    PARAMETER VALUES
========================================================================================
*/
// ### PARAMETERS
params.read_folder             = ""
params.read_pattern               = "**.{fq,fastq,fq.gz,fastq.gz}"
params.sequencing_summary_path = "${projectDir}/sequencing_summary*.txt"
params.backbone                   = "BBCS"
// Backbone file is used for custom backbones.
params.backbone_file  = ""
// reference indexes are expected to be in reference folder
params.reference = ""
params.output_dir = "$HOME/Data/CyclomicsSeq"
params.backbone_barcode = false


// method selection
params.consensus_method        = "cycas"

// Pipeline performance metrics
params.min_repeat_count = 3

if (params.backbone == "BB41") {
    backbone_file = "$projectDir/backbones/BB41.fasta"
}
else if (params.backbone == "BB41T") {
    backbone_file = "$projectDir/backbones/BB41T.fasta"
}
else if (params.backbone == "BB42") {
    backbone_file = "$projectDir/backbones/BB42.fasta"
}
else if (params.backbone == "BB43") {
    backbone_file = "$projectDir/backbones/BB43.fasta"
}
else if (params.backbone == "BB43b") {
    backbone_file = "$projectDir/backbones/BB43b.fasta"
}
else if (params.backbone == "BB22") {
    backbone_file = "$projectDir/backbones/BB22.fasta"
}
else if (params.backbone == "BB25") {
    backbone_file = "$projectDir/backbones/BB25.fasta"
}
else if (params.backbone == "BBCS") {
    backbone_file = "$projectDir/backbones/BBCS.fasta"
}
else if (params.backbone == "BBCR") {
    backbone_file = "$projectDir/backbones/BBCR.fasta"
}
else {
    backbone_file = params.backbone_file
}


// ### Printout for user
log.info """
    ===================================================
    Cyclomics/CyclomicsSeq_Consensus : Consensus generation
    ===================================================
    Inputs:
        input_reads              : $params.read_folder
        read_pattern             : $params.read_pattern
        reference                : $params.reference
        backbone                 : $params.backbone
        backbone_file            : $params.backbone_file
        output folder            : $params.output_dir
        Cmd line                 : $workflow.commandLine
    Method:  
        consensus_method        : $params.consensus_method

    Other:
        tbd
"""

include {
    Cycas
    Cyclotron
    Cygnus
    CygnusPrimed
    TideHunter
} from "./subworkflows"

/*
========================================================================================
    Workflow
========================================================================================
*/
workflow {
    log.info """Cyclomics consensus pipeline started"""

    // Process inputs:
    // add the trailing slash if its missing 
    if (params.read_folder.endsWith("/")){
        read_pattern = "${params.read_folder}${params.read_pattern}"
    }
    else {
        read_pattern = "${params.read_folder}/${params.read_pattern}"
    }
    read_dir_ch = Channel.fromPath( params.read_folder, type: 'dir', checkIfExists: true)
    
    // Create an item where we have the path and the sample ID.
    //  If the path is a directory with fastq's in it directly,
    // The sample ID will be that directory name. eg. fastq_pass/ could be the sample ID.
    // Alternatively, the sample ID could be barcode01/barcode02 etc.
    read_fastq = Channel.fromPath(read_pattern, checkIfExists: true) \
        | map(x -> [x.Parent.simpleName, x.simpleName,x])
    

    // Based on the selected method collect the other inputs and start pipelines.
    if (params.consensus_method == "cycas") {
        log.info """Cycas consensus generation method selected."""
        backbone  = Channel.fromPath(backbone_file, checkIfExists: true)
        reference = Channel.fromPath(params.reference, checkIfExists: true)
        // Cycas(read_dir_ch, backbone, reference)
        log.warn """Please run the default implementation of cyclomicsseq. \n\n exitting....."""
    }
    else if (params.consensus_method == "Cyclotron") {
        log.info """Cyclotron consensus generation method selected."""
        backbone  = Channel.fromPath(backbone_file, checkIfExists: true)
        Cyclotron(read_fastq, backbone)
    }
    else if (params.consensus_method == "Cygnus") {
        log.info """Cygnus consensus generation method selected."""
        Cygnus(read_fastq)
    }
    else if (params.consensus_method == "Cygnus_primed") {
        log.info """Cygnus_primed consensus generation method selected with primer rotation."""
        backbone  = Channel.fromPath(backbone_file, checkIfExists: true)
        CygnusPrimed(read_fastq, backbone, backbone, params.backbone_barcode)
    }
    else if (params.consensus_method == "tidehunter") {
        log.info """TideHunter consensus generation method selected."""
        backbone  = Channel.fromPath(backbone_file, checkIfExists: true)
        TideHunter(read_fastq, backbone)
    }
    else {
        log.warn "Unknown consensus generation method selected"
    }

    
}