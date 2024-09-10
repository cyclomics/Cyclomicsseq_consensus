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
params.read_folder    = ""
params.read_pattern   = "**.{fq,fastq,fq.gz,fastq.gz}"
params.backbone       = ""
params.primer_file    = ""
params.reference      = ""
// Backbone file is used for custom backbones.
params.backbone_file  = ""

params.output_dir = "$HOME/Data/CyclomicsSeq"
params.backbone_barcode = false


// method selection
params.summarize_input  = true
params.summarize_output  = true
params.skip_alignment = false
params.consensus_method = "Cycas"

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

log.debug "parameters obj : ${params}"
log.debug "workflow object: ${workflow}"

// ### Printout for user workflow.revision is the github version tag.as
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
        Version  : $workflow.revision
        Revision : $workflow.commitId
"""

include {
    Cyclotron
    TideHunter
} from "./subworkflows"

// include {
//     Cycas
// } from "./nextflow_utils/consensus/modules/cycas"

include {
    CycasConsensus
    CygnusConsensus
    CygnusAlignedConsensus
    CygnusPrimedConsensus
    Cygnus2AlignedConsensus
} from "./nextflow_utils/consensus/subworkflows"

include {
    PrepareGenome
} from "./nextflow_utils/parse_convert/subworkflows"

include {
    Minimap2Index
    Minimap2Align as Minimap2AlignByID
} from "./nextflow_utils/parse_convert/modules/minimap"

include {
    MergeFasta
} from "./nextflow_utils/parse_convert/modules/seqkit"

include {
    SummerizeReadsStdout as SummarizePerSampleID_in
    SummerizeReadsStdout as SummarizePerSampleID_out
} from "./nextflow_utils/reporting/modules/seqkit"

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

    // Create an item where we have the path and the sample ID and the file ID.
    //  If the path is a directory with fastq's in it directly,
    // The sample ID will be that directory name. eg. fastq_pass/ could be the sample ID.
    // Alternatively, the sample ID could be barcode01/barcode02 etc.
    read_fastq = Channel.fromPath(read_pattern, checkIfExists: true) \
        | map(x -> [x.Parent.simpleName, x.simpleName,x])
    read_fastq.dump(tag: "input-data")

    if (params.summarize_input){
        summary_in = SummarizePerSampleID_in(read_fastq.groupTuple())
        summary_in.subscribe { x ->
            log.info "\nSummary per sample of the input:\n$x"
        }
    }

    // Based on the selected method collect the other inputs and start pipelines.
    if (params.consensus_method == "Cycas") {
        if (params.reference == "" || backbone_file == "") {
            log.error \
            """Please provide reference genome and backbone file for Cycas method.
            reference genome can be provided with --reference and was: '${params.reference}' 
            backbone file can be provided with --backbone_file or --backbone and were: '${params.backbone_file}' or '${params.backbone}'
            """
            // we need some delay to display the error message above (in ms). 
            sleep(200)
            exit 1
        }
        log.info """Cycas consensus generation method selected."""
        backbone  = Channel.fromPath(backbone_file, checkIfExists: true)
        reference = Channel.fromPath(params.reference, checkIfExists: true)
        PrepareGenome(reference, params.reference, backbone)
        // .collect() to turn into repeating value channel.
        reference_mmi = PrepareGenome.out.mmi_combi.collect()
        
        CycasConsensus(read_fastq, reference_mmi)
        consensus = CycasConsensus.out
        // Drop the metadata jsons for now
        consensus = consensus.map{ it -> it.take(3)}
    }
    else if (params.consensus_method == "Cyclotron") {
        log.info """Cyclotron consensus generation method selected."""
        backbone  = Channel.fromPath(backbone_file, checkIfExists: true)
        consensus = Cyclotron(read_fastq, backbone)

    }

    else if (params.consensus_method == "Cygnus") {
        log.info """Cygnus consensus generation method selected."""
        consensus = Cygnus(read_fastq)

    }
    
    else if (params.consensus_method == "Cygnus_primed") {
        log.info """Cygnus_primed consensus generation method selected with primer rotation."""
        if (params.primer_file == ""){
            log.error \
            """Please provide a primer file for Cygnus_primed method.
            primer file can be provided with --primer_file and was: '${params.primer_file}'
            """
            // we need some delay to display the error message above (in ms). 
            sleep(200)
            exit 1
        }
        primer = Channel.fromPath(params.primer_file, checkIfExists: true)
        // We need a value channel to repeat the usage.
        primer = primer.collect()
        CygnusPrimedConsensus(read_fastq, primer)
        consensus = CygnusPrimedConsensus.out
        consensus.dump()
        // consensus.view()
    }

    else if (params.consensus_method == "Cygnus_aligned") {
        log.info """Cygnus_aligned consensus generation method selected."""
        log.info """We will align against the provided primer."""
        if (params.reference == "" || backbone_file == "") {
            log.error \
            """Please provide reference genome and backbone file for Cygnus)aligned method.
            reference genome can be provided with --reference and was: '${params.reference}' 
            backbone file can be provided with --backbone_file or --backbone and were: '${params.backbone_file}' or '${params.backbone}'
            """
            // we need some delay to display the error message above (in ms). 
            sleep(200)
            exit 1
        }
        backbone  = Channel.fromPath(backbone_file, checkIfExists: true)
        reference = Channel.fromPath(params.reference, checkIfExists: true)
        PrepareGenome(reference, params.reference, backbone)
        // .collect() to turn into repeating value channel.
        reference_mmi = PrepareGenome.out.mmi_combi.collect()
        consensus = CygnusAlignedConsensus(read_fastq, reference_mmi)
    }
    else if (params.consensus_method == "Cygnus2_aligned") {
        log.info """Cygnus2_aligned consensus generation method selected."""
        log.info """We will align against the provided primer."""
        if (params.reference == "" || backbone_file == "") {
            log.error \
            """Please provide reference genome and backbone file for Cygnus)aligned method.
            reference genome can be provided with --reference and was: '${params.reference}' 
            backbone file can be provided with --backbone_file or --backbone and were: '${params.backbone_file}' or '${params.backbone}'
            """
            // we need some delay to display the error message above (in ms). 
            sleep(200)
            exit 1
        }
        backbone  = Channel.fromPath(backbone_file, checkIfExists: true)
        reference = Channel.fromPath(params.reference, checkIfExists: true)
        PrepareGenome(reference, params.reference, backbone)
        // .collect() to turn into repeating value channel.
        reference_mmi = PrepareGenome.out.mmi_combi.collect()
        consensus = Cygnus2AlignedConsensus(read_fastq, reference_mmi)
    }

    else if (params.consensus_method == "Tidehunter") {
        log.info """TideHunter consensus generation method selected."""
        backbone  = Channel.fromPath(backbone_file, checkIfExists: true)
        TideHunter(read_fastq, backbone)
    }

    else {
        log.warn "Unknown consensus generation method selected"
        sleep(200)
        exit 1
    }

    // Publish the consensus fastqs grouped in folders by their sample id
    consensus.map { it ->
        it[2].copyTo("${params.output_dir}/results/${it[0]}/${it[2].name}")
        }
    
    if (params.summarize_output){
        summary_out = SummarizePerSampleID_out(consensus.groupTuple())
        summary_out.subscribe { x ->
            log.info "\nSummary per sample of the output:\n$x"
        }
    }

    if (params.reference != "" && !params.skip_alignment){
        log.warn "Aligning to provided reference."
        sleep(200)
        
        // If we have a reference, we then might as well align it all per sample ID.
        reference = Channel.fromPath(params.reference)
        Minimap2Index(reference)
        reference_mmi_final = Minimap2Index.out
        consensus_by_id = consensus.groupTuple().map{ it -> [it[0], it[0], it[2]]}
        consensus_by_id.dump(tag: 'consensus-pre-alignment')
        Minimap2AlignByID(consensus_by_id, reference_mmi_final)
        consensus_aligned = Minimap2AlignByID.out
        // Publish the aligned consensus, but as a single file
        consensus_aligned.map { it ->
        it[2].copyTo("${params.output_dir}/results_aligned/${it[2].name}")    // Copy the file to the target directory
        }
    }
}


workflow.onComplete{
    log.info ("\nDone. The results are available in following folder --> $params.output_dir\n")
}
