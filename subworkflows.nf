
include {
    RunCygnusConsensus
    RunCygnusConsensusUnpublished
    RunCyclotronConsensus
    RotateBySequence
} from "./processes"

// include {
//     Extract5PrimeFasta
// } from "./nextflow_utils/parse_convert/modules/seqkit"

include {
    Extract5PrimeFasta
    Extract3PrimeRemainderFasta
    Cut5P3PAdapter
 } from "./nextflow_utils/sequence_analysis/modules/cutadapt" 


workflow Cycas {
    take:
        read_fastqs
        backbone
        reference

    emit:
        consensus
    main:
        log.info """Cycas subworkflow pipeline started"""

        consensus = ""
}

workflow Cyclotron {
    take:
        read_fastqs
        backbone

    emit:
        consensus
    main:
        length_prof = 15
        log.info """Cyclotron subworkflow pipeline started"""

        pre_consensus_array = read_fastqs.combine(backbone)
        consensus_with_bb = RunCyclotronConsensus(pre_consensus_array)

        backbone_seq_5p = Extract5PrimeFasta(backbone, length_prof)
        backbone_seq_3p = Extract3PrimeRemainderFasta(backbone, length_prof)

        consensus = Cut5P3PAdapter(consensus_with_bb, backbone_seq_5p, backbone_seq_3p)
}


workflow Cygnus {
    take:
        read_fastqs
    emit:
        consensus
    main:
        log.info """Cygnus subworkflow pipeline started"""
        consensus = ""
        RunCygnusConsensus(read_fastqs)
}

workflow CygnusPrimed {
    take:
        read_fastqs
        primers

    emit:
        consensus
    main:
        log.info """Cygnus Primed subworkflow pipeline started"""
        
        pre_consensus = RunCygnusConsensusUnpublished(read_fastqs)
        // primers.view()
        // pre_consensus.view()
        pre_consensus_combined = pre_consensus.combine(primers)
        // pre_consensus_combined.view()
        consensus = RotateBySequence(pre_consensus_combined)
        consensus = ""
}

workflow TideHunter {
    take:
        read_fastqs
        backbone

    emit:
        consensus
    main:
        log.info """Tide subworkflow pipeline started"""

        consensus = ""
}
