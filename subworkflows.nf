
include {
    RunCygnusConsensus
    RunCygnusConsensusUnpublished
    RunCyclotronConsensus
    RotateBySequence
} from "./processes"

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
        log.info """Cyclotron subworkflow pipeline started"""
        // backbone.view()
        pre_consensus_array = read_fastqs.combine(backbone)
        // pre_consensus_array.view()
        consensus = RunCyclotronConsensus(pre_consensus_array)
        // if (params.dev_mode) {
            // backbone.view()
            // pre_consensus_array.first().view()
            // consensus.first().view()
        // }
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
