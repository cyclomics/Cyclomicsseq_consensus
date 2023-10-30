
include {
    RunCygnusConsensus
    RunCygnusConsensusUnpublished
    RunCyclotronConsensus
    RotateBySequence
    Cut5P3PAdapter
    Cut3PAdapter
    Cut5PAdapter
    Cut5P3PFullBackbone
    DemuxBackboneBarcodes
    RegroupReads
} from "./processes"

// include {
//     Extract5PrimeFasta
// } from "./nextflow_utils/parse_convert/modules/seqkit"

include {
    ExtractFullFasta
    Extract5PrimeFasta
    Extract3PrimeRemainderFasta
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
        consensus_trimmed

    main:
        length_prof = 15
        log.info """Cyclotron subworkflow pipeline started"""

        pre_consensus_array = read_fastqs.combine(backbone)
        consensus_with_bb = RunCyclotronConsensus(pre_consensus_array)

        // Create value channels out of the output by calling the first function that is a value factory method.
        backbone_seq_3p = Extract3PrimeRemainderFasta(backbone, length_prof).first()

        // backbone_seq_3p.view()
        // TODO: Xadapter flags the adapter in the cutadapt cmd??
        // consensus = Cut3PAdapter(consensus_with_bb, backbone_seq_3p)

        adapter_sequence = ExtractFullFasta(backbone).first()
        consensus_trimmed =  Cut5P3PFullBackbone(consensus_with_bb, adapter_sequence)

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
        adapter
        backbone_demux

    emit:
        consensus
    main:
        log.info """Cygnus Primed subworkflow pipeline started"""
        
        pre_consensus = RunCygnusConsensusUnpublished(read_fastqs)
        // primers.view()
        // pre_consensus.view()
        // pre_consensus.map{ it -> [it[1], it[2].size()]}.view()
        // empty gzipped files are not size 0.
        pre_consensus = pre_consensus.filter{ it -> it[2].size() >1024}
        // pre_consensus.map{ it -> [it[1], it[2].size()]}.view()

        pre_consensus_combined = pre_consensus.combine(primers)
        // pre_consensus_combined.view()
        consensus_rotated = RotateBySequence(pre_consensus_combined)
        // consensus_rotated.view()

        if (backbone_demux == true) {
            log.info """Cygnus Primed: Detecting barcodes in Backbone"""

            demuxed = DemuxBackboneBarcodes(consensus_rotated.combine(adapter))

            // Demux makes files per barcode found, so we need to transform to a tuple ber barcode
            // transform the tuple structure
            demux_flapmapped = demuxed.map {it -> it[2].collect {fq -> [it[0], it[1], fq]}}.flatMap()
            // rename the ID to the barcode, move the old ID to the sample column
            demux_reIDed = demux_flapmapped.map {it -> [it[2].simpleName, it[0], it[2]]}
            // demux_flapmapped.view()
            // demux_reIDed.view()
            per_barcode = demux_reIDed.groupTuple(by:0)
            // per_barcode.view()
            
            per_barcode_merged = RegroupReads(per_barcode)

            trimmable = per_barcode_merged
        }
        else {
            log.info """Cygnus Primed: Keeping folder barcode structure"""

            trimmable = consensus_rotated
        }
        // println("DEV: Trimmable")
        // trimmable.view()

        adapter_sequence = ExtractFullFasta(adapter).first()
        consensus_trimmed = Cut5P3PFullBackbone(trimmable, adapter_sequence)

        log.info """Cygnus Primed: Removing empty files."""
        // consensus_trimmed = consensus_trimmed.filter{ it -> it[2].size() >1024}

        consensus = consensus_trimmed
        // adapter.view()
        // consensus_trimmed.combine(adapter).view()
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
        exit(1)
}
