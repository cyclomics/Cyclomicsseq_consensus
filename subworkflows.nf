

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

        consensus = ""
}

workflow Cygnus {
    take:
        read_fastqs

    emit:
        consensus
    main:
        log.info """Cygnus subworkflow pipeline started"""
        consensus = ""
}

workflow CygnusPrimed {
    take:
        read_fastqs
        primers

    emit:
        consensus
    main:
        log.info """Cygnus Primed subworkflow pipeline started"""
        pre_consensus = Cygnus(read_fastqs)
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
