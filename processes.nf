process RunCygnusConsensus{
    publishDir "${params.output_dir}/${sample}", mode: 'move'
    input:
        tuple val(sample), val(ID), path(fq)
    output:
        tuple val(sample), val(ID), path("*.fastq.gz")
    script:
        """
        python3 $params.cygnus_location $fq ${fq.simpleName}_consenus.fastq.gz
        """
}
process RunCygnusConsensusUnpublished{
    input:
        tuple val(sample), val(ID), path(fq)
    output:
        tuple val(sample), val(ID), path("*.fastq.gz")
    script:
        """
        python3 $params.cygnus_location $fq ${fq.simpleName}_consenus.fastq.gz
        """
}

process RotateBySequence{
    publishDir "${params.output_dir}/${sample}", mode: 'move'
    input:
        tuple val(sample), val(ID), path(fq), path(primers)
    output:
        tuple val(sample), val(ID), path("*.fastq.gz")
    script:
        """
        python3 $params.rotate_sequence_location $fq $primers ${fq.simpleName}_rotated_consenus.fastq.gz
        """
}

process RunCyclotronConsensus{
    publishDir "${params.output_dir}/${sample}", mode: 'move'
    input:
        tuple val(sample), val(ID), path(fq), path(backbone)
    output:
        tuple val(sample), val(ID), path("*.fastq.gz")
    script:
        """
        echo hello world
        python3 $params.cyclotron_location $fq $backbone ${fq.simpleName}_consenus.fastq.gz
        """
}

