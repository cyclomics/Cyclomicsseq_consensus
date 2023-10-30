process RunCygnusConsensus{
    publishDir "${params.output_dir}/consensus/${sample}", mode: 'copy'
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
    publishDir "${params.output_dir}/consensus/${sample}", mode: 'copy'
    input:
        tuple val(sample), val(ID), path(fq), path(primers)
    output:
        tuple val(sample), val(ID), path("*.fastq")
    script:
        """
        python3 $params.rotate_sequence_location $fq $primers ${fq.simpleName}_rotated_consenus.fastq
        """
}

process RunCyclotronConsensus{
    publishDir "${params.output_dir}/consensus/${sample}", mode: 'copy'
    input:
        tuple val(sample), val(ID), path(fq), path(backbone)
    output:
        tuple val(sample), val(ID), path("*.fastq")
    script:
        """
        python3 $params.cyclotron_location $fq $backbone ${fq.simpleName}_consenus.fastq
        """
}

process Cut5P3PAdapter {
    publishDir "${params.output_dir}/trimmed/${sample}", mode: 'copy'

    input:
        tuple val(sample), val(ID), path(fq)
        val(adapter_5p)
        val(adapter_3p)
    output:
        tuple val(ID), path("${fq.simpleName}_cutadapt.fq.gz")
    script:
    """
    cutadapt -a $adapter_3p -g $adapter_5p --revcomp -m 5 -o ${fq.simpleName}_cutadapt.fq.gz $fq
    """
}

process Cut3PAdapter {
    publishDir "${params.output_dir}/trimmed/${sample}", mode: 'copy'

    input:
        tuple val(sample), val(ID), path(fq)
        val(adapter_3p)

    output:
        tuple val(ID), path("${fq.simpleName}_cutadapt.fq.gz")
    script:
    """
    cutadapt -a $adapter_3p --revcomp -m 5 -o ${fq.simpleName}_cutadapt.fq.gz $fq
    """
}

process Cut5PAdapter {
    publishDir "${params.output_dir}/trimmed/${sample}", mode: 'copy'

    input:
        tuple val(sample), val(ID), path(fq)
        val(adapter_3p)
        
    output:
        tuple val(ID), path("${fq.simpleName}_cutadapt.fq.gz")
    script:
    """
    cutadapt -g $adapter_3p --revcomp -m 5 -o ${fq.simpleName}_cutadapt.fq.gz $fq
    """
}

process Cut5P3PFullBackbone {
    publishDir "${params.output_dir}/trimmed/${sample}", mode: 'copy'

    input:
        tuple val(sample), val(ID), path(fq)
        val(sequence)
        
    output:
        tuple val(ID), path("${fq.simpleName}_cutadapt.fq.gz")
    script:
    """
    cutadapt -g $sequence -a $sequence --revcomp -m 5 -o ${fq.simpleName}_cutadapt.fq.gz $fq
    """
}

process DemuxBackboneBarcodes {
    // publishDir "${params.output_dir}/backbone_demuxed/${sample}", mode: 'copy'

    input:
        tuple val(ID), val(sample), path(fq), path(backbone_fa)

    output:
         tuple val(ID), val(sample), path("backbone_demux/*.fastq")

    script:
    """
    mkdir backbone_demux
    detect_backbone_barcode.py $backbone_fa $fq backbone_demux
    for i in \$(ls backbone_demux); do
        mv backbone_demux/\$i backbone_demux/\$(basename \$i).${ID}_${sample}.fastq
    done
    """

}

process RegroupReads {
    publishDir "${params.output_dir}/demux_regroup/${ID}", mode: 'copy'

    input:
        tuple val(ID), val(sample), path(fq)

    output:
         tuple val(ID), val(sample), path("*.fastq")

    script:
    """
    seqkit seq **.fastq > ${ID}.fastq
    """
}