process CUSTOM_GETCHROMSIZES {
    label 'ALL'

    container "biocontainers/samtools:1.16.1--h6899075_1"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path ("*.sizes"), emit: sizes
    tuple val(meta), path ("*.fai")  , emit: fai
    tuple val(meta), path ("*.gzi")  , emit: gzi, optional: true


    script:
    """
    samtools faidx $fasta
    cut -f 1,2 ${fasta}.fai > ${fasta}.sizes
    """

    stub:
    """
    touch ${fasta}.fai
    touch ${fasta}.sizes
    """
}