process CUSTOM_GETCHROMSIZES {
    label 'ALL'

    container "quay.io/biocontainers/samtools:1.16.1--h6899075_1"

    input:
    path(fasta)

    output:
    tuple val("genome"), path ("*.sizes"), emit: sizes
    tuple val("genome"), path ("*.fai")  , emit: fai
    tuple val("genome"), path ("*.gzi")  , emit: gzi, optional: true

    script:
    """
    samtools faidx ${fasta}
    cut -f 1,2 ${fasta}.fai > ${fasta}.sizes
    """

    stub:
    """
    touch ${fasta}.fai
    touch ${fasta}.sizes
    """
}