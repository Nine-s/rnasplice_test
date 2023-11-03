process SAMTOOLS {
    label 'samtools'
    publishDir params.outdir

    container "biocontainers/samtools:v1.7.0_cv4"

    input:
    tuple val(sample_name), path(sam_file)
    
    output:
    tuple val(sample_name), path("${sample_name}.sorted.bam"), emit: bam 
    
    script:
    """
    samtools view -bS ${sam_file} -@ ${params.threads} | samtools sort -o ${sample_name}.sorted.bam -T tmp  -@ ${params.threads} 
    """
    
}