process SAMTOOLS {
    label 'ALL'
    publishDir params.outdir

    container "biocontainers/samtools:v1.7.0_cv4"

    input:
    tuple val(sample_name), path(sam_file), val(condition)
    
    output:
    tuple val(sample_name), path("${sample_name}*bam"), val(condition), emit: bam 
    
    script:
    """
    samtools view -bS ${sam_file} -@ ${params.threads} | samtools sort -o ${sample_name}.sorted.bam -T tmp  -@ ${params.threads} 
    """
    
}