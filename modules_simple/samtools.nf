process SAMTOOLS {
    label 'samtools'
    publishDir params.outdir

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
    'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(sample_name), path(sam_file)
    
    output:
    tuple val(sample_name), path("${sample_name}.sorted.bam"), emit: bam 
    
    script:
    """
    samtools view -bS ${sam_file} -@ ${params.threads} | samtools sort -o ${sample_name}.sorted.bam -T tmp  -@ ${params.threads} 
    """
    
}