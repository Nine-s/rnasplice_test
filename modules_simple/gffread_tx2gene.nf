process GFFREAD_TX2GENE {
    tag "$gtf"
    label 'ALL'

    //conda "bioconda::gffread=0.12.1"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/gffread:0.12.1--h8b12597_0' :
    //    'biocontainers/gffread:0.12.1--h8b12597_0' }"
    //container 'biocontainers/gffread:0.12.1--h8b12597_0'
    container "zavolab/gffread:0.11.7"
    
    input:
    path gtf

    output:
    path "*.tx2gene.tsv" , emit: tx2gene

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: '--table transcript_id,gene_id,gene_name'
    def prefix = task.ext.prefix ?: "${gtf.baseName}"
    """
    gffread $args $gtf | sort -u 1> ${prefix}.tx2gene.tsv
    """
}
