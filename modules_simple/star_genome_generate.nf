process STAR_GENOMEGENERATE {
    container "mgibio/star:2.7.0f"
    label 'star'
    publishDir params.outdir
    
    input:
    path(reference)
    path(annotation)

    output:
    path("star/*")

    script:
    """
    mkdir star
    STAR \\
            --runMode genomeGenerate \\
            --genomeDir star/ \\
            --genomeFastaFiles ${reference} \\
            --sjdbGTFfile ${annotation} \\
            --runThreadN ${params.threads} \\
	
    """
}