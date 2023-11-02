process STAR_GENOMEGENERATE {
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