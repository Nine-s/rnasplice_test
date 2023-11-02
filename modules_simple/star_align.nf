
process STAR_ALIGN {
    label 'star'
    publishDir params.outdir

    input:
    tuple val(sample_name), path(reads)
    path(index)
    path(annotation)

    output:
    tuple val(sample_name), path("${sample_name}*.sam"), emit: sam 

    script:
    """
    if [[ (${params.strand} == "firststrand") || (${params.strand} == "secondstrand") ]]; then
		STAR \\
            --genomeDir . \\
            --readFilesIn ${reads[0]} ${reads[1]}  \\
            --runThreadN ${params.threads} \\
            --outFileNamePrefix ${sample_name}. \\
            --sjdbGTFfile ${annotation} \\
			--alignSoftClipAtReferenceEnds No \\
			--outFilterIntronMotifs RemoveNoncanonical \\
			--outSAMattrIHstart 0
	elif [[ ${params.strand} == "unstranded" ]]; then
		STAR \\
            --genomeDir . \\
            --readFilesIn ${reads[0]} ${reads[1]}  \\
			--alignSoftClipAtReferenceEnds No \\
			--outSAMstrandField intronMotif \\
			--outFilterIntronMotifs RemoveNoncanonical \\
        	--runThreadN ${params.threads} \\
        	--outFileNamePrefix ${sample_name}. \\
        	--sjdbGTFfile ${annotation} \\
			--outSAMattrIHstart 0
	else  
		echo ${params.strand} > error_strandness.txt
		echo "strandness cannot be determined" >> error_strandness.txt
	fi
   """
}
