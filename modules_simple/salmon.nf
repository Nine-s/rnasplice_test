process SALMON_QUANT {
    label 'samtools'
    publishDir params.outdir
    
    input:
    tuple val(sample_name), path(bam_file)
    path(ref_transcripts)
    
    output:
    tuple val(sample_name), path("${sample_name}.sorted.bam"), emit: sample_bam 
    
    script:
    """
    if [[ (params.strand == "firststrand") ]]; then 
	    salmon quant -t ${ref_transcripts} -l ISR -a ${bam_file} -o salmon_quant -p ${params.threads} 
    elif [[ (params.strand == "secondstrand") ]]; then 
        salmon quant -t ${ref_transcripts} -l ISF -a ${bam_file} -o salmon_quant -p ${params.threads}
	elif [[ params.strand == "unstranded" ]]; then
		salmon quant -t ${ref_transcripts} -l IU -a ${bam_file} -o salmon_quant -p ${params.threads}
	else  
		echo params.strand > error_strandness.txt
		echo "strandness cannot be determined" >> error_strandness.txt
	fi
    """
    
}


