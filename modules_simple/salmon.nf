process SALMON_QUANT {
    label 'samtools'
    publishDir params.outdir
    container "biocontainers/salmon:v0.12.0ds1-1b1-deb_cv1"
    
    input:
    tuple val(sample_name), path(bam_file), val(condition)
    path(ref_transcripts)
    
    output:
    tuple val(sample_name), path("salmon_quant"), val(condition), emit: sample_bam 
    
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


