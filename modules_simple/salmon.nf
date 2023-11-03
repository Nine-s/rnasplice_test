process SALMON_QUANT {
    label 'samtools'
    publishDir params.outdir
    container "biocontainers/salmon:v0.12.0ds1-1b1-deb_cv1"

    input:
    tuple val(sample_name), path(reads), val(condition)
    path(index)
    
    output:
    //tuple val(sample_name), path("transcripts_quant"), val(condition), emit: quantification 
    path("transcripts_quant"), emit: transcripts

    script:
    """
    if [[ (${params.strand} == "firststrand") ]]; then 
	    salmon quant -i ${index} -l ISR -1 ${reads[0]} -2 ${reads[1]} --validateMappings -o transcripts_quant -p ${params.threads} 
    elif [[ (${params.strand} == "secondstrand") ]]; then 
        salmon quant -i ${index} -l ISF -1 ${reads[0]} -2 ${reads[1]} --validateMappings -o transcripts_quant -p ${params.threads}
	elif [[ ${params.strand} == "unstranded" ]]; then
		salmon quant -i ${index} -l IU -1 ${reads[0]} -2 ${reads[1]} --validateMappings -o transcripts_quant -p ${params.threads}
	else  
		echo ${params.strand} > error_strandness.txt
		echo "strandness cannot be determined" >> error_strandness.txt
	fi
    """
    
}


