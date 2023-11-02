
process TRIMGALORE {

    container "dukegcb/trim-galore:0.4.4"

    // Input FastQ files for read 1 and read 2
    input:
        tuple val(name), path(reads), val(condition)
    
    output:
        tuple val(name), path("${name}*val*.f*q"), val(condition) emit: reads
    script:
    """
    trim_galore --paired --output_dir . ${reads[0]} ${reads[1]}
    """
}