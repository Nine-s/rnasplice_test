
process TRIMGALORE {

    container "dukegcb/trim-galore:0.4.4"

    // Input FastQ files for read 1 and read 2
    input:
        tuple val(name), path(reads), val(condition)
    
    output:
        tuple val(name), path("*val*.f*q*.gz"), val(condition), emit: reads
    script:
    """
    trim_galore --paired ${reads[0]} ${reads[1]}
    """
}