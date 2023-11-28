process MULTIQC {
    label 'ALL'
    publishDir params.outdir
    container "staphb/multiqc:1.8"

    input:
    path(salmon_transcripts)
    path(salmon_json_info)
    tuple val(name_trim), path("*val*.f*q*.gz"), val(condition)
    path(trim_log)                        
    path(trim_html)                           
    path(trim_zip)                        
    tuple val(sample_name), path("${sample_name}*.sam"), val(condition)
    path(star_log_final)   
    path(star_log_out)   
    path(fastqc)
    //path(outdir)

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    script:
    //def parent directory for each file

    """
    echo ${salmon}
    echo ${trim_galore}
    echo ${star}
    echo ${fastqc}
    # cp every files in the folders to current location
    mkdir -p multiqc_report
    multiqc -o multiqc_report .
    """
}
        // echo ${salmon}
        // echo ${trim_galore}
        // echo ${star}
        // echo ${fastqc}