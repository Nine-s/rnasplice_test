process MULTIQC {
    label 'ALL'
    publishDir params.outdir
    container "staphb/multiqc:1.8"

    input:
    path(salmon_transcripts)
    path(salmon_json_info)
    tuple val(name_trim), path("*val*.f*q*.gz"), val(condition)
    path(trim_log)                                             
    tuple val(sample_name), path("${sample_name}*.sam"), val(condition)
    path(star_log_final)   
    path(star_log_out)   
    path(fastqc)

    output:
    path "*multiqc_report"

    script:
    """
    mkdir multiqc_report
    multiqc .
    """
}
    // path "*multiqc_report.html", emit: report
    // path "*_data"              , emit: data
    // path "*_plots"             , optional:true, emit: plots

    // mkdir -p multiqc_report
    // multiqc -o multiqc_report .


    //
        // echo ${salmon}
        // echo ${trim_galore}
        // echo ${star}
        // echo ${fastqc}