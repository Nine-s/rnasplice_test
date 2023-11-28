process MULTIQC {
    label 'ALL'
    publishDir params.outdir
    container "staphb/multiqc:1.8"

    input:
    // path(salmon)
    // path(trim_galore)
    // path(star)
    // path(fastqc)
    path(outdir)

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    script:
        """

        # Create the output directory
        mkdir -p multiqc_report


        # Run MultiQC to generate the report
        multiqc -o multiqc_report ${params.outdir}
        """
}
        // echo ${salmon}
        // echo ${trim_galore}
        // echo ${star}
        // echo ${fastqc}