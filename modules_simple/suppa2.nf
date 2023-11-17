
process runSUPPA {
    label 'ALL'

    // Input GTF file, BAM files, and event definition file
    input:
    path gtf_file
    file bam_files

    // Output directory for SUPPA results
    output:
    path "suppa_output"

    script:
    """
    ${params.basedir}/bin/suppa.py generateEvents -i \${bam_file} -o events.csv -e \${gtf_file}
    ${params.basedir}/bin/suppa.py psi -i events.csv -o psi.csv
    """
}