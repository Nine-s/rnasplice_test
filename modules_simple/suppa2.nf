
process runSUPPA {

    // Input GTF file, BAM files, and event definition file
    input:
    path gtf_file
    file bam_files

    // Output directory for SUPPA results
    output:
    path "suppa_output"

    script:
    """
    suppa.py generateEvents -i \${bam_file} -o events.csv -e \${gtf_file}
    suppa.py psi -i events.csv -o psi.csv
    """
}