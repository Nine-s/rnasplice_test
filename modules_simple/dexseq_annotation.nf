
process DEXSEQ_ANNOTATION {
    tag "$gtf"
    label 'process_medium'

    container 'quay.io/biocontainers/htseq:2.0.2--py310ha14a713_0'

    input:
    path gtf         // path gtf file
    val aggregation  // val params.aggregation

    output:
    path "*.gff"        , emit: gff

    script:

    def aggregation = aggregation ? '' : '-r no'

    """
    dexseq_prepare_annotation.py $gtf ${prefix}.gff $aggregation
    """
}
