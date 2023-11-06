process DEXSEQ_COUNT {
    label 'process_medium'

    //conda "bioconda::htseq=2.0.2"
    container "biocontainers/htseq:v0.11.2-1-deb-py3_cv1"

    input:
    tuple val(name), path(bam), val(condition)
    path (gff)
    val alignment_quality                   // val params.alignment_quality

    output:
    tuple val(name), path("*.clean.count.txt"), emit: dexseq_clean_txt

    script:

    def read_type = '-p yes'

    def alignment_quality = "-a ${alignment_quality}"

    def strandedness = ''
    if (params.strand == 'forward') {
        strandedness = '-s yes'
    } else if (params.strand == 'reverse') {
        strandedness = '-s reverse'
    } else if (params.strand == 'unstranded') {
        strandedness = '-s no'
    }

    script:
    """
    ${params.basedir}/bin/dexseq_count.py $gff $read_type -f bam $bam -r pos ${name}.clean.count.txt $alignment_quality $strandedness
    """
}