process RMATS_PREP {
    tag "$cond1-$cond2"
    label 'process_high'

    //container 'quay.io/biocontainers/mulled-v2-8ea76ff0a6a4c7e5c818fd4281abf918f92eeeae:121e48ab4817ec619c157a346458efca1ccf3c0a-0'
    //container "zavolab/rmats:4.0.2" 
    //container "zavolab/rmats2sashimiplot:2.0.3"
    container "xinglab/rmats:v4.2.0"

    input:
    path gtf                                     // /path/to/genome.gtf
    tuple val(contrast), val(cond1), val(cond2), val(meta1), val(meta2), path(bam1), path(bam2), path(bam1_text), path(bam2_text)
    val rmats_read_len                           // val params.rmats_read_len
    val rmats_splice_diff_cutoff                 // val params.rmats_splice_diff_cutoff
    val rmats_novel_splice_site                  // val params.rmats_novel_splice_site
    val rmats_min_intron_len                     // val params.rmats_min_intron_len
    val rmats_max_exon_len                       // val params.rmats_max_exon_len

    output:
    tuple val(contrast), path("$cond1-$cond2/rmats_temp/*") , emit: rmats_temp
    path "$cond1-$cond2/rmats_prep.log"                     , emit: log

    script:
    def prefix = task.ext.prefix ?: "$cond1-$cond2"

    // Default strandedness to fr-unstranded - also if user supplies "unstranded"
    def strandedness = 'fr-unstranded'

    // Change strandedness based on user samplesheet input
    if (params.strand == 'forward') {
        strandedness  = 'fr-secondstrand'
    } else if (params.strand == 'reverse') {
        strandedness  = 'fr-firststrand'
    }

    def read_type = 'paired'

    // Whether user wants to run with novel splice sites flag
    def novel_splice_sites = rmats_novel_splice_site ? '--novelSS' : ''

    // Additional args for when running with --novelSS flag
    // User defined else defauls to 50, 500
    def min_intron_len = ''
    def max_exon_len   = ''
    if (rmats_novel_splice_site) {
        min_intron_len = rmats_min_intron_len ? "--mil ${rmats_min_intron_len}" : '--mil 50'
        max_exon_len   = rmats_max_exon_len ? "--mel ${rmats_max_exon_len}" : '--mel 500'
    }

    """
    mkdir -p $prefix/rmats_temp

    mkdir -p $prefix/rmats_prep

    python /rmats/rmats.py \\
        --gtf $gtf \\
        --b1 $bam1_text \\
        --b2 $bam2_text \\
        --od $prefix/rmats_prep \\
        --tmp $prefix/rmats_temp \\
        -t $read_type \\
        --libType $strandedness \\
        --readLength $rmats_read_len \\
        --variable-read-length \\
        --nthread 1 \\
        --cstat $rmats_splice_diff_cutoff \\
        --task prep \\
        $novel_splice_sites \\
        $min_intron_len \\
        $max_exon_len \\
        --allow-clipping \\
        1> $prefix/rmats_prep.log

    """

}
//        --tstat 1 \\