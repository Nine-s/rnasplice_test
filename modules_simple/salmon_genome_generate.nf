process SALMON_GENOMEGENERATE {
    label 'ALL'
    
    container "zavolab/salmon:1.1.0"
    //container "biocontainers/salmon:v0.12.0ds1-1b1-deb_cv1"
    //container "biocontainers/salmon:v0.7.2ds1-2b1-deb_cv1"

    input:
    path genome_fasta
    path transcript_fasta

    output:
    path "salmon"      , emit: index

    script:
    def get_decoy_ids = "grep '^>' $genome_fasta | cut -d ' ' -f 1 | cut -d \$'\\t' -f 1 > decoys.txt"
    def gentrome      = "gentrome.fa"
    // if (genome_fasta.endsWith('.gz')) {
    //     get_decoy_ids = "grep '^>' <(gunzip -c $genome_fasta) | cut -d ' ' -f 1 | cut -d \$'\\t' -f 1 > decoys.txt"
    //     gentrome      = "gentrome.fa.gz"
    // }
    """
    $get_decoy_ids
    sed -i.bak -e 's/>//g' decoys.txt
    cat $transcript_fasta $genome_fasta > $gentrome

    salmon index -t ${gentrome} -d decoys.txt -p ${params.threads} -i salmon
    """
//--decoys decoys.txt
}