process DEXSEQ_DTU {
    label 'process_high'

    //container "zavolab/r_dge_dtu:3.5.1"
    container "biocontainers/htseq:v0.11.2-1-deb-py3_cv1"

    input:
    path drimseq_sample_data
    path drimseq_d_counts
    path drimseq_contrast_data
    val ntop

    output:
    path "DEXSeqDataSet.*.rds"  , emit: dexseq_exon_dataset_rds
    path "DEXSeqResults.*.rds"  , emit: dexseq_exon_results_rds
    path "DEXSeqResults.*.tsv"  , emit: dexseq_exon_results_tsv
    path "perGeneQValue.*.rds"  , emit: dexseq_gene_results_rds
    path "perGeneQValue.*.tsv"  , emit: dexseq_gene_results_tsv

    script:
    """
    ${params.basedir}/bin/run_dexseq_dtu.R $drimseq_sample_data $drimseq_contrast_data $drimseq_d_counts $ntop
    """

}