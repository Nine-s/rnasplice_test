process DEXSEQ_EXON {
    label 'process_high'

    container "zavolab/r_dge_dtu:3.5.1"

    input:
    path ("dexseq_clean_counts/*")    // path: dexseq_clean_counts
    path gff                          // path: dexseq_gff
    path samplesheet                  // path: samplesheet
    path contrastsheet                // path: contrastsheet
    val ntop                          // val: n_dexseq_plot

    output:
    path "DEXSeqDataSet.*.rds"  , emit: dexseq_exon_dataset_rds
    path "DEXSeqResults.*.rds"  , emit: dexseq_exon_results_rds
    path "perGeneQValue.*.rds"  , emit: dexseq_gene_results_rds
    path "DEXSeqResults.*.csv"  , emit: dexseq_exon_results_csv
    path "perGeneQValue.*.csv"  , emit: dexseq_gene_results_csv
    path "plotDEXSeq.*.pdf"     , emit: dexseq_plot_results_pdf
    path "versions.yml"         , emit: versions

    script:
    //${params.baseDir}
    """
    ${params.basedir}/bin/run_dexseq_exon.R dexseq_clean_counts $gff $samplesheet $contrastsheet $ntop

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-dexseq:  \$(Rscript -e "library(DEXSeq); cat(as.character(packageVersion('DEXSeq')))")
    END_VERSIONS
    """
}
