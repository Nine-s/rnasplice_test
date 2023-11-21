//include { PREPARE_GENOME                                       } from '../subworkflows/local/prepare_genome'
include { CAT_FASTQ } from './modules_simple/cat.nf'
include { FASTQC } from './modules_simple/fastqc.nf'
include { TRIMGALORE } from './modules_simple/trimgalore.nf'
include { CUSTOM_GETCHROMSIZES              } from './modules_simple/getchromsizes.nf'
include { GFFREAD_TX2GENE  } from './modules_simple/gffread_tx2gene.nf'
include { TXIMPORT         } from './modules_simple/tximport.nf'
include { UNTAR            } from './modules_simple/untar.nf'

include { STAR_GENOMEGENERATE } from './modules_simple/star_genome_generate.nf'
include { STAR_ALIGN } from './modules_simple/star_align.nf'
include { SAMTOOLS } from './modules_simple/samtools.nf'

include { SALMON_GENOMEGENERATE } from './modules_simple/salmon_genome_generate.nf'
include { SALMON_QUANT  } from './modules_simple/salmon.nf'

include { BEDTOOLS_GENOMECOV } from './modules_simple/bedtools_genomecov'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_TO_BIGWIG_FORWARD } from '../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig/main'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_TO_BIGWIG_REVERSE } from '../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig/main'

include { MULTIQC } from './modules_simple/multiqc.nf'

include { DEXSEQ_ANNOTATION } from './modules_simple/dexseq_annotation.nf'
include { DRIMSEQ_FILTER  } from './modules_simple/drimseq_filter.nf'
include { DEXSEQ_DTU      } from './modules_simple/dexseq_dtu.nf'
include { DEXSEQ_COUNT } from './modules_simple/dexseq_count.nf'
include { DEXSEQ_EXON } from './modules_simple/dexseq_exon.nf'
//include { TX2GENE_TXIMPORT as TX2GENE_TXIMPORT_STAR_SALMON     } from 
//include { DIFFSPLICE } from './modules_simple/modules_simple/suppa_splicing_analysis.nf'
//'../subworkflows/local/tx2gene_tximport'

workflow{

read_pairs_ch = Channel
    .fromPath( params.csv_input )
    .splitCsv(header: true, sep: ',')
    .map {row -> tuple(row.sample, [row.fastq_1, row.fastq_2], row.condition)}
    .view()

//CAT_FASTQ(read_pairs_ch)

FASTQC( read_pairs_ch )

TRIMGALORE( read_pairs_ch )



SALMON_GENOMEGENERATE ( params.genome, params.transcripts_fasta )
SALMON_QUANT(TRIMGALORE.out.reads, SALMON_GENOMEGENERATE.out.index)

STAR_GENOMEGENERATE(params.genome, params.annotation_gtf)

STAR_ALIGN(TRIMGALORE.out.reads, STAR_GENOMEGENERATE.out, params.annotation_gtf )

SAMTOOLS( STAR_ALIGN.out.sam )

//
//// STEP 5: Create bigWig coverage files 
//
CUSTOM_GETCHROMSIZES(params.genome)
BEDTOOLS_GENOMECOV(SAMTOOLS.out.bam)
BEDCLIP(BEDTOOLS_GENOMECOV.out.bedgraph_forward, CUSTOM_GETCHROMSIZES.out.sizes)
BEDGRAPHTOBIGWIG(BEDTOOLS_GENOMECOV.out.bedgraph_forward, CUSTOM_GETCHROMSIZES.out.sizes)
//TODO add reverse
//
//9. Summarize QC (MultiQC)
//

//MULTIQC()

//
// 11. Differential exon usage with DEXSeq or edgeR
//

DEXSEQ_ANNOTATION(params.annotation_gtf)//, params.aggregation)

//get gff for this part
DEXSEQ_COUNT (
    SAMTOOLS.out.bam,
    DEXSEQ_ANNOTATION.out.gff,
    params.alignment_quality
)

dexseq_clean_counts = MERGE_RESULTS_DEXSEQ(DEXSEQ_COUNT.out.dexseq_clean_txt.collect())

DEXSEQ_EXON (
    dexseq_clean_counts,
    DEXSEQ_ANNOTATION.out.gff,
    params.csv_input,
    params.csv_contrastsheet,
    params.n_dexseq_plot
)

// 12. Differential transcript usage

// https://github.com/nf-core/rnasplice/blob/dev/subworkflows/local/drimseq_dexseq_dtu.nf

GFFREAD_TX2GENE ( params.annotation_gtf )

salmon_results = MERGE_RESULTS_SALMON(SALMON_QUANT.out.transcripts.collect())

TXIMPORT ( salmon_results, GFFREAD_TX2GENE.out.tx2gene )

DRIMSEQ_FILTER( TXIMPORT.out.txi_dtu, TXIMPORT.out.tximport_tx2gene, params.csv_input, params.min_samps_gene_expr, params.min_samps_feature_expr, params.min_samps_feature_prop, params.min_feature_expr, params.min_feature_prop, params.min_gene_expr )

DEXSEQ_DTU(DRIMSEQ_FILTER.out.drimseq_samples_tsv, DRIMSEQ_FILTER.out.drimseq_counts_tsv, params.csv_contrastsheet, params.n_dexseq_plot)




//
// 12. Event-based splicing analysis:
//

// RMATS (
//                 ch_contrastsheet,
//                 ch_genome_bam_conditions,
//                 PREPARE_GENOME.out.gtf,
//                 is_single_condition,
//                 params.rmats_read_len,
//                 params.rmats_splice_diff_cutoff,
//                 params.rmats_novel_splice_site,
//                 params.rmats_min_intron_len,
//                 params.rmats_max_exon_len,
//                 params.rmats_paired_stats
//             )

//             ch_versions = ch_versions.mix(RMATS.out.versions)

//  ch_suppa_tpm = params.suppa_tpm ? PREPARE_GENOME.out.suppa_tpm : TX2GENE_TXIMPORT_STAR_SALMON.out.suppa_tpm

//             // Run SUPPA
//             SUPPA_STAR_SALMON (
//                 PREPARE_GENOME.out.gtf,
//                 ch_suppa_tpm,
//                 ch_samplesheet,
//                 ch_contrastsheet,
//                 params.suppa_per_local_event,
//                 params.generateevents_boundary,
//                 params.generateevents_threshold,
//                 params.generateevents_exon_length,
//                 params.generateevents_event_type,
//                 params.generateevents_pool_genes,
//                 params.psiperevent_total_filter,
//                 params.diffsplice_local_event,
//                 params.diffsplice_method,
//                 params.diffsplice_area,
//                 params.diffsplice_lower_bound,
//                 params.diffsplice_alpha,
//                 params.diffsplice_tpm_threshold,
//                 params.diffsplice_nan_threshold,
//                 params.diffsplice_gene_correction,
//                 params.diffsplice_paired,
//                 params.diffsplice_median,
//                 params.clusterevents_local_event,
//                 params.clusterevents_dpsithreshold,
//                 params.clusterevents_eps,
//                 params.clusterevents_metric,
//                 params.clusterevents_min_pts,
//                 params.clusterevents_method,
//                 params.clusterevents_sigthreshold ?: false,
//                 params.clusterevents_separation ?: false,
//                 params.suppa_per_isoform
//             )

}

process MERGE_RESULTS_SALMON {
    publishDir params.outdir
    container "zavolab/salmon:1.1.0"

    input:
    path out_folders
    
    output:
    path("salmon")
    
    script:
    """
    mkdir salmon
    mv  ${out_folders} salmon
    """
}


process MERGE_RESULTS_DEXSEQ {
    publishDir params.outdir
    container "zavolab/salmon:1.1.0"

    input:
    path(out_files)
    
    output:
    path("dexseq_clean_counts")
    
    script:
    """
    mkdir dexseq_clean_counts
    mv  ${out_files} dexseq_clean_counts
    """
}