//include { GFFREAD_TX2GENE.out.tx2gene                               } from '../subworkflows/local/GFFREAD_TX2GENE.out.tx2genee { CAT_FASTQ } from './modules_simple/cat.nf'
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

include { BEDTOOLS_GENOMECOV } from './modules_simple/bedtoolsgenomecov.nf'
include { BEDGRAPHTOBIGWIG as BEDGRAPH_TO_BIGWIG_FORWARD } from './modules_simple/bedgraphtobigwig.nf'
include { BEDGRAPHTOBIGWIG as BEDGRAPH_TO_BIGWIG_REVERSE } from './modules_simple/bedgraphtobigwig.nf'
include { BEDCLIP as BEDCLIP_FORWARD } from './modules_simple/bedclip.nf'
include { BEDCLIP as BEDCLIP_REVERSE } from './modules_simple/bedclip.nf'

include { MULTIQC } from './modules_simple/multiqc.nf'

include { DEXSEQ_ANNOTATION } from './modules_simple/dexseq_annotation.nf'
include { DRIMSEQ_FILTER  } from './modules_simple/drimseq_filter.nf'
include { DEXSEQ_DTU      } from './modules_simple/dexseq_dtu.nf'
include { DEXSEQ_COUNT } from './modules_simple/dexseq_count.nf'
include { DEXSEQ_EXON } from './modules_simple/dexseq_exon.nf'

include { CREATE_BAMLIST         } from './modules_simple/create_bamlist.nf'
include { RMATS_PREP             } from './modules_simple/rmats_prep.nf'
include { RMATS_POST             } from './modules_simple/rmats_post.nf'

// include { SPLIT_FILES as SPLIT_FILES_IOE } from './modules_simple/suppa_split_files.nf' 
// include { SPLIT_FILES as SPLIT_FILES_TPM } from '../../modules/local/suppa_split_files.nf'
// include { SPLIT_FILES as SPLIT_FILES_IOI } from '../../modules/local/suppa_split_files.nf'
// include { PSIPEREVENT   } from '../../modules/local/suppa_psiperevent.nf'
// include { PSIPERISOFORM } from '../../modules/local/suppa_psiperisoform.nf'
// include { DIFFSPLICE as DIFFSPLICE_IOE } from '../../modules/local/suppa_diffsplice.nf'
// include { DIFFSPLICE as DIFFSPLICE_IOI } from '../../modules/local/suppa_diffsplice.nf'
// include { GENERATE_EVENTS as GENERATE_EVENTS_IOE } from '../../modules/local/suppa_generateevents.nf'
// include { GENERATE_EVENTS as GENERATE_EVENTS_IOI } from '../../modules/local/suppa_generateevents.nf'
// include { CLUSTERGROUPS as CLUSTERGROUPS_IOE } from '../../modules/local/suppa_clustergroups.nf'
// include { CLUSTEREVENTS as CLUSTEREVENTS_IOE } from '../../modules/local/suppa_clusterevents.nf'

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

BEDCLIP_FORWARD(BEDTOOLS_GENOMECOV.out.bedgraph_forward, CUSTOM_GETCHROMSIZES.out.sizes)
BEDGRAPH_TO_BIGWIG_FORWARD(BEDCLIP_FORWARD.out.bedgraph, CUSTOM_GETCHROMSIZES.out.sizes)

BEDCLIP_REVERSE(BEDTOOLS_GENOMECOV.out.bedgraph_reverse, CUSTOM_GETCHROMSIZES.out.sizes)
BEDGRAPH_TO_BIGWIG_REVERSE(BEDCLIP_REVERSE.out.bedgraph, CUSTOM_GETCHROMSIZES.out.sizes)

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
// 13. Event-based splicing analysis:
//

// ch_genome_bam_conditions = SAMTOOLS.out.bam.map { meta, bam -> [meta.condition, meta, bam] }.groupTuple(by:0)


ch_genome_bam_conditions = SAMTOOLS.out.bam.map { name, bam, condition -> [condition, name, bam] }.groupTuple(by:0)


    ch_contrasts = 
    Channel.fromPath(params.csv_contrastsheet).splitCsv(header:true)

    ch_contrasts = ch_contrasts
        .map { it -> [it['treatment'], it] }
        .combine ( ch_genome_bam_conditions, by: 0 )
        .map { it -> it[1] + ['meta1': it[2], 'bam1': it[3]] }

    ch_contrasts = ch_contrasts
        .map { it -> [it['control'], it] }
        .combine ( ch_genome_bam_conditions, by: 0 )
        .map { it -> it[1] + ['meta2': it[2], 'bam2': it[3]] }


        //
        // Create input channel
        //

        ch_bam = ch_contrasts.map { [ it.contrast, it.treatment, it.control, it.bam1, it.bam2 ] }

        //
        // Create input bam list file
        //

        CREATE_BAMLIST (
            ch_bam
        )

        ch_bamlist = CREATE_BAMLIST.out.bamlist

        //
        // Join bamlist with contrasts
        //

        ch_contrasts = ch_contrasts
            .map { it -> [it['contrast'], it] }
            .combine ( ch_bamlist, by: 0 )
            .map { it -> it[1] + ['bam1_text': it[2], 'bam2_text': it[3]] }

        //
        // Create input channels
        //

        ch_contrasts_bamlist = ch_contrasts.map { [ it.contrast, it.treatment, it.control, it.meta1, it.meta2, it.bam1, it.bam2, it.bam1_text, it.bam2_text ] }

        //
        // Run rMATS prep step
        //

        RMATS_PREP (
            params.annotation_gtf,
            ch_contrasts_bamlist,
            params.rmats_read_len,
            params.rmats_splice_diff_cutoff,
            params.rmats_novel_splice_site,
            params.rmats_min_intron_len,
            params.rmats_max_exon_len
        )

        ch_rmats_temp = RMATS_PREP.out.rmats_temp

        //
        // Join rmats temp with contrasts
        //

        ch_contrasts = ch_contrasts
            .map { it -> [it['contrast'], it] }
            .combine ( ch_rmats_temp, by: 0 )
            .map { it -> it[1] + ['rmats_temp': it[2]] }

        //
        // Create input channels
        //

        ch_contrasts_bamlist = ch_contrasts.map { [ it.contrast, it.treatment, it.control, it.meta1, it.meta2, it.bam1, it.bam2, it.bam1_text, it.bam2_text, it.rmats_temp ] }

        //
        // Run rMATs post step
        //

        RMATS_POST (
            GFFREAD_TX2GENE.out.tx2gene,
            ch_contrasts_bamlist,
            params.rmats_read_len,
            params.rmats_splice_diff_cutoff,
            params.rmats_novel_splice_site,
            params.rmats_min_intron_len,
            params.rmats_max_exon_len,
            params.rmats_paired_stats
        )

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