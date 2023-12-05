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
include { DEXSEQ_COUNT } from './modules_simple/dexseq_count.nf'
include { DEXSEQ_EXON } from './modules_simple/dexseq_exon.nf'
include { DRIMSEQ_FILTER  } from './modules_simple/drimseq_filter.nf'
include { DEXSEQ_DTU      } from './modules_simple/dexseq_dtu.nf'

workflow{

read_pairs_ch = Channel
    .fromPath( params.csv_input )
    .splitCsv(header: true, sep: ',')
    .map {row -> tuple(row.sample, [row.fastq_1, row.fastq_2], row.condition)}
    .view()

//CAT_FASTQ(read_pairs_ch)

FASTQC( read_pairs_ch )
//genome
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

//MULTIQC(SALMON_QUANT.out, TRIMGALORE.out, STAR_ALIGN.out, FASTQC.out)
//MULTIQC(SALMON_QUANT.out.collect(), TRIMGALORE.out.collect(), STAR_ALIGN.out.collect(), FASTQC.out.collect())
MULTIQC(SALMON_QUANT.out.json_info.collect(), TRIMGALORE.out.log.collect(), STAR_ALIGN.out.log_final.collect(), FASTQC.out.zip.collect())


//SALMON_QUANT.out.transcripts.collect(), TRIMGALORE.out.reads.collect(), STAR_ALIGN.out.collect(), FASTQC.out.zip

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

MERGE_RESULTS_SALMON(SALMON_QUANT.out.transcripts.collect())

TXIMPORT ( MERGE_RESULTS_SALMON.out.merged, GFFREAD_TX2GENE.out.tx2gene )

DRIMSEQ_FILTER( TXIMPORT.out.txi_dtu, TXIMPORT.out.tximport_tx2gene, params.csv_input, params.min_samps_gene_expr, params.min_samps_feature_expr, params.min_samps_feature_prop, params.min_feature_expr, params.min_feature_prop, params.min_gene_expr )

DEXSEQ_DTU(DRIMSEQ_FILTER.out.drimseq_samples_tsv, DRIMSEQ_FILTER.out.drimseq_counts_tsv, params.csv_contrastsheet, params.n_dexseq_plot)

}

process MERGE_RESULTS_SALMON {
    publishDir params.outdir
    container "zavolab/salmon:1.1.0"

    input:
    path out_folders
    
    output:
    path("salmon"), emit:merged
    
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