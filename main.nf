#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SkewX: A Nextflow pipeline for skewed X inactivation analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/QGouil/SkewX
    Publication: https://doi.org/10.1101/gr.279396.124
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Nanopore reads PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.input) {
    ch_input = Channel.fromPath(params.input, checkIfExists: true)
} else { exit 1, 'Input sample sheet not specified!' }
if (params.reference) {
    ch_reference = Channel.fromPath("${params.reference}", checkIfExists: true).map{
        it -> tuple(id: it.baseName, it, "${it}.fai")
    }
} else { exit 1, "Reference FASTA not specified!" }
if (params.cgi_bedfile) {
    ch_cgibed = Channel.fromPath(params.cgi_bedfile, checkIfExists: true)
        .map{
            it -> tuple(id: it.baseName, it)
        }
} else { exit 1, "CGI Bed file not specified!"}
if (!(params.deepvariant_model in ["WGS", "WES", "PACBIO", "ONT_R104", "HYBRID_PACBIO_ILLUMINA"])) {
    exit 1, "DeepVariant model must be one of WGS, WES, PACBIO, ONT_R104, or HYBRID_PACBIO_ILLUMINA"
}
if (!(params.stage in ["raw", "phased", "haplotagged"])) {
    exit 1, "Stage must be one of: raw, phased, haplotagged"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Define processes and modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include {INPUT_CHECK} from './subworkflows/local/input_check.nf'
include {MINIMAP2} from "./modules/local/minimap2/main.nf"
include {SAMTOOLS_MERGE} from "./modules/local/samtools/merge/main.nf"
include {SAMTOOLS_INDEX as SAMTOOLS_INDEX_SAMPLES} from "./modules/local/samtools/index/main.nf"
include {SAMTOOLS_INDEX as SAMTOOLS_INDEX_MERGED} from "./modules/local/samtools/index/main.nf"
include {SAMTOOLS_INDEX as SAMTOOLS_INDEX_HAPLOTAG} from "./modules/local/samtools/index/main.nf"
include {SAMTOOLS_INDEX as SAMTOOLS_INDEX_HAPLOTAG_MERGED} from "./modules/local/samtools/index/main.nf"
include {FILTER_PASS} from "./modules/local/bcftools/view_pass/main.nf"
include {WHATSHAP_PHASE} from "./modules/local/whatshap/phase/main.nf"
include {WHATSHAP_STATS} from "./modules/local/whatshap/stats/main.nf"
include {WHATSHAP_HAPLOTAG} from "./modules/local/whatshap/haplotag/main.nf"
include {WHATSHAP_HAPLOTAG as WHATSHAP_HAPLOTAG_MERGED} from "./modules/local/whatshap/haplotag/main.nf"
include {MOSDEPTH} from "./modules/local/mosdepth/main.nf"
include {MOSDEPTH as MOSDEPTH_MERGED} from "./modules/local/mosdepth/main.nf"
include {SAMTOOLS_VIEWHP} from "./modules/local/samtools/view_hp/main.nf"
include {R_CLUSTERBYMETH} from "./modules/local/R/cluster_by_meth/main.nf"
include {reporting} from "./subworkflows/reporting.nf"
include {separated_deepvariant} from "./subworkflows/local/deepvariant/main.nf"
include { INFER_KARYOTYPE } from './modules/local/py/infer_karyotype/main'

//
// WORKFLOW: Run main SkewX analysis pipeline
//
workflow SKEWX {

    /*
    =====================================================
    INPUT PREPARATION
    =====================================================
    */

    ch_checked_input = INPUT_CHECK(ch_input)

    ch_separate_samples = ch_checked_input
        .map { individual, sample, bam, vcf ->
            tuple([id: individual, sample: sample], bam)
        }

    ch_vcf_per_individual = ch_checked_input
        .map { individual, sample, bam, vcf ->
            tuple(individual, file(vcf))
        }

    if (params.ubam) {
        ch_aligned = MINIMAP2(ch_separate_samples, ch_reference)
    } else {
        ch_aligned = ch_separate_samples
    }

    ch_samples = SAMTOOLS_INDEX_SAMPLES(ch_aligned)

    ch_grouped_bams = ch_samples
        .map { meta, bam, bai ->
            tuple(meta.id, meta.sample, bam, bai)
        }
        .groupTuple()
        .map { individual, samples, bams, bais ->
            tuple([id: individual, sample: samples], bams, bais)
        }

    (ch_multiple_bams, ch_single_bams) = ch_grouped_bams.branch {
        multi: it[1].size() > 1
        single: it[1].size() == 1
    }

    ch_single_bams = ch_single_bams.map { meta, bams, bais ->
        tuple(meta, bams[0], bais[0])
    }

    SAMTOOLS_MERGE(ch_multiple_bams)
        | SAMTOOLS_INDEX_MERGED
        | mix(ch_single_bams)
        | set { ch_merged_bam }

    /*
    =====================================================
    STAGE ROUTING
    =====================================================
    */

    if (params.stage == "haplotagged") {

        log.info "Stage: haplotagged — skipping DeepVariant and Whatshap"

        ch_samples_haplotag = ch_merged_bam
        ch_whatshap_stats_blocks = Channel.empty()

    } else {

        if (params.stage == "raw") {

            ch_reference_rep_merged = ch_merged_bam
                .combine(ch_reference.collect())
                .map { meta, merged_bam, merged_bam_idx, meta_ref, ref, ref_idx ->
                    tuple(meta_ref, ref, ref_idx)
                }

            dv_args = Channel.from([
                regions: params.deepvariant_region,
                model_type: params.deepvariant_model,
                num_shards: params.deepvariant_num_shards
            ])

            (ch_vcf, _) = separated_deepvariant(
                dv_args,
                ch_merged_bam,
                ch_reference
            )
            .multiMap {
                vcf: tuple(it[0], it[1], it[2], it[5])
            }

            ch_vcf_pass = FILTER_PASS(ch_vcf)

            ch_vcf_phased = WHATSHAP_PHASE(
                ch_vcf_pass,
                ch_reference_rep_merged
            )

        }

        if (params.stage == "phased") {

            log.info "Stage: phased — running Whatshap phasing"

            // Prepare input VCF + BAM
            ch_vcf_input = ch_merged_bam
                .map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }
                .join(ch_vcf_per_individual, by: 0)
                .map { id, meta, bam, bai, vcf ->
                    tuple(meta, bam, bai, vcf, "${vcf}.tbi")
                }

            // Prepare reference (same format as raw stage)
            ch_reference_rep_merged = ch_merged_bam
                .combine(ch_reference.collect())
                .map { meta, merged_bam, merged_bam_idx, meta_ref, ref, ref_idx ->
                    tuple(meta_ref, ref, ref_idx)
                }

            // Run Whatshap phasing
            ch_vcf_phased = WHATSHAP_PHASE(
                ch_vcf_input,
                ch_reference_rep_merged
            )            

        }

        ch_whatshap_stats_blocks = WHATSHAP_STATS(
            ch_vcf_phased.map { meta, bam, bam_idx, vcf, vcf_idx ->
                tuple(meta, vcf, vcf_idx)
            }
        )

        (ch_tmp_samples, ch_reference_rep) = ch_samples
            .map { meta, bam, bam_idx ->
                tuple(meta.id, meta.sample, bam, bam_idx)
            }
            .combine(
                ch_vcf_phased.map {
                    meta, merged_bam, merged_bam_idx, vcf, vcf_idx ->
                        tuple(meta.id, meta.sample, vcf, vcf_idx)
                },
                by: 0
            )
            .combine(ch_reference.collect())
            .multiMap { id, single_sample, bam, bai,
                        samples, vcf, vcf_idx,
                        ref_id, ref, ref_idx ->

                samples: tuple([id: id, sample: single_sample], bam, bai, vcf, vcf_idx)
                ref: tuple(ref_idx, ref, ref_idx)
            }

        WHATSHAP_HAPLOTAG(ch_tmp_samples, ch_reference_rep)
            | SAMTOOLS_INDEX_HAPLOTAG
            | set { ch_samples_haplotag }

    }

    /*
    =====================================================
    DOWNSTREAM ANALYSIS (COMMON TO ALL STAGES)
    =====================================================
    */

    // MOSDEPTH runs on full BAM (chrX_Y_18_21)
    // .bed     → INFER_KARYOTYPE
    // .dist    → reporting (coverage plots)
    // .summary → INFER_KARYOTYPE
    MOSDEPTH(ch_samples_haplotag)

    // Karyotype 
    ch_karyotype = INFER_KARYOTYPE(
        MOSDEPTH.out.bed
            .join(MOSDEPTH.out.summary, by: 0)
            .map { meta, bed_gz, bed_gz_csi, summary_txt ->
                tuple(meta, bed_gz, bed_gz_csi, summary_txt)
            }
    )

    ch_karyotype_qc = ch_karyotype.karyotype_tsv
        .map { meta, tsv ->
            def lines = tsv.text.readLines()
            def header = lines[0].split('\t')
            def values = lines[1].split('\t')

            def qc_idx = header.indexOf("qc_flag")
            def qc_flag = values[qc_idx]

            tuple(meta + [qc_flag: qc_flag])
        }

    ch_branch = ch_karyotype_qc.branch {
        pass:    it.qc_flag == "pass"
        skipped: it.qc_flag.startsWith("skipped")
        flagged: it.qc_flag.startsWith("flagged")
    }


    ch_samples_haplotag_pass = ch_samples_haplotag
        .join(ch_branch.pass.map{ meta -> tuple(meta.id, true) }, by: 0)
        .map { id, meta, bam, bai, _ -> tuple(meta, bam, bai) }

    (ch_tmp_samples_haplotag, ch_cgibed_rep) = ch_samples_haplotag_pass
        .combine(ch_cgibed.collect())
        .multiMap { it ->
            samples_haplotag: tuple(it[0], it[1], it[2])
            cgibed: tuple(it[3], it[4])
        }

    ch_hpreads = SAMTOOLS_VIEWHP(ch_tmp_samples_haplotag, ch_cgibed_rep)

    ch_clustered_reads = R_CLUSTERBYMETH(ch_hpreads, ch_cgibed_rep)


    if (params.stage != "haplotagged") {
        book = reporting(
            MOSDEPTH.out.dist,
            ch_samples_haplotag,
            ch_whatshap_stats_blocks,
            ch_clustered_reads,
            ch_cgibed,
            ch_karyotype.karyotype_tsv,
            ch_karyotype.karyotype_plot
        )
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    SKEWX ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
