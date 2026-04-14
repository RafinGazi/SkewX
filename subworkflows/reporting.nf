include {MOSDEPTH_PLOTDIST} from "../modules/local/mosdepth/plotdist/main.nf"
include {NANOCOMP} from "../modules/local/nanocomp/main.nf"
include {REPORT_INDIVIDUAL} from "../modules/local/report/main.nf"
include {REPORT_SKIPPED} from "../modules/local/report/main.nf"
include {REPORT_BOOK} from "../modules/local/report/main.nf"

workflow reporting {

    take:
        mosdepth_report_results
        haplotagged_samples
        whatshap_stats_blocks
        clustered_reads
        cgi_bed
        karyotype_tsv
        karyotype_plot
        karyotype_all_plots
        ch_pass
        ch_skipped
        ch_flagged

    main:

    /*
    ---------------- PASS ----------------
    */

    ch_global_dist_bysample = mosdepth_report_results
        .map { tuple(it[0], it[1]) }

    ch_mosdepth_dist_report = MOSDEPTH_PLOTDIST(ch_global_dist_bysample)
    ch_nanocomp = NANOCOMP(haplotagged_samples)

    ch_combined_pass = ch_nanocomp
        .map { tuple(it[0].id, it[0].sample, it[1]) }
        .join(ch_mosdepth_dist_report.map { tuple(it[0].id, it[0].sample, it[1]) })
        .join(whatshap_stats_blocks.map { tuple(it[0].id, it[0].sample, it[1], it[2]) })
        .join(clustered_reads.map { tuple(it[0].id, it[0].sample, it[1], it[2]) }.groupTuple())
        .join(karyotype_tsv.map { tuple(it[0].id, it[1]) })
        .join(karyotype_plot.map { tuple(it[0].id, it[1]) })
        .join(karyotype_all_plots.map { tuple(it[0].id, it[1]) })

        .map { it ->
        tuple(
            [id: it[0], sample: it[1]],
            it[2] + [it[4]],    // htmls (nanocomp + mosdepth)
            it[6],              // whatshap stats
            it[7],              // whatshap blocks
            it[9],              // clustered reads
            it[10],             // skew tsv
            it[11],             // karyotype tsv
            it[12],             // karyotype plot
            it[13] instanceof List ? it[13] : [it[13]]              // karyotype all plots
        )
    }

        .combine(cgi_bed.map { it[1] })
        .combine(channel.fromPath("${projectDir}/assets/report-templates/individual_report.qmd", checkIfExists: true))

    ch_pass_reports = REPORT_INDIVIDUAL(ch_combined_pass)


    /*
    ---------------- SKIPPED ----------------
    */

    ch_skipped = ch_skipped.mix(ch_flagged)

    ch_combined_skipped = ch_skipped
        .map { meta -> tuple(meta.id, meta.sample, meta.qc_flag) }
        .join(karyotype_tsv.map { tuple(it[0].id, it[1]) })
        .join(karyotype_plot.map { tuple(it[0].id, it[1]) })
        .join(karyotype_all_plots.map { tuple(it[0].id, it[1]) })
        .map { id, sample, qc_flag, karyo_tsv, karyo_plot, karyo_all_plots ->
            tuple(
                [id: id, sample: sample, qc_flag: qc_flag],
                karyo_tsv,
                karyo_plot,
                karyo_all_plots instanceof List ? karyo_all_plots : [karyo_all_plots]
            )
        }
        .combine(channel.fromPath("${projectDir}/assets/report-templates/skipped_report.qmd", checkIfExists: true))

    ch_skipped_reports = REPORT_SKIPPED(ch_combined_skipped)


    /*
    ---------------- BOOK ----------------
    */

    ch_all_qmds = ch_pass_reports.qmds.mix(ch_skipped_reports.qmds)

    ch_book_template_files = channel.fromPath([
        "${projectDir}/assets/report-templates/_quarto_template.yml",
        "${projectDir}/assets/report-templates/index.qmd"
    ], checkIfExists: true).collect()

    book = REPORT_BOOK(
        ch_book_template_files,
        ch_all_qmds.collect(),
        ch_pass_reports.htmls.collect().ifEmpty([]),
        ch_pass_reports.whatshap_stats.collect().ifEmpty([]),
        ch_pass_reports.whatshap_blocks.collect().ifEmpty([]),
        ch_pass_reports.clustered_reads.collect().ifEmpty([]),
        ch_pass_reports.skew_tsv.collect().ifEmpty([]),
        ch_pass_reports.karyotype_tsv.mix(ch_skipped_reports.karyotype_tsv).collect().ifEmpty([]),
        ch_pass_reports.karyotype_plot.mix(ch_skipped_reports.karyotype_plot).collect().ifEmpty([]),
        ch_pass_reports.karyotype_all_plots.mix(ch_skipped_reports.karyotype_all_plots).collect().ifEmpty([]),
        cgi_bed.map { it[1] }
    )

    emit:
        book
}