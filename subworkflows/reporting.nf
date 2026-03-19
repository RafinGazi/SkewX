include {MOSDEPTH_PLOTDIST} from "../modules/local/mosdepth/plotdist/main.nf"
include {NANOCOMP} from "../modules/local/nanocomp/main.nf"
include {REPORT_INDIVIDUAL} from "../modules/local/report/main.nf"
include {REPORT_BOOK} from "../modules/local/report/main.nf"

workflow reporting {
    take:
        mosdepth_report_results // channel containing by-sample and merged mosdepth report results
        haplotagged_samples     // channel containing by-sample and merged sample haplotagged bams
        whatshap_stats_blocks   // channel containing phased vcfs block stats
        clustered_reads         // channel containing clustered reads and skew information
        cgi_bed                 // single-item channel containing CGI bed file
        karyotype_tsv           // channel containing karyotype TSV per individual
        karyotype_plot          // channel containing karyotype coverage plot per individual
        cohort_tsv              // channel containing all karyotype's cohort tsv

    main:
        // prepare mosdepth coverage report
        ch_global_dist_bysample = mosdepth_report_results
            .map{ it -> tuple(it[0], it[1]) }

        ch_mosdepth_dist_report = MOSDEPTH_PLOTDIST(ch_global_dist_bysample)

        // prepare nanocomp reports
        ch_nanocomp = NANOCOMP(haplotagged_samples)

        // combine all inputs for individual report
        ch_combined_qc_reports = ch_nanocomp
            .map{ it -> tuple(it[0].id, it[0].sample, it[1]) }
            .join(ch_mosdepth_dist_report.map{ it -> tuple(it[0].id, it[0].sample, it[1]) })
            .join(whatshap_stats_blocks.map{ it -> tuple(it[0].id, it[0].sample, it[1], it[2]) })
            .join(clustered_reads.map{ it -> tuple(it[0].id, it[0].sample, it[1], it[2]) }.groupTuple())
            .join(karyotype_tsv.map{ it -> tuple(it[0].id, it[1]) })
            .join(karyotype_plot.map{ it -> tuple(it[0].id, it[1]) })
            .join(cohort_tsv.map{ it -> tuple(it[0].id, it[1]) })
            .map{ it -> tuple(
                [id: it[0], sample: it[1]],
                it[2] + [it[4]],   // htmls
                it[6],             // whatshap_stats
                it[7],             // whatshap_blocks
                it[9],             // clustered_reads
                it[10],            // skew_tsv
                it[11],            // karyotype_tsv
                it[12],             // karyotype_plot
                it[13]              // cohort_tsv
            )}
            .combine(cgi_bed.map{ it -> it[1] })
            .combine(channel.fromPath("${projectDir}/assets/report-templates/individual_report.qmd", checkIfExists: true))

        ch_reporting_files = REPORT_INDIVIDUAL(ch_combined_qc_reports)

        // create channel for templates
        ch_book_template_files = channel.fromPath([
            "${projectDir}/assets/report-templates/_quarto_template.yml",
            "${projectDir}/assets/report-templates/index.qmd"
        ], checkIfExists: true).collect()

        book = REPORT_BOOK(
            ch_book_template_files,
            ch_reporting_files.qmds.collect(),
            ch_reporting_files.htmls.collect(),
            ch_reporting_files.whatshap_stats.collect(),
            ch_reporting_files.whatshap_blocks.collect(),
            ch_reporting_files.clustered_reads.collect(),
            ch_reporting_files.skew_tsv.collect(),
            ch_reporting_files.karyotype_tsv.collect(),
            ch_reporting_files.karyotype_plot.collect(),
            ch_reporting_files.cohort_tsv.collect(),
            cgi_bed.map{ it -> it[1] }
        )

    emit:
        book
}