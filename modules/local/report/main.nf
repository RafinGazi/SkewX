process REPORT_INDIVIDUAL {

    tag "$meta.id"
    label "process_single"

    input:
    tuple val(meta),
          path(htmls),
          path(whatshap_stats),
          path(whatshap_blocks),
          path(clustered_reads_tsv),
          path(skew_tsv),
          path(karyotype_tsv),
          path(karyotype_plot),
          path(cohort_tsv),
          path(cohort_plot),
          path(cgi_bed),
          path(report_template),
          path(skew_phased_vcf),
          path(skew_metrics)

    output:
    path("${meta.id}_report.qmd"), emit: qmds
    path(htmls),                   emit: htmls
    path("_${whatshap_stats.baseName}.qmd"), emit: whatshap_stats
    path(whatshap_blocks),         emit: whatshap_blocks
    path(clustered_reads_tsv),     emit: clustered_reads
    path(skew_tsv),                emit: skew_tsv
    path(karyotype_tsv),           emit: karyotype_tsv
    path(karyotype_plot),          emit: karyotype_plot
    path(skew_phased_vcf),         emit: skew_phased_vcf 

    script:
    """
    # copy individual template
    cp "${report_template}" "${meta.id}_report.qmd"

    # substitute individual id into report
    sed -i "s/ext_meta_id/${meta.id}/g" "${meta.id}_report.qmd"

    # sub whatshap stats blocks file path into report
    sed -i "s/ext_blocks_stats_file/${whatshap_blocks}/g" "${meta.id}_report.qmd"

    # sub path to CGI bed file into each report
    sed -i "s/ext_CGI_bed_file/${cgi_bed}/g" "${meta.id}_report.qmd"

    # sub tissue names into report
    sed -i "s/ext_all_tissues_list/${meta.sample}/g" "${meta.id}_report.qmd"

    # sub karyotype tsv path into report
    sed -i "s/ext_karyotype_tsv/${karyotype_tsv}/g" "${meta.id}_report.qmd"

    # sub karyotype plot path into report
    sed -i "s/ext_karyotype_plot/${karyotype_plot}/g" "${meta.id}_report.qmd"

    # skew phased vcf path into report
    sed -i "s|ext_skew_phased_vcf|${skew_phased_vcf}|g" "${meta.id}_report.qmd"

    # skew metrics path into report
    sed -i "s|ext_skew_metrics|${skew_metrics}|g" "${meta.id}_report.qmd"
    
    # sub cohort tsv path into report
    sed -i "s/ext_cohort_tsv/${cohort_tsv}/g" "${meta.id}_report.qmd"

    # sub cohort plot path into report
    sed -i "s/ext_cohort_plot/${cohort_plot}/g" "${meta.id}_report.qmd"

    # turn text files into qmd for code formatting
    echo '```' | cat - ${whatshap_stats} > "_${whatshap_stats.baseName}.qmd"
    echo '```' >> "_${whatshap_stats.baseName}.qmd"
    """

}

process REPORT_BOOK {

    label "process_low"
    stageInMode "copy"
    publishDir "${params.outdir}", mode: "copy"
    conda "${moduleDir}/../R/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://qgouil/skewx/skewx-r:0.2' :
        'ghcr.io/qgouil/skewx-r:0.2' }"

    input:
        path(book_template_files)
        path(qmds)
        path(mosdepth_htmls)
        path(whatshap_stats)
        path(whatshap_blocks)
        path(clustered_reads)
        path(skews)
        path(karyotype_tsvs)
        path(karyotype_plots)
        path(cohort_tsvs)     
        path(cohort_plots)    
        path(cgi_bed)

    output:
        path("_book")  

    script:
    """
    # initialize _quarto.yml
    cp _quarto_template.yml _quarto.yml

    # append chapters (each patient) to _quarto.yml
    for rep in *_report.qmd
    do
        echo "    - \$rep" >> _quarto.yml
    done

    # add downloadthis quarto extension
    unzip ${projectDir}/assets/report-templates/_extensions.zip -d ./

    export XDG_CACHE_HOME="./quarto_cache"
    mkdir -p \$XDG_CACHE_HOME
    quarto render
    """
}