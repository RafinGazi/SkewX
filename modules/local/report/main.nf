process REPORT_INDIVIDUAL {

    tag "$meta.id"

    input:
    tuple val(meta),
          path(htmls),
          path(whatshap_stats),
          path(whatshap_blocks),
          path(clustered_reads),
          path(skew_tsv),
          path(karyotype_tsv),
          path(karyotype_plot),
          path(cgi_bed),
          path(report_template)

    output:
    path("${meta.id}_report.qmd"), emit: qmds
    path(htmls), emit: htmls
    path("_${whatshap_stats.baseName}.qmd"), emit: whatshap_stats
    path(whatshap_blocks), emit: whatshap_blocks
    path(clustered_reads), emit: clustered_reads
    path(skew_tsv), emit: skew_tsv
    path(karyotype_tsv), emit: karyotype_tsv
    path(karyotype_plot), emit: karyotype_plot

    script:
    """
    cp "${report_template}" "${meta.id}_report.qmd"

    # identifiers
    sed -i "s/ext_meta_id/${meta.id}/g" "${meta.id}_report.qmd"
    sed -i "s/ext_all_tissues_list/${meta.sample}/g" "${meta.id}_report.qmd"

    # file paths for R
    sed -i "s|ext_blocks_stats_file|${whatshap_blocks}|g" "${meta.id}_report.qmd"
    sed -i "s|ext_CGI_bed_file|${cgi_bed}|g" "${meta.id}_report.qmd"
    sed -i "s|ext_karyotype_tsv|${karyotype_tsv}|g" "${meta.id}_report.qmd"
    sed -i "s|ext_karyotype_plot|${karyotype_plot}|g" "${meta.id}_report.qmd"

    # format whatshap stats include
    echo '```' | cat - ${whatshap_stats} > "_${whatshap_stats.baseName}.qmd"
    echo '```' >> "_${whatshap_stats.baseName}.qmd"
    """
}


process REPORT_SKIPPED {

    tag "$meta.id"

    input:
    tuple val(meta),
          path(karyotype_tsv),
          path(karyotype_plot),
          path(report_template)

    output:
    path("${meta.id}_skipped_report.qmd"), emit: qmds

    script:
    """
    cp "${report_template}" "${meta.id}_skipped_report.qmd"

    sed -i "s/ext_meta_id/${meta.id}/g" "${meta.id}_skipped_report.qmd"
    sed -i "s|ext_karyotype_tsv|${karyotype_tsv}|g" "${meta.id}_skipped_report.qmd"
    sed -i "s|ext_karyotype_plot|${karyotype_plot}|g" "${meta.id}_skipped_report.qmd"
    """
}


process REPORT_BOOK {

    label "process_low"
    stageInMode "copy"
    publishDir "${params.outdir}", mode: "copy"
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
        path(skew_tsvs)
        path(karyotype_tsvs)
        path(karyotype_plots)
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