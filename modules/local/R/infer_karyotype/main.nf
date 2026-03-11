process INFER_KARYOTYPE {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.outdir}/karyotype", mode: "copy"

    conda "${moduleDir}/../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://qgouil/skewx/skewx-r:0.2' :
        'ghcr.io/qgouil/skewx-r:0.2' }"

    input:
    tuple val(meta), path(bed_gz), path(bed_gz_csi), path(summary_txt)

    output:
    tuple val(meta), path("${meta.id}_karyotype.tsv"),          emit: karyotype_tsv
    tuple val(meta), path("${meta.id}_karyotype_coverage.png"), emit: karyotype_plot

    script:
    """
    Rscript ${projectDir}/bin/infer_karyotype.R \\
        ${bed_gz} \\
        ${summary_txt} \\
        ${meta.id} \\
        ${meta.id}_karyotype.tsv
    """
}