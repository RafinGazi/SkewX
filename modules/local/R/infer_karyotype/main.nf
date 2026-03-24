process INFER_KARYOTYPE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/karyotype", mode: "copy"
    input:
    tuple val(meta), path(bed_gz), path(bed_gz_csi), path(summary_txt)
    output:
    tuple val(meta), path("${meta.id}_karyotype.tsv"),          emit: karyotype_tsv
    tuple val(meta), path("${meta.id}_karyotype_coverage.png"), emit: karyotype_plot
    script:
    """
    /usr/bin/python3 ${projectDir}/bin/infer_karyotype_v2.py \\
        ${bed_gz} \\
        ${summary_txt} \\
        ${meta.id} \\
        ${meta.id}_karyotype.tsv
    """
}
