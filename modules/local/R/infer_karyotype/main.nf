process INFER_KARYOTYPE {
    tag "${meta.id}"

    container "ghcr.io/qgouil/skewx-r:0.2"

    input:
    tuple val(meta), path(vcf), path(vcf_tbi)

    output:
    tuple val(meta), path("${meta.id}_karyotype.tsv"), emit: karyotype

    script:
    """
    Rscript ${projectDir}/bin/infer_karyotype.R ${vcf} ${meta.id}
    """
}
