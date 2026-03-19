process COHORT_KARYOTYPE_QC {

    tag "cohort_karyotype_qc"

    publishDir "${params.outdir}/karyotype", mode: "copy"

    container null

    input:
    path(karyotype_tsvs)

    output:
    path("cohort_qc.tsv"), emit: qc_tsv
    path("cohort_rx_ry.png"), emit: qc_plot

    script:
    """
    python3 ${projectDir}/bin/cohort_karyotype_qc.py \
        ${karyotype_tsvs.join(' ')} \
        --out_prefix cohort
    """
}