process SKEW_PHASE {

    tag "$meta.id"
    label "process_low"

    conda "bioconda::bcftools=1.17"
    container "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"

    input:
    tuple val(meta), path(vcf)
    path skew_file

    output:
    tuple val(meta),
          path("${meta.id}.skew_phased.vcf.gz"),
          path("${meta.id}.skew_phased.vcf.gz.tbi"),
          path("${meta.id}.skew_phased.vcf.metrics.txt")

    script:
    """
    skew_phase.py \
        --vcf $vcf \
        --skew $skew_file \
        --out ${meta.id}.skew_phased.vcf

    bgzip ${meta.id}.skew_phased.vcf
    tabix -p vcf ${meta.id}.skew_phased.vcf.gz
    """
}