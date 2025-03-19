// -*- groovy -*-

// Just direct mapping of the reads to the reference, for IGV manual control
process MAP_READS {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)/meta.name}"

    input:
    tuple val(ref_meta), path(ref_fa), path(ref_fa_fai), path(ref_fa_mmi)
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path(bam), path(bam_bai)

    script:
    bam = "reads_vs_${ref_meta.name}.bam"
    bam_bai = bam+".bai"
    """
minimap2 -ax map-ont -t${task.cpus} ${params.minimap2_cs_tag} ${ref_fa_mmi} ${fastq} | \
    samtools sort -@${task.cpus} -o ${bam} -
samtools index ${bam}
    """
}

process BAMSTATS_MAPPED_READS {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)/meta.name}"

    input:
    tuple val(meta), path(bam), path(bam_bai)

    output:
    tuple val(meta), path(bamstats)

    script:
    bamstats = bam + ".bamstats"
    """
samtools stats ${bam} > ${bamstats}
    """
}

process FLAGSTAT_MAPPED_READS {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)/meta.name}"

    input:
    tuple val(meta), path(bam), path(bam_bai)

    output:
    tuple val(meta), path(flagstat)

    script:
    flagstat = bam + ".flagstat"
    """
samtools flagstat ${bam} > ${flagstat}
    """
}
