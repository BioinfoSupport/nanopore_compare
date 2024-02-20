// Downloads the model or gets the one already present in the container
process GET_MEDAKA_MODEL {
    storeDir "${params.store_dir}/medaka/${medaka_model}"
    
    input:
    val(medaka_model)

    output:
    path("*")

    script:
    """
HOME=. res=`medaka tools resolve_model --model ${medaka_model}`
if [ ! -d .medaka/data ]; then
    cp \$res .
else
    mv .medaka/data/* .
    rm -r .medaka
fi
    """
}

// Utility conversion
process BAM_TO_FASTA {
    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.fasta"), path("*.fasta.fai")

    script:
    """
samtools fasta ${bam} > ${bam.getBaseName()}.fasta
samtools faidx ${bam.getBaseName()}.fasta
    """
}

process BAM_TO_FASTQ_GZ {
    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.fastq.gz")

    script:
    """
samtools fastq ${bam} | bgzip > ${bam.getBaseName()}.fastq.gz
    """
}
