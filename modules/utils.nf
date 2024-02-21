// -*- groovy -*-

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

// Get medaka model into the store_dir/medaka based on the model name
process GET_CLAIR3_MODEL {
    storeDir "${params.store_dir}/clair"
    
    input:
    val(model)

    output:
    path("${model}")

    script:
    """
curl -O https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/${model}.tar.gz
tar xvzf ${model}.tar.gz
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
samtools fastq -T'*' ${bam} | bgzip -@${task.cpus} > ${meta.name}.fastq.gz
    """
}

process CRAM_TO_FASTQ_GZ {
    input:
    tuple val(meta), path(cram)

    output:
    tuple val(meta), file("*.fastq.gz")

    script:
    """
samtools fastq -T'*' ${cram} | bgzip -@${task.cpus} > ${meta.name}.fastq.gz
    """
}
