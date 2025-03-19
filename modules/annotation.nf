
process PROKKA_ANNOTATE {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)/meta.name}"
    
    input:
    tuple val(meta), path(fasta), path(fai), path(mmi)

    output:
    tuple val(meta), path("prokka_annotation")

    script:
    """
prokka --outdir prokka_annotation --prefix ${meta.name} --cpus ${task.cpus} ${params.prokka_args} ${fasta}
    """
}


process ADD_GENOME_JSON {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)/meta.name}"
    cpus 1

    input:
    tuple val(meta), path("prokka_annotation")

    output:
    tuple val(meta), path("${meta.name}.genome.json")
    path("prokka_annotation/${meta.name}.fsa.fai")

    script:
    """
samtools faidx prokka_annotation/${meta.name}.fsa
cat > ${meta.name}.genome.json <<EOF
{
  "id": "${meta.name}",
  "name": "${meta.name}",
  "fastaURL": "prokka_annotation/${meta.name}.fsa",
  "indexURL": "prokka_annotation/${meta.name}.fsa.fai",
  "tracks": [
    {
      "name": "Prokka Genes",
      "format": "gff",
      "url": "prokka_annotation/${meta.name}.gff"
    },
    {
      "name": "Liftover ${params.ref_id}",
      "format": "gff",
      "url": "${params.ref_id}.lifted.gff"
    }
  ]
}
EOF
    """
}

process LIFTOVER_REF_ANNOTATIONS {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)/meta.name}"
    cpus 1
    
    input:
    tuple val(ref_meta), path(ref_fasta), path(ref_fai)
    tuple val(meta), path(fasta), path(fai), path(mmi)
    path(ref_gff)

    output:
    path("${ref_meta.name}.lifted.gff")
    
    script:
    """
minimap2 -un -cx asm5 ${ref_fasta} ${fasta} >liftover.paf
paf2chain --input liftover.paf >liftover.chain
CrossMap gff liftover.chain ${ref_gff} lifted.gff
(
  head -1 ${ref_gff}
  cat lifted.gff
) > ${ref_meta.name}.lifted.gff
    """

}

