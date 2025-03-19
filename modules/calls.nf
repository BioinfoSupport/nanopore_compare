// -*- groovy -*-

// Calling by medaka from reads vs sequence
process MEDAKA_CALL {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)/meta.name}",
	saveAs: { it.replaceFirst(/ref/, ref_meta.name) }
    cpus 2

    input:
    tuple val(ref_meta), path(ref_fa), path(ref_fa_fai), path(ref_fa_mmi)
    tuple val(meta), path(fastq)
    path(medaka_model_path)

    output:
    tuple val(meta), path("medaka.reads_vs_ref.vcf.gz"), path("medaka.reads_vs_ref.vcf.gz.csi")
    
    script:
    """
medaka_variant -i ${fastq} -r ${ref_fa} -m ${medaka_model_path} \
        -t ${task.cpus} -o medaka.reads_vs_ref
## Fix of strange annotation mistake -- medaka duplicates some vcf lines?
## Also attempt to "standardize" the calls
echo '${meta.name}' >sample.txt
bcftools sort medaka.reads_vs_ref/medaka.annotated.vcf | \
    bcftools reheader -s sample.txt - | \
    bcftools norm -d exact -a - -o medaka.reads_vs_ref.vcf.gz
bcftools index medaka.reads_vs_ref.vcf.gz
    """
    stub:
    """
mkdir medaka.reads_vs_ref
echo ${ref_id} ${meta} | bgzip > medaka.reads_vs_ref.vcf.gz
touch medaka.reads_vs_ref.vcf.gz.csi
    """
}


// Call using minmap2/paftools form two assemblies
process MM2_CALL {
    // publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)/meta.name}",
    // 	saveAs: { it.replaceFirst(/ref/, ref_meta.name) }

    input:
    tuple val(ref_meta), path(ref_fa), path(ref_fa_fai)
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("mm2.cons_vs_ref.vcf.gz"), path("mm2.cons_vs_ref.vcf.gz.csi")
    tuple val(meta), path("mm2.cons_vs_ref.paf.gz")
    
    script:
    """
minimap2 -cx asm5 --cs -z1000000 ${ref_fa} ${fastq} | \
    sort -k6,6 -k8,8n | paftools.js call -s ${meta.name} -l0 -L0 -q0 - | \
    bgzip > mm2.cons_vs_ref.paf.gz

minimap2 -cx asm5 --cs -z1000000 ${ref_fa} ${fastq} | \
    sort -k6,6 -k8,8n | paftools.js call -s ${meta.name} -f ${ref_fa} -l0 -L0 -q0 - | \
    bcftools sort - |
    bcftools norm -a - -o mm2.cons_vs_ref.vcf.gz
bcftools index mm2.cons_vs_ref.vcf.gz
"""
    stub:
    """
echo ${ref_id} ${meta} | bgzip > mm2.cons_cs_ref.paf.gz
echo ${ref_id} ${meta} | bgzip > mm2.cons_cs_ref.vcf.gz
touch mm2.cons_cs_ref.vcf.gz.csi
"""
}



// Call using minmap2/paftools form two assemblies
process MM2_CALL_ANNOTATE {
    label 'Rscript'
    
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)/meta.name}",
	saveAs: { it.replaceFirst(/ref/, ref_meta.name) }

    input:
    tuple val(ref_meta), path(ref_fa), path(ref_fa_fai)
    tuple val(meta), path(vcf), path(bam)

    output:
    tuple val(meta), path("mm2.cons_vs_ref.annotated.vcf.gz"), path("mm2.cons_vs_ref.annotated.vcf.gz.csi")
    
    script:
    """

annotate_consensus_vcf.R ${vcf} ${bam} mm2.cons_vs_ref.annotated.vcf
bgzip mm2.cons_vs_ref.annotated.vcf
bcftools index mm2.cons_vs_ref.annotated.vcf.gz
    """
}

// Call using minmap2/paftools form two assemblies
process MM2_CALL_STDOUT {

    input:
    tuple val(ref_meta), path("ref.fasta"), path("ref.fasta.fai")
    tuple val(meta), path(fastq)

    output:
    stdout
    
    script:
    """
minimap2 -cx asm5 --cs -z1000000 ref.fasta ${fastq} | \
    sort -k6,6 -k8,8n | paftools.js call -s ${meta.name} -f ref.fasta -l0 -L0 -q0 - | \
    bcftools sort -
"""
}


// Calling by medaka from reads vs sequence
process CLAIR3_CALL {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)/meta.name}",
	saveAs: { it.replaceFirst(/clair\/merge_output/, "clair_vs_"+ref_meta.name) }

    input:
    tuple val(ref_meta), path(ref_fa), path(ref_fa_fai)
    tuple val(meta), path(bam), path(bam_bai)
    path(clair_model_path)

    output:
    tuple val(meta), path("clair/merge_output.vcf.gz"), path("clair/merge_output.vcf.gz.tbi")
    
    script:
    """
/opt/bin/run_clair3.sh \
  --bam_fn=${bam} \
  --ref_fn=${ref_fa} \
  --threads=${task.cpus} \
  --platform="${params.clair3_platform}" \
  --include_all_ctgs \
  --haploid_precise \
  --model_path="${clair_model_path}" \
  --sample_name="${meta.name}" \
  --output="`pwd`/clair"               ## absolute output path prefix
    """
}


// Calling by medaka from reads vs sequence
process CLAIR3_CALL_SENSITIVE {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)/meta.name}",
	saveAs: { it.replaceFirst(/clair\/merge_output/, "clair_sensitive_vs_"+ref_meta.name) }

    input:
    tuple val(ref_meta), path(ref_fa), path(ref_fa_fai)
    tuple val(meta), path(bam), path(bam_bai)
    path(clair_model_path)

    output:
    tuple val(meta), path("clair/merge_output.vcf.gz"), path("clair/merge_output.vcf.gz.tbi")
    
    script:
    """
/opt/bin/run_clair3.sh \
  --bam_fn=${bam} \
  --ref_fn=${ref_fa} \
  --threads=${task.cpus} \
  --platform="${params.clair3_platform}" \
  --include_all_ctgs \
  --haploid_sensitive \
  --model_path="${clair_model_path}" \
  --sample_name="${meta.name}" \
  --output="`pwd`/clair"               ## absolute output path prefix
    """
}


// Calling SV by sniffles
process SNIFFLES_CALL {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)/meta.name}"

    input:
    tuple val(ref_meta), path(ref_fa), path(ref_fa_fai)
    tuple val(meta), path(bam), path(bam_bai)

    output:
    tuple val(meta), path("sniffles.vcf.gz")
    
    script:
    """
sniffles \
  --input ${bam} \
  --reference ${ref_fa} \
  --threads ${task.cpus} \
  --vcf sniffles.vcf
bgzip sniffles.vcf
    """
}


process ANNOTATE_VCF_BY_REGIONS {
    input:
    tuple path(vcf), path(vcf_csi)
    path(regions)

    output:
    tuple path("region_annotated.vcf.gz"), path("region_annotated.vcf.gz.csi")

    script:
    """
echo '##INFO=<ID=NOTE,Number=1,Type=String,Description="Region manual note">' > annot.hdr
## Fix potential bad regions from igv
awk 'BEGIN{OFS="\t"}; NF<4{\$4="Empty note"}; {print}' ${regions} > fixed.${regions}

bcftools annotate -a fixed.${regions} -c CHROM,FROM,TO,NOTE -h annot.hdr ${vcf} -o region_annotated.vcf.gz
bcftools index region_annotated.vcf.gz
    """
}

// Align the final consensus to the reference.
// Can be further used for IGV display or various remaping operations
process REMAP_CONSENSUS_TO_REFERENCE {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)/meta.name}",
	saveAs: { it.replaceFirst(/ref/, ref_meta.name) }

    input:
    tuple val(ref_meta), path(ref_fa), path(ref_fa_fai)
    tuple val(meta), path(consensus_fastx)

    output:
    tuple val(meta), path("*.bam"), path ("*.bai")

    script:
    """
minimap2 -cx asm5 --cs -z1000000 -a ${ref_fa} ${consensus_fastx} | \
        samtools sort -o consensus_vs_ref.bam -
samtools index consensus_vs_ref.bam
    """
    stub:
    """
touch consensus_vs_ref.bam consensus_vs_ref.bam.bai
    """
}


// Merging Medaka haploid calls and minimapo2 calls from consensus
// In coinciding case Medaka call is stored
// Short consensus calls are filtered out
// TODO: Filter out MEDAKA calls within large consensus deletions
process MERGE_MEDAKA_AND_CONS_CALLS {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)/meta.name}",
	saveAs: { it.replaceFirst(/ref/, ref_meta.name+suffix) }

    input:
    tuple val(ref_meta), path(ref_fasta), path(ref_fai)
    tuple val(meta), path(medaka_vcf), path(medaka_csi), path(cons_vcf), path(cons_csi)
    val(suffix)
    
    output:
    tuple val(meta), path("joined_vs_ref.vcf.gz"), path("joined_vs_ref.vcf.gz.csi")
    
    script:
    // TODO Move thresholds to parameters?
    """
bcftools query -f'%CHROM\\t%POS0\\t%END\\t%ID\\n' -i 'STRLEN(REF)>50' ${cons_vcf} >long_sv.bed
if [ -s long_sv.bed ]
then
    bcftools view -T ^long_sv.bed ${medaka_vcf} -o medaka_not_in_sv.vcf.gz
else
    bcftools view ${medaka_vcf} -o medaka_not_in_sv.vcf.gz
fi
bcftools index medaka_not_in_sv.vcf.gz

bcftools concat -a -D medaka_not_in_sv.vcf.gz ${cons_vcf} | \
    bcftools filter -e 'QNAME="." & QUAL<${params.medaka_lowq}' -s MEDAKA_LQ | \
    bcftools filter -e 'QNAME="." & DP<=3' -s MEDAKA_LDP | \
    bcftools filter -e 'QNAME!="." & abs(strlen(ALT)-strlen(REF))<50' -s CONS_Short \
        -o joined_vs_ref.vcf.gz
bcftools index joined_vs_ref.vcf.gz
    """
}



// Region annotation bz regions bed should be made at table level -- not quite good for VCF
process MERGE_JOINED_CALL_FILES {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)}",
	saveAs: { it.replaceFirst(/ref/, ref_meta.name+suffix) }

    input:
    path 'vcf*.vcf.gz'
    path regions_annotation // Optional
    tuple val(ref_meta), path(ref_fasta), path(ref_fai)
    val(suffix) // Indicate lifted for publishing

    output:
    tuple path('joined_vs_ref.merged.vcf.gz'), path('joined_vs_ref.merged.vcf.gz.csi')
    tuple path('joined_vs_ref.merged.discordant.vcf.gz'), path('joined_vs_ref.merged.discordant.vcf.gz.csi'), optional: true

    script:
    """
for f in vcf*.vcf.gz; do
    bcftools index \$f
done
if [ -f vcf.vcf.gz ]; then # Single file given, no merge needed
    cp vcf.vcf.gz joined_vs_ref.merged0.vcf.gz
    mv vcf.vcf.gz.csi joined_vs_ref.merged0.vcf.gz.csi
else
    # Sort samples by name for convenience
    for f in vcf*.vcf.gz; do bcftools query -l \$f; done | sort > samples.txt

    bcftools merge -F x vcf*.vcf.gz | bcftools view -S samples.txt - | \
        bcftools view -e 'AN=0' - -o joined_vs_ref.merged0.vcf.gz
    bcftools index joined_vs_ref.merged0.vcf.gz
fi

## Annotate the merged file
cat >annot.hdr <<EOF
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">
##INFO=<ID=NOTE,Number=1,Type=String,Description="Region manual note">
EOF
bcftools annotate -h annot.hdr -c CHROM,FROM,TO,GENE -a ${file(params.gene_annotation_bed)} -o joined_vs_ref.merged1.vcf.gz joined_vs_ref.merged0.vcf.gz
bcftools index joined_vs_ref.merged1.vcf.gz
""" + ( regions_annotation.name != 'NO_FILE' ?
	"""
    ## Fix potential bad regions from igv
    awk -F"\t" 'BEGIN{OFS="\t"}; NF<4{\$4="EmptyNote"}; {gsub(/ /,"_",\$4);print}' ${regions_annotation} > fixed.regions.bed
    bcftools annotate -a fixed.regions.bed -c CHROM,FROM,TO,NOTE -h annot.hdr joined_vs_ref.merged1.vcf.gz -o joined_vs_ref.merged.vcf.gz
    bcftools index  joined_vs_ref.merged.vcf.gz
        """ : """
    mv joined_vs_ref.merged1.vcf.gz joined_vs_ref.merged.vcf.gz
    mv joined_vs_ref.merged1.vcf.gz.csi joined_vs_ref.merged.vcf.gz.csi
        """ ) +
    """
## Filtering concordant variants
bcftools filter -m+ -s CONCORDANT -e 'AN=N_SAMPLES && N_ALT=1' joined_vs_ref.merged.vcf.gz -o joined_vs_ref.merged.discordant.vcf.gz
bcftools index joined_vs_ref.merged.discordant.vcf.gz
    """
}


process CONVERT_VCF_TO_TABLE {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)}", saveAs: { saveas }

    input:
    path(vcf)
    val(saveas)

    output:
    path("common_variants.tsv")
    
    script:
    """
echo -ne 'CHROM\tPOS\tREF\tALT\tFILTER\tQUAL\tINFO/DP\tGENE\tNOTE' > common_variants.tsv
bcftools query -l ${vcf} | while read s; do echo -ne "\t\$s"; done >> common_variants.tsv
echo >> common_variants.tsv
bcftools query -u -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%QUAL\t%INFO/DP\t%GENE\t%NOTE[\t%TGT]\n' ${vcf} >> common_variants.tsv
"""
}
