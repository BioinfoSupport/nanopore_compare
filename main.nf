// -*- groovy -*-

// Required parameter
params.sample_files = []
// Old style, it is prefered to use sample_files
params.sample_fastqs = []
// params.sample_files = ["data/nanopore/barcode01.guppy.pass.fastq.gz","data/nanopore/barcode02.guppy.pass.fastq.gz"]

// genome size is not really needed
// no-lat-contigs to prevent extra contigs?
params.flye_flags = "--genome-size 6m --asm-coverage 50 --no-alt-contigs"
params.medaka_model = "r1041_e82_400bps_sup_g615"
params.medaka_variant_model = "r1041_e82_400bps_sup_variant_g615"
// https://github.com/nanoporetech/medaka/commit/11e9614c4c142a1f2da2b6735eca4911f787b399
// FASTQ dorada model `dna_r10.4.1_e8.2_sup@v3.5.1`
// corresponds to `r1041_e82_400bps_sup_g615 model`

params.exclude_samples = []
params.consensus_ref_sample = ""

params.ref_fa = "data/saureus/Saureus8325.fasta"
params.ref_id = "SA8325"
params.ref_gff = "data/saureus/Saureus8325.gff"
params.gene_annotation_bed = "${workflow.projectDir}/assets/NO_FILE"

params.store_dir = "store_dir"
params.output_dir = "out"

params.my_medaka_compare = false
params.medaka_polish = true

params.save_medaka_consensus = true

params.regions_annotation = "${workflow.projectDir}/assets/NO_FILE"
// Some marginal calls appear
params.medaka_lowq = "20"

params.prokka_args = "--usegenus --genus Staphylococcus --species aureus --strain 8325-4"

// CS tag is probably useless. '-y' is required to copy modification MM and ML tags
params.minimap2_cs_tag = "--cs -y"

params.do_flye_meta = false

params.publish_mode = "symlink"

params.help = false


// Optional additional callers
params.clair3 = false
params.clair3_platform = "ont"
params.clair3_model = "r1041_e82_400bps_sup_g615"

params.sniffles = false



import groovy.json.JsonOutput


include {
    GET_MEDAKA_MODEL;
    GET_MEDAKA_MODEL as GET_MEDAKA_VARIANT_MODEL;
    GET_CLAIR3_MODEL;
    CRAM_TO_FASTQ_GZ
    DIR_TO_FASTQ_GZ
} from "./modules/utils.nf"

include {
    MEDAKA_CALL as MEDAKA_CALL_VS_REF;
    MEDAKA_CALL as MEDAKA_CALL_VS_CONSENSUS
    MM2_CALL as MM2_CALL_VS_REF;
    MM2_CALL as MM2_CALL_VS_CONSENSUS
    CLAIR3_CALL
    CLAIR3_CALL_SENSITIVE
    SNIFFLES_CALL
    REMAP_CONSENSUS_TO_REFERENCE;
    REMAP_CONSENSUS_TO_REFERENCE as REMAP_CONSENSUS_TO_REFERENCE_CONS
    MM2_CALL_ANNOTATE;
    MM2_CALL_ANNOTATE as MM2_CALL_ANNOTATE_CONS;
    MERGE_MEDAKA_AND_CONS_CALLS;
    MERGE_MEDAKA_AND_CONS_CALLS as MERGE_MEDAKA_AND_CONS_CALLS_CONS_NOLIFT
    MERGE_MEDAKA_AND_CONS_CALLS as MERGE_MEDAKA_AND_CONS_CALLS_CONS_LIFTED
    MERGE_JOINED_CALL_FILES;
    MERGE_JOINED_CALL_FILES as MERGE_JOINED_CALL_FILES_CONS;
    MERGE_JOINED_CALL_FILES as MERGE_JOINED_CALL_FILES_CONS_LIFTED;
    CONVERT_VCF_TO_TABLE;
    CONVERT_VCF_TO_TABLE as CONVERT_VCF_TO_TABLE_CONS} from  './modules/calls.nf'

include {
    MAP_READS as MAP_READS_TO_REF
    MAP_READS as MAP_READS_TO_CONS
    BAMSTATS_MAPPED_READS as BAMSTATS_MAPPED_READS_TO_REF
    FLAGSTAT_MAPPED_READS
} from './modules/mapping.nf'

include {
    PROKKA_ANNOTATE
    LIFTOVER_REF_ANNOTATIONS 
    ADD_GENOME_JSON
    GFF_TO_GENE_BED
} from './modules/annotation.nf'

////////////////////////////////////////////////////////////////////////
// Consensus generation
////////////////////////////////////////////////////////////////////////

// Assembles draft consensus assembly
process FLYE_ASSEMBLE {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)/meta.name}"
//        pattern: "flye/{assembly_info.txt,assembly_graph.gfa}",
//        saveAs: { it.replaceFirst(/\//, ".") }
    
    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("flye/assembly.fasta")
    tuple val(meta), path("flye/*")

    script:
    """
flye ${params.flye_flags} --threads ${task.cpus} \
    --out-dir flye --nano-hq ${fastq}
    """
    stub:
    """
mkdir flye
touch flye/assembly.fasta flye/assembly_info.txt flye/assembly_graph.gfa
    """
}

process FLYE_ASSEMBLE_FULL {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)/meta.name}",
        saveAs: { it.replaceFirst(/flye\//, "flye_full/") }
    
    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("flye/assembly.fasta")
    tuple val(meta), path("flye/*")

    script:
    """
flye --meta --keep-haplotypes --threads ${task.cpus} \
    --out-dir flye --nano-hq ${fastq}
    """
    stub:
    """
mkdir flye
touch flye/assembly.fasta flye/assembly_info.txt flye/assembly_graph.gfa
    """
}


// Polish the draft consensus assembly, and save it as fastq, so qualities are avilable
process MEDAKA_POLISH {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)/meta.name}",
             pattern: "medaka/consensus.fastq",
             saveAs: { "consensus.fastq" },
             enabled: params.save_medaka_consensus
    
    input:
    tuple val(meta), path(draft_fasta), path(fastq)
    path(medaka_model_path)

    output:
    tuple val(meta), path("medaka/consensus.fastq")
    tuple val(meta), path("medaka/*")

    script:
    """
medaka_consensus -q -i ${fastq} -d ${draft_fasta} -o medaka \
    -t ${task.cpus} -m ${medaka_model_path}
    """

    stub:
    """
mkdir medaka
touch medaka/consensus.fastq
    """
}


// Utility conversion
// Adds the minimap index as it will be used actively
// TODO Potentially, stage as hard links instead using stageInMode for medaka/minimap processes?
process FASTQ_TO_FASTA_MMI {
    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("seq.fasta"), path("seq.fasta.fai"), path("seq.fasta.map-ont.mmi")

    script:
    """
sed -n '1~4s/^@/>/p;2~4p' ${fastq} > seq.fasta
samtools faidx seq.fasta
# Using medaka mini_align defaults -- move to options?
minimap2 -I 16G -x map-ont -d seq.fasta.map-ont.mmi seq.fasta
    """
}

// Add mmi to the original reference
// Quite hacky
process FASTA_ADD_MMI {
    storeDir "${file(params.ref_fa).getParent()}"
    
    input:
    tuple val(meta), path(ref), path(fai)

    output:
    tuple val(meta), path(ref), path(fai), path("${mmi}")

    script:
    mmi = ref+".map-ont.mmi"
    """
# Using medaka mini_align defaults -- move to options?
minimap2 -I 16G -x map-ont -d ${mmi} ${ref}
    """
}

process FASTA_ADD_FAI {
    storeDir "${file(params.ref_fa).getParent()}"
    
    input:
    tuple val(meta), path(ref)

    output:
    tuple val(meta), path(ref), path(fai)
    
    script:
    fai = ref+".fai"
    """
samtools faidx ${ref}
    """
}

// Generate the polished consensus and return it
workflow POLISHED_CONSENSUS {
    take:
    fastq
    medaka_model_path

    main:

    draft = FLYE_ASSEMBLE(fastq)
    if (params.do_flye_meta) {
        FLYE_ASSEMBLE_FULL(fastq)
    }
    if (params.medaka_polish) {
    // Have to join matching on meta to have draft and fastq in sync
        polished = MEDAKA_POLISH(draft[0].join(fastq), medaka_model_path)
        result = polished[0]
    } else {
        result = draft[0]
    }

    emit:
    result
}

////////////////////////////////////////////////////////////////////////





workflow CALL_VS_REF {
    take:
    reference_fasta
    cons_fastq
    fastq
    medaka_model_path

    main:

    reference_fasta_nommi = reference_fasta.map { it.subList(0,3) }
    
    vcfs = MEDAKA_CALL_VS_REF(reference_fasta, fastq, medaka_model_path)
    vcfs_mm2 = MM2_CALL_VS_REF(reference_fasta_nommi, cons_fastq)

    bam = REMAP_CONSENSUS_TO_REFERENCE(reference_fasta_nommi, cons_fastq)

    vcfs_mm2_annotated = MM2_CALL_ANNOTATE(
        reference_fasta_nommi,
        vcfs_mm2[0].map({[it[0],it[1]]})
            .join(
        bam.map({[it[0],it[1]]})))
    // TODO: Add vcfs based on consensus seq
    //       Annotate by conensus base quality
    // TODO: Merge properly

    joined_vcfs = MERGE_MEDAKA_AND_CONS_CALLS(reference_fasta_nommi, vcfs.join(vcfs_mm2_annotated), "")
    
    // Bad emit?
    emit:
    vcfs[0]
    vcfs_mm2_annotated
    joined_vcfs
}










////////////////////////////////////////////////////////////////////////
// Calling versus the consensus of one of the samples
////////////////////////////////////////////////////////////////////////

process ADD_ASSEMBLY_PAF_FOR_LIFTOVER {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)/cons_meta.name}",
        pattern: "*.paf",
        saveAs: { it.replaceFirst(/ref/,ref_meta.name) }

    input:
    tuple val(ref_meta), path(ref_fa), path(ref_fa_fai)
    tuple val(cons_meta), path(cons_fastx), path(cons_fastx_fai), path(cons_mmi)

    output:
    tuple val(cons_meta), path(cons_fastx), path(cons_fastx_fai), path(cons_mmi), path("consensus_vs_ref.mm2.paf")

    script:
    """
minimap2 -un -cx asm5 --cs ${ref_fa} ${cons_fastx} > consensus_vs_ref.mm2.paf
    """
}

process ADD_ASSEMBLY_CHAIN_FOR_LIFTOVER {
    publishDir mode: 'copy', path: "${file(params.output_dir)/cons_meta.name}",
        saveAs: { it.replaceFirst(/ref/,ref_meta.name) }

    input:
    tuple val(ref_meta), path(ref_fa), path(ref_fa_fai)
    tuple val(cons_meta), path(cons_fastx), path(cons_fastx_fai), path(cons_mmi)

    output:
    tuple val(cons_meta), path("consensus.fasta"), path("consensus.fasta.fai"), path("consensus.fasta.map-ont.mmi"), path("ref_vs_consensus.chain")

    script:
    """
minimap2 -un -cx asm5 ${cons_fastx} ${ref_fa} > ref_vs_consensus.mm2.paf
ln ${cons_fastx} consensus.fasta
ln ${cons_fastx_fai} consensus.fasta.fai
ln ${cons_mmi} consensus.fasta.map-ont.mmi
paf2chain --input ref_vs_consensus.mm2.paf > ref_vs_consensus.chain
    """
}

// Problem: --no-comp-alleles produces bad result in case of "restoring" SNP. What is teh way arond this?
process LIFTOVER_MEDAKA_CALLS {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)/meta.name}",
        saveAs: { it.replaceFirst(/ref/, paf_meta.name) }
    
    input:
    tuple val(paf_meta), path(fastx), path(fai), path(mmi), path(chain)
    tuple val(ref_meta), path(ref_fasta), path(_ref_fasta_fai)
    tuple val(meta), path(vcf), path(vcf_csi)

    output:
    tuple val(meta), path("medaka.reads_vs_ref.lifted.vcf.gz"), path("medaka.reads_vs_ref.lifted.vcf.gz.csi"), path("medaka.reads_vs_ref.lifted.vcf.unmap")
    
    script:
    """
CrossMap vcf --no-comp-alleles ${chain} ${vcf} ${ref_fasta} medaka.reads_vs_ref.lifted.vcf
bcftools sort medaka.reads_vs_ref.lifted.vcf -o medaka.reads_vs_ref.lifted.vcf.gz
bcftools index medaka.reads_vs_ref.lifted.vcf.gz
    """
}

process LIFTOVER_MM2_CALLS {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)/meta.name}",
        saveAs: { it.replaceFirst(/ref/, paf_meta.name) }

    input:
    tuple val(paf_meta), path(fastx), path(fai), path(mmi), path(chain)
    tuple val(ref_meta), path(ref_fasta), path(ref_fasta_fai)
    tuple val(meta), path(vcf), path(vcf_csi)

    output:
    tuple val(meta), path("mm2.cons_vs_ref.lifted.vcf.gz"), path("mm2.cons_vs_ref.lifted.vcf.gz.csi"), path("mm2.cons_vs_ref.lifted.vcf.unmap")
    
    script:
    """
CrossMap vcf --no-comp-alleles ${chain} ${vcf} ${ref_fasta} mm2.cons_vs_ref.lifted.vcf
bcftools sort mm2.cons_vs_ref.lifted.vcf -o mm2.cons_vs_ref.lifted.vcf.gz
bcftools index mm2.cons_vs_ref.lifted.vcf.gz
    """
}


// TODO Report missing coordinates
workflow CALL_VS_CONS {
    take:
    consensus_reference_paf_and_fasta
    cons_fastq
    fastq
    medaka_model_path
    ref
    
    main:

    // Remove paf for those who does not need
    reference_fasta_mmi = consensus_reference_paf_and_fasta.map({[it[0], it[1], it[2], it[3]]})
    // and also remove the mmi index
    reference_fasta = consensus_reference_paf_and_fasta.map({[it[0], it[1], it[2]]})
    
    vcfs = MEDAKA_CALL_VS_CONSENSUS(reference_fasta_mmi, fastq, medaka_model_path)
    vcfs_mm2 = MM2_CALL_VS_CONSENSUS(reference_fasta, cons_fastq)

    // Danger step -- can lose positions
    bam = REMAP_CONSENSUS_TO_REFERENCE_CONS(reference_fasta, cons_fastq)

    // Double danger -- losing positions
    vcfs_mm2_annotated = MM2_CALL_ANNOTATE_CONS(
        reference_fasta,
        vcfs_mm2[0].map({[it[0],it[1]]})
            .join(
        bam.map({[it[0],it[1]]})))


    // Before liftover
    joined_vcfs_nolift = MERGE_MEDAKA_AND_CONS_CALLS_CONS_NOLIFT(reference_fasta,
         vcfs[0].join(vcfs_mm2_annotated[0]), "")
    // Lifted files
    lift_vcfs = LIFTOVER_MEDAKA_CALLS(consensus_reference_paf_and_fasta, ref, vcfs[0]) // | map( {it.subList(0,3)} )
    lift_vcfs_mm2 = LIFTOVER_MM2_CALLS(consensus_reference_paf_and_fasta, ref, vcfs_mm2_annotated[0]) // | map( {it.subList(0,3)} )

    joined_vcfs = MERGE_MEDAKA_AND_CONS_CALLS_CONS_LIFTED(reference_fasta,
           lift_vcfs.map( {it.subList(0,3)} )
           .join(lift_vcfs_mm2.map( {it.subList(0,3)} )), ".lifted")

    emit:
    lift_vcfs[0]
    lift_vcfs_mm2[0]
    joined_vcfs
    joined_vcfs_nolift
}






workflow ANNOTATE_CONSENSUS {
    take:
    ref
    cons

    main:

    LIFTOVER_REF_ANNOTATIONS(ref, cons, file(params.ref_gff, checkIfExists:true))
    PROKKA_ANNOTATE(cons) | ADD_GENOME_JSON
    
//    emit:
}



////////////////////// Final VCF merging processes

process MERGE_VCFS {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)}",
        saveAs: {
        it.replaceFirst(/merged_vcf/,
        'medaka.reads_vs_'+params.consensus_ref_sample+'.lifted.merged') }

    input:
    path 'vcf*.vcf.gz'

    output:
    tuple path('merged_vcf.vcf.gz'), path('merged_vcf.vcf.gz.csi')

    script:
    """
for f in vcf*; do
    bcftools index \$f
done
if [ -f vcf.vcf.gz ]; then # Single file given, no merge needed
    cp vcf.vcf.gz merged_vcf.vcf.gz
    mv vcf.vcf.gz.csi merged_vcf.vcf.gz.csi
else
    bcftools merge vcf*.vcf.gz | bcftools filter -e 'QUAL<30' -s LOWQUAL | bcftools filter -e 'DP<10' -s 'LOWDEPTH' -o m.vcf.gz
    bcftools view -s `bcftools query -l m.vcf.gz|sort|tr '\n' ','|sed 's/,\$//'` m.vcf.gz -o merged_vcf.vcf.gz
    rm m.vcf.gz
    bcftools index merged_vcf.vcf.gz
fi
"""
}

process MERGE_READS_VS_REF_VCFS {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)}",
        saveAs: { it.replaceFirst(/merged_reads_vs_ref_vcf/,
          'medaka.reads_vs_'+params.ref_id+'.merged') }

    input:
    path 'vcf*.vcf.gz'

    output:
    tuple path('merged_reads_vs_ref_vcf.vcf.gz'), path('merged_reads_vs_ref_vcf.vcf.gz.csi')

    script:
    """
for f in vcf*; do
    bcftools index \$f
done
if [ -f vcf.vcf.gz ]; then # Single file given, no merge needed
    cp vcf.vcf.gz merged_reads_vs_ref_vcf.vcf.gz
    mv vcf.vcf.gz.csi merged_reads_vs_ref_vcf.vcf.gz.csi
else
    bcftools merge vcf*.vcf.gz | bcftools filter -e 'QUAL<30' -s LOWQUAL | bcftools filter -e 'DP<10' -s 'LOWDEPTH' -o m.vcf.gz
    bcftools view -s `bcftools query -l m.vcf.gz|sort|tr '\n' ','|sed 's/,\$//'` m.vcf.gz -o merged_reads_vs_ref_vcf.vcf.gz
    rm m.vcf.gz
    bcftools index merged_reads_vs_ref_vcf.vcf.gz
fi
"""
}

process MERGE_READS_VS_REF_VCFS_MM {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)}",
        saveAs: { it.replaceFirst(/merged_mm_vs_ref_vcf/,
          'mm2.cons_vs_'+params.ref_id+'.merged') }

    input:
    path 'vcf*.vcf.gz'

    output:
    tuple path('merged_mm_vs_ref_vcf.vcf.gz'), path('merged_mm_vs_ref_vcf.vcf.gz.csi')

    script:
    """
for f in vcf*; do
    bcftools index \$f
done
if [ -f vcf.vcf.gz ]; then # Single file given, no merge needed
    cp vcf.vcf.gz merged_mm_vs_ref_vcf.vcf.gz
    mv vcf.vcf.gz.csi merged_mm_vs_ref_vcf.vcf.gz.csi
else
    bcftools merge vcf*.vcf.gz | bcftools filter -e 'QUAL<30' -s LOWQUAL | bgzip > m.vcf.gz
    bcftools view -s `bcftools query -l m.vcf.gz|sort|tr '\n' ','|sed 's/,\$//'` m.vcf.gz -o merged_mm_vs_ref_vcf.vcf.gz
    rm m.vcf.gz
    bcftools index merged_mm_vs_ref_vcf.vcf.gz
fi
"""
}

process MERGE_VCFS_MM {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)}",
        saveAs: {
        it.replaceFirst(/merged_mm_vcf/,
        'mm2.cons_vs_'+params.consensus_ref_sample+'.lifted.merged') }

    input:
    path 'vcf*.vcf.gz'

    output:
    tuple path('merged_mm_vcf.vcf.gz'), path('merged_mm_vcf.vcf.gz.csi')

    script:
    """
for f in vcf*; do
    bcftools index \$f
done
if [ -f vcf.vcf.gz ]; then # Single file given, no merge needed
    cp vcf.vcf.gz merged_mm_vcf.vcf.gz
    mv vcf.vcf.gz.csi merged_mm_vcf.vcf.gz.csi
else
    bcftools merge vcf*.vcf.gz | bcftools filter -e 'QUAL<30' -s LOWQUAL | bgzip > m.vcf.gz
    bcftools view -s `bcftools query -l m.vcf.gz|sort|tr '\n' ','|sed 's/,\$//'` m.vcf.gz -o merged_mm_vcf.vcf.gz
    rm m.vcf.gz
    bcftools index merged_mm_vcf.vcf.gz
fi
"""
}

process MERGE_CLAIR3_VCFS {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)}",
        saveAs: {it.replaceFirst(/ref/, params.ref_id)}

    input:
    path 'vcf*.vcf.gz'

    output:
    tuple path('clair_vs_ref.merged.vcf.gz'), path('clair_vs_ref.merged.vcf.gz.csi')

    script:
    """
for f in vcf*.vcf.gz; do
    bcftools index \$f
done
if [ -f vcf.vcf.gz ]; then # Single file given, no merge needed
    cp vcf.vcf.gz clair_vs_ref.merged.vcf.gz
    mv vcf.vcf.gz.csi clair_vs_ref.merged.vcf.gz.csi
else
    bcftools merge vcf*.vcf.gz | bgzip > m.vcf.gz
    bcftools view -s `bcftools query -l m.vcf.gz|sort|tr '\n' ','|sed 's/,\$//'` m.vcf.gz -o clair_vs_ref.merged.vcf.gz
    rm m.vcf.gz
    bcftools index clair_vs_ref.merged.vcf.gz
fi
"""
}

process MERGE_CLAIR3_SENSITIVE_VCFS {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)}",
        saveAs: {it.replaceFirst(/ref/, params.ref_id).replaceFirst(/clair/, "clair_sensitive_")}

    input:
    path 'vcf*.vcf.gz'

    output:
    tuple path('clair_vs_ref.merged.vcf.gz'), path('clair_vs_ref.merged.vcf.gz.csi')

    script:
    """
for f in vcf*.vcf.gz; do
    bcftools index \$f
done
if [ -f vcf.vcf.gz ]; then # Single file given, no merge needed
    cp vcf.vcf.gz clair_vs_ref.merged.vcf.gz
    mv vcf.vcf.gz.csi clair_vs_ref.merged.vcf.gz.csi
else
    bcftools merge vcf*.vcf.gz | bgzip > m.vcf.gz
    bcftools view -s `bcftools query -l m.vcf.gz|sort|tr '\n' ','|sed 's/,\$//'` m.vcf.gz -o clair_vs_ref.merged.vcf.gz
    rm m.vcf.gz
    bcftools index clair_vs_ref.merged.vcf.gz
fi
"""
}






process MEDAKA_GENERATE_CONSENSUS_HDF {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)/meta.name}",
             saveAs: { it.replaceFirst(/ref/, ref_meta.name) }
    cpus 2 // Check -- this is really teh max here? Do not override in parameters later

    input:
    tuple val(ref_meta), path(ref_fa), path(ref_fa_fai), path(ref_fa_mmi)
    tuple val(meta), path(bam), path(bam_bai)
    path(medaka_model_path)

    output:
    tuple val(meta), path("${meta.name}_vs_ref.hdf")
    
    script:
    """
medaka consensus "${bam}" "${meta.name}_vs_ref.hdf" \
        --model "${medaka_model_path}" --batch_size "100" --threads "${task.cpus}"
    """
}
process MEDAKA_REPOLISH {
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)}"

input:
    tuple val(ref_meta), path(ref_fa), path(ref_fa_fai), path(ref_fa_mmi)
    tuple val(ref_consensus_meta), path(ref_hdf)

output:
    tuple val(ref_meta), path(polished_fasta), path("${polished_fasta}.fai")    

script:
    polished_fasta = "${ref_meta.name}_repolished.fasta"
    """
medaka stitch --threads ${task.cpus} ${ref_hdf} ${ref_fa} ${polished_fasta}
samtools faidx ${polished_fasta}
    """    
}
// process MEDAKA_VARIANTS {
//     publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)}",
//              saveAs: { it.replaceFirst(/ref/, ref_meta.name) }
//     cpus 2

//    input:
//     tuple val(ref_meta), path(ref_fa), path(ref_fa_fai), path(ref_fa_mmi)
//     path "consensus_hdf_*.hdf"

//     output:
//     path "medaka_variants_vs_ref.vcf"
    
//     script:
//     """
//     medaka variant --verbose --ambig_ref "${ref_fa}" consensus_hdf_*.hdf "medaka_variants_vs_ref.vcf"
//     """
// }
process MY_MEDAKA_COMPARE {
    cpus 1
    container "/home/users/b/bezrukov/medaka_compare/medaka_compare.sif"
    publishDir mode: "${params.publish_mode}", path: "${file(params.output_dir)/meta.name}",
             saveAs: { it.replaceFirst(/ref/, ref_meta.name) }

    input:
    tuple val(ref_meta), path(ref_fa), path(ref_fa_fai), path(ref_fa_mmi)
    tuple val(ref_consensus_meta), path(ref_hdf)
    tuple val(meta), path(sample_hdf)

    output:
    path "${meta.name}_medaka_compare_vs_ref.vcf"
    
    script:
    """
    medaka_compare ${ref_fa} ${ref_hdf} ${sample_hdf} > ${meta.name}_medaka_compare_vs_ref.vcf
    """
    
}


////////////////////////////////////////////////////////////////////////
////    Main workflow
////////////////////////////////////////////////////////////////////////
workflow {

    //// Normalize and check the provided parameters ///////////////////
    if (params.help) {
        println("""
Nanopore pipeline for microbial genome
""")
        exit(0)
    }
    // Combine new and old style sample file parameter
    sample_files = params.sample_files + params.sample_fastqs
    if (sample_files == []) {
        println("Please, provide --sample_files parameter")
        exit(1)
    }

    
    //// Input file preparation ////////////////////////////////////////
    channel.fromPath(sample_files, checkIfExists:true)
        .map({ // Generate sample names from sample filename
        meta = [name: it.getSimpleName()]
        [meta, it]
            })
        .filter({ !params.exclude_samples.contains(it[0].name) })
        .branch( {
        cram: it[1] ==~ ".*\\.(cram|bam)\$"
        fastq: it[1] ==~ ".*\\.fastq\\.gz\$"
        dir: file(it[1]).exists() && file(it[1]).isDirectory()
        unknown: true
            } )
        .set { sample_files_split_by_type }

    sample_files_split_by_type.unknown.view { "Sample file ${it[1]} has unsupported extension" }
    
    sample_files_split_by_type.fastq
        .concat(sample_files_split_by_type.cram | CRAM_TO_FASTQ_GZ)
        .concat(sample_files_split_by_type.dir  | DIR_TO_FASTQ_GZ )
        .set { fastq } // "fastq" is now the channel for all input files in fastq with metadata


    //// Reference preparation /////////////////////////////////////////
    channel.value([[name: params.ref_id], file(params.ref_fa)]) |
        FASTA_ADD_FAI | set { ref }  // Value channel with name, fasta, and fai
    ref |
        FASTA_ADD_MMI | set { ref_mmi } // Value channel with name, fasta, fai, and mmi index

    
    //// Model download ////////////////////////////////////////////////
    medaka_model_path = GET_MEDAKA_MODEL(params.medaka_model) // | view { "medaka model ${it}" }
    medaka_variant_model_path = GET_MEDAKA_VARIANT_MODEL(params.medaka_variant_model) // | view { "medaka variant model ${it}" }
    if (params.clair3) {
        clair_model_path = GET_CLAIR3_MODEL(params.clair3_model) // | view { "clair model ${it}" }
    }


    //// Consensus /////////////////////////////////////////////////////
    POLISHED_CONSENSUS(fastq, medaka_model_path)
    // POLISHED_CONSENSUS.out.cons_bam | view

    //// Calling vs Reference //////////////////////////////////////////
    MAP_READS_TO_REF(ref_mmi, fastq)

    // QC statistics
    MAP_READS_TO_REF.out | BAMSTATS_MAPPED_READS_TO_REF
    MAP_READS_TO_REF.out | FLAGSTAT_MAPPED_READS

    CALL_VS_REF(ref_mmi, POLISHED_CONSENSUS.out, fastq, medaka_variant_model_path) // | view

    CALL_VS_REF.out[0].map({it[1]}).collect() | MERGE_READS_VS_REF_VCFS
    CALL_VS_REF.out[1].map({it[1]}).collect() | MERGE_READS_VS_REF_VCFS_MM

    // Setup bed files for annotation
    gene_annotation_bed_file = file(params.gene_annotation_bed, checkIfExists:true)
    if (gene_annotation_bed_file.name == 'NO_FILE') {
	gene_annotation_bed_file = GFF_TO_GENE_BED(file(params.ref_gff, checkIfExists:true))
    }
    regions_annotation_bed_file = file(params.regions_annotation, checkIfExists:true)

    
    MERGE_JOINED_CALL_FILES(CALL_VS_REF.out[2].map({it[1]}).collect(),
            gene_annotation_bed_file,
            regions_annotation_bed_file,
            ref, "")
    // Here we take the discordant vcf
    CONVERT_VCF_TO_TABLE(MERGE_JOINED_CALL_FILES.out[1].map({it[0]}), "common_variants.tsv")

    //// CLAIR3 consensus calling /////////////////////////////////////
    if (params.clair3) {
        CLAIR3_CALL(ref, MAP_READS_TO_REF.out, clair_model_path)
        CLAIR3_CALL.out.map({it[1]}).collect() | MERGE_CLAIR3_VCFS

        CLAIR3_CALL_SENSITIVE(ref, MAP_READS_TO_REF.out, clair_model_path)
        CLAIR3_CALL_SENSITIVE.out.map({it[1]}).collect() | MERGE_CLAIR3_SENSITIVE_VCFS
    }

    //// SNIFFLES consensus calling ///////////////////////////////
    if (params.sniffles) {
        SNIFFLES_CALL(ref, MAP_READS_TO_REF.out)
    }

    //// End: Calling vs Reference //////////////////////////////////////////

    if (params.consensus_ref_sample != "") {
    
        //// Calling vs selected sample consensus //////////////////////////
        cons_ref_fasta = POLISHED_CONSENSUS.out.first({
            it[0].name == params.consensus_ref_sample}) | 
            FASTQ_TO_FASTA_MMI
        cons_ref_paf_fasta = ADD_ASSEMBLY_CHAIN_FOR_LIFTOVER(ref, cons_ref_fasta)

        ANNOTATE_CONSENSUS(ref, cons_ref_fasta)

        MAP_READS_TO_CONS(cons_ref_fasta, fastq)

        // It is not fully clear if it is better to use variant or consensus medaka model here
        if (params.my_medaka_compare) {
            MEDAKA_GENERATE_CONSENSUS_HDF(cons_ref_fasta, MAP_READS_TO_CONS.out, medaka_variant_model_path)
            // MEDAKA_GENERATE_CONSENSUS_HDF.out | view
            ref_hdf = MEDAKA_GENERATE_CONSENSUS_HDF.out.first({it[0].name==params.consensus_ref_sample})
            other_hdf = MEDAKA_GENERATE_CONSENSUS_HDF.out.filter({it[0].name!=params.consensus_ref_sample})
            // ref_hdf | view
            // other_hdf|view
            MY_MEDAKA_COMPARE(cons_ref_fasta, ref_hdf, other_hdf)
        }

        
        other_fastq = fastq.filter({
            it[0].name != params.consensus_ref_sample})
        other_cons_fastq = POLISHED_CONSENSUS.out.filter({
            it[0].name != params.consensus_ref_sample})

        CALL_VS_CONS(cons_ref_paf_fasta, other_cons_fastq, other_fastq, medaka_variant_model_path, ref) // | view
    
        CALL_VS_CONS.out[0].map({it[1]}).collect() | MERGE_VCFS
        CALL_VS_CONS.out[1].map({it[1]}).collect() | MERGE_VCFS_MM

        MERGE_JOINED_CALL_FILES_CONS_LIFTED(CALL_VS_CONS.out[2].map({it[1]}).collect(),
            gene_annotation_bed_file,
            regions_annotation_bed_file,
            cons_ref_fasta.map({it.subList(0,3)}), ".lifted")
        MERGE_JOINED_CALL_FILES_CONS(CALL_VS_CONS.out[3].map({it[1]}).collect(),
	     gene_annotation_bed_file,
             regions_annotation_bed_file,
             cons_ref_fasta.map({it.subList(0,3)}), "")

        // Here we take the common vcf, as tehy are vs reference sample already
        CONVERT_VCF_TO_TABLE_CONS(MERGE_JOINED_CALL_FILES_CONS.out[0].map({it[0]}), "variants_vs_"+params.consensus_ref_sample+".tsv")

        //// End: Calling vs selected sample consensus //////////////////////////
    }
    
}

// Dump the workflow parameters and current directory git status
workflow.onComplete = {
    file = new File(params.output_dir, "workflow.log").newWriter()
    file << "Info:\n"
    file << "$workflow" << "\n\n"
    file << "git status:\n"
    file << "git status -s -b".execute().text
    file << "git hash:\n"
    file << "git rev-parse --short HEAD".execute().text
    file << "\nDONE!\n"
    file.close()
    
    file = new File(params.output_dir, "params.json").newWriter()
    json_str = JsonOutput.toJson(params)
    json_indented = JsonOutput.prettyPrint(json_str)
    file << json_indented << "\n"
    file.close()    
}


