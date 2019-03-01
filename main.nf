#!/usr/bin/env nextflow

log.info "===================================================================="
log.info "                          Mutect                                    "
log.info "===================================================================="


Channel.fromPath(params.samples)
    .ifEmpty { exit 1, "samples file not found: ${params.samples}" }
    .splitCsv(sep: '\t')
    .map{ patientId, sampleId, status, fastq1, fastq2 -> [patientId, sampleId, status, file(fastq1).baseName, [file(fastq1), file(fastq2)]] }
    .set { samples }

// to download reads files from SRA 
int threads = Runtime.getRuntime().availableProcessors()

params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
if (params.fasta) {
    Channel.fromPath(params.fasta)
           .ifEmpty { exit 1, "fasta annotation file not found: ${params.fasta}" }
           .into { fasta_bwa; fasta_mutect }
}
params.fai = params.genome ? params.genomes[ params.genome ].fai ?: false : false
if (params.fai) {
    Channel.fromPath(params.fai)
           .ifEmpty { exit 1, "fasta index file not found: ${params.fai}" }
           .set { fai_mutect }
}
params.dict = params.genome ? params.genomes[ params.genome ].dict ?: false : false
if (params.dict) {
    Channel.fromPath(params.dict)
           .ifEmpty { exit 1, "dict annotation file not found: ${params.dict}" }
           .set { dict_mutect }
}

params.bwa_index_amb = params.genome ? params.genomes[ params.genome ].bwa_index_amb ?: false : false
if (params.bwa_index_amb) {
    Channel.fromPath(params.bwa_index_amb)
           .ifEmpty { exit 1, "bwa_index_amb annotation file not found: ${params.bwa_index_amb}" }
           .set { bwa_index_amb }
}
params.bwa_index_ann = params.genome ? params.genomes[ params.genome ].bwa_index_ann ?: false : false
if (params.bwa_index_ann) {
    Channel.fromPath(params.bwa_index_ann)
           .ifEmpty { exit 1, "bwa_index_ann annotation file not found: ${params.bwa_index_ann}" }
           .set { bwa_index_ann }
}
params.bwa_index_bwt = params.genome ? params.genomes[ params.genome ].bwa_index_bwt ?: false : false
if (params.bwa_index_bwt) {
    Channel.fromPath(params.bwa_index_bwt)
           .ifEmpty { exit 1, "bwa_index_bwt annotation file not found: ${params.bwa_index_bwt}" }
           .set { bwa_index_bwt }
}
params.bwa_index_pac = params.genome ? params.genomes[ params.genome ].bwa_index_pac ?: false : false
if (params.bwa_index_pac) {
    Channel.fromPath(params.bwa_index_pac)
           .ifEmpty { exit 1, "bwa_index_pac annotation file not found: ${params.bwa_index_pac}" }
           .set { bwa_index_pac }
}
params.bwa_index_sa = params.genome ? params.genomes[ params.genome ].bwa_index_sa ?: false : false
if (params.bwa_index_sa) {
    Channel.fromPath(params.bwa_index_sa)
           .ifEmpty { exit 1, "bwa_index_sa annotation file not found: ${params.bwa_index_sa}" }
           .set { bwa_index_sa }
}

bwa_index = fasta_bwa.merge(bwa_index_amb, bwa_index_ann, bwa_index_bwt, bwa_index_pac, bwa_index_sa)
bwa = samples.combine(bwa_index)

process BWA {
    tag "$reads"
    container 'kathrinklee/bwa:latest'

    input:
    set val(patientId), val(sampleId), val(status), val(name), file(reads),
    file(fasta), file(amb), file(ann), file(bwt), file(pac), file(sa) from bwa

    output:
    set val(patientId), val(sampleId), val(status), val(name), file("${name}.sam") into sam

    """
    bwa mem -M -R '@RG\\tID:${name}\\tSM:${name}\\tPL:Illumina' $fasta $reads > ${name}.sam
    """
}

process BWA_sort {
    tag "$sam"
    container 'lifebitai/samtools:latest'

    input:
    set val(patientId), val(sampleId), val(status), val(name), file(sam) from sam

    output:
    set val(patientId), val(sampleId), val(status), val(name), file("${name}-sorted.bam") into bam_sort

    """
    samtools sort -o ${name}-sorted.bam -O BAM $sam
    """
}

process MarkDuplicates {
    tag "$bam_sort"
    container 'broadinstitute/gatk:latest'

    input:
    set val(patientId), val(sampleId), val(status), val(name), file(bam_sort) from bam_sort

    output:
    set val(patientId), val(sampleId), val(status), val(name), file("${name}.bam") into bam_index

    """
    gatk MarkDuplicates -I $bam_sort -M metrics.txt -O ${name}.bam
    """
}


process IndexBam {
  tag "$bam"
  container 'lifebitai/samtools:latest'

  input:
  set val(patientId), val(sampleId), val(status), val(name), file(bam) from bam_index

  output:
  set val(patientId), val(sampleId), val(status), val(name), file("${name}.bam"), file("${name}.bam.bai") into bam_mutect

  script:
  """
  cp $bam bam.bam
  mv bam.bam ${name}.bam
  samtools index ${name}.bam
  """
}

bamsNormal = Channel.create()
bamsTumour = Channel.create()

bam_mutect.choice( bamsTumour, bamsNormal ) { it[2] =~ 1 ? 0 : 1 }

combined_bam = bamsNormal.combine(bamsTumour, by: 0)

ref_mutect = fasta_mutect.merge(fai_mutect, dict_mutect)
mutect = combined_bam.combine(ref_mutect)


process Mutect2 {
    tag "$bam"
    container 'broadinstitute/gatk:latest'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    set val(patientId), val(sampleId), val(status), val(name), file(bam), file(bai),
    val(tumourSampleId), val(tumourStatus), val(tumourName), file(tumourBam), file(tumourBai),
    file(fasta), file(fai), file(dict) from mutect

    output:
    file("${tumourSampleId}_vs_${sampleId}.vcf") into results

    script:
    """
    gatk Mutect2 \
    -R ${fasta}\
    -I ${tumourBam}  -tumor ${tumourName} \
    -I ${bam} -normal ${name} \
    -O ${tumourSampleId}_vs_${sampleId}.vcf

    #gatk --java-options "-Xmx\${task.memory.toGiga()}g" \
    #-L \${intervalBed} \
    """
}