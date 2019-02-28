#!/usr/bin/env nextflow

log.info "===================================================================="
log.info "                          Mutect                                    "
log.info "===================================================================="


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


/*
 * Create a channel for input read files
 * Dump can be used for debugging purposes, e.g. using the -dump-channels operator on run
 */
if (params.accession) {
    Channel.fromPath(params.accession)
    .ifEmpty { exit 1, "Text file containing SRA id's not found: ${params.accession}" }
    .set { sraIDs }
} else if (params.reads) {
  Channel
      .fromPath(params.reads)
      .map { file -> tuple(file.baseName, file) }
      .ifEmpty { exit 1, "Cannot find any reads matching: ${reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
      .dump(tag:'input')
      .into { reads_samplename; reads_bwa }
  } else if (params.reads_folder){
    reads="${params.reads_folder}/${params.reads_prefix}_{1,2}.${params.reads_extension}"
    Channel
        .fromFilePairs(reads, size: 2)
        .ifEmpty { exit 1, "Cannot find any reads matching: ${reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .dump(tag:'input')
        .into { reads_samplename; reads_bwa }
  } else if (params.bam) {
  Channel.fromPath(params.bam)
         .map { file -> tuple(file.baseName, file) }
         .ifEmpty { exit 1, "BAM file not found: ${params.bam}" }
         .set { bam_bqsr }
} else {
  exit 1, "Please specify either --reads singleEnd.fastq, --reads_folder pairedReads, --bam myfile.bam or --accession SRAids.txt"
}


if (params.accession) {
    sraIDs.splitText().map { it -> it.trim() }.set { singleSRAId }

    process fastqDump {
        tag "$id"
        container = 'lifebitai/kallisto-sra'
        cpus threads

        input:
        val id from singleSRAId

        output:
        set val(id), file('*.fastq.gz') into reads_bwa

        script:
        """
        parallel-fastq-dump --sra-id $id --threads ${task.cpus} --gzip
        """ 
    }
}


bwa_index = fasta_bwa.merge(bwa_index_amb, bwa_index_ann, bwa_index_bwt, bwa_index_pac, bwa_index_sa)
bwa = reads_bwa.combine(bwa_index)

process BWA {
    tag "$reads"
    container 'kathrinklee/bwa:latest'

    input:
    set val(name), file(reads), file(fasta), file(amb), file(ann), file(bwt), file(pac), file(sa) from bwa

    output:
    set val(name), file("${name}.sam") into sam

    """
    bwa mem -M -R '@RG\\tID:${name}\\tSM:${name}\\tPL:Illumina' $fasta $reads > ${name}.sam
    """
}


process BWA_sort {
    tag "$sam"
    container 'lifebitai/samtools:latest'

    input:
    set val(name), file(sam) from sam

    output:
    set val(name), file("${name}-sorted.bam") into bam_sort

    """
    samtools sort -o ${name}-sorted.bam -O BAM $sam
    """
}

process MarkDuplicates {
    tag "$bam_sort"
    container 'broadinstitute/gatk:latest'

    input:
    set val(name), file(bam_sort) from bam_sort

    output:
    set val(name), file("${name}.bam") into bam_index

    """
    gatk MarkDuplicates -I $bam_sort -M metrics.txt -O ${name}.bam
    """
}

process IndexBam {
  tag "$bam"
  container 'lifebitai/samtools:latest'

  input:
  set val(name), file(bam) from bam_index

  output:
  set val(name), file("${name}.bam"), file("${name}.bam.bai") into bam_mutect

  script:
  """
  cp $bam bam.bam
  mv bam.bam ${name}bam
  samtools index ${name}.bam
  """
}


ref_mutect = fasta_mutect.merge(fai_mutect, dict_mutect)
mutect = bam_mutect.combine(ref_mutect)

process Mutect2 {
    tag "$bam"
    container 'broadinstitute/gatk:latest'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    set val(name), file(bam), file(bai), file(fasta), file(fai), file(dict) from mutect

    output:
    file("${name}.vcf") into results

    script:
    """
    gatk Mutect2 \
    -I $bam \
    -O ${name}.vcf \
    -R $fasta \
    -tumor $name 
    """
}

