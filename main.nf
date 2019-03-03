#!/usr/bin/env nextflow

log.info "===================================================================="
log.info "                          Mutect                                    "
log.info "===================================================================="


Channel.fromPath(params.samples)
    .ifEmpty { exit 1, "samples file not found: ${params.samples}" }
    .splitCsv(sep: '\t')
    .map{ patientId, sampleId, status, fastq1, fastq2 -> [patientId, sampleId, status, file(fastq1).baseName, [file(fastq1), file(fastq2)]] }
    .into { samples; reads }

// to download reads files from SRA 
int threads = Runtime.getRuntime().availableProcessors()

params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
if (params.fasta) {
    Channel.fromPath(params.fasta)
           .ifEmpty { exit 1, "fasta annotation file not found: ${params.fasta}" }
           .into { fasta_bwa; fasta_baserecalibrator; fasta_mutect }
}
params.fai = params.genome ? params.genomes[ params.genome ].fai ?: false : false
if (params.fai) {
    Channel.fromPath(params.fai)
           .ifEmpty { exit 1, "fasta index file not found: ${params.fai}" }
           .into { fai_mutect; fai_baserecalibrator }
}
params.dict = params.genome ? params.genomes[ params.genome ].dict ?: false : false
if (params.dict) {
    Channel.fromPath(params.dict)
           .ifEmpty { exit 1, "dict annotation file not found: ${params.dict}" }
           .into { dict_mutect; dict_baserecalibrator }
}
params.dbsnp_gz = params.genome ? params.genomes[ params.genome ].dbsnp_gz ?: false : false
if (params.dbsnp_gz) {
    Channel.fromPath(params.dbsnp_gz)
           .ifEmpty { exit 1, "dbsnp annotation file not found: ${params.dbsnp_gz}" }
           .set { dbsnp_gz}
}
params.dbsnp_idx_gz = params.genome ? params.genomes[ params.genome ].dbsnp_idx_gz ?: false : false
if (params.dbsnp_idx_gz) {
    Channel.fromPath(params.dbsnp_idx_gz)
           .ifEmpty { exit 1, "dbsnp_idx_gz annotation file not found: ${params.dbsnp_idx_gz}" }
           .set { dbsnp_idx_gz}
}
params.golden_indel_gz = params.genome ? params.genomes[ params.genome ].golden_indel_gz ?: false : false
if (params.golden_indel_gz) {
    Channel.fromPath(params.golden_indel_gz)
           .ifEmpty { exit 1, "golden_indel_gz annotation file not found: ${params.golden_indel_gz}" }
           .set { golden_indel_gz }
}
params.golden_indel_idx_gz = params.genome ? params.genomes[ params.genome ].golden_indel_idx_gz ?: false : false
if (params.golden_indel_idx_gz) {
    Channel.fromPath(params.golden_indel_idx_gz)
           .ifEmpty { exit 1, "golden_indel_idx_gz annotation file not found: ${params.golden_indel_idx_gz}" }
           .set { golden_indel_idx_gz }
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

process gunzip_dbsnp {
    tag "$dbsnp_gz"

    input:
    file dbsnp_gz from dbsnp_gz
    file dbsnp_idx_gz from dbsnp_idx_gz

  	output:
  	file "*.vcf" into dbsnp, dbsnp_variantrecalibrator_snps, dbsnp_variantrecalibrator_indels
  	file "*.vcf.idx" into dbsnp_idx, dbsnp_idx_variantrecalibrator_snps, dbsnp_idx_variantrecalibrator_indels

    script:
    if ( "${dbsnp_gz}".endsWith(".gz") ) {
     """
     gunzip -d --force $dbsnp_gz
   	 gunzip -d --force $dbsnp_idx_gz
     """
   } else {
     """
     cp $dbsnp_gz dbsnp.vcf
     cp $dbsnp_idx_gz dbsnp.vcf.idx
     """
   }
}

process gunzip_golden_indel {
    tag "$golden_indel_gz"

    input:
    file golden_indel_gz from golden_indel_gz
    file golden_indel_idx_gz from golden_indel_idx_gz

    output:
    file "*.vcf" into golden_indel, golden_indel_variantrecalibrator_indels
    file "*.vcf.idx" into golden_indel_idx, golden_indel_idx_variantrecalibrator_indels

    script:
    if ( "${golden_indel_gz}".endsWith(".gz") ) {
        """
        gunzip -d --force $golden_indel_gz
        gunzip -d --force $golden_indel_idx_gz
        """
    } else {
        """
        cp $golden_indel_gz golden_indel.vcf
        cp $golden_indel_idx_gz golden_indel.vcf.idx
        """
    }
}

process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}
    container 'lifebitai/fastqc'

    when:
    !params.skip_fastqc

    input:
    set val(patientId), val(sampleId), val(status), val(name), file(reads) from reads

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}

bwa_index = fasta_bwa.merge(bwa_index_amb, bwa_index_ann, bwa_index_bwt, bwa_index_pac, bwa_index_sa)
bwa = samples.combine(bwa_index)

process BWA {
    tag "$name"
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
    set val(patientId), val(sampleId), val(status), val(name), file("${name}-sorted.bam") into bam_sort, bam_sort_qc

    """
    samtools sort -o ${name}-sorted.bam -O BAM $sam
    """
}

process RunBamQCmapped {
    tag "$bam"

    container 'maxulysse/sarek:latest'

    input:
    set val(patientId), val(sampleId), val(status), val(name), file(bam) from bam_sort_qc

    output:
    file("${name}") into bamQCmappedReport

    when: !params.skip_multiqc

    script:
    // TODO: add --java-mem-size=${task.memory.toGiga()}G
    """
    qualimap \
    bamqc \
    -bam ${bam} \
    --paint-chromosome-limits \
    --genome-gc-distr HUMAN \
    -nt ${task.cpus} \
    -skip-duplicated \
    --skip-dup-mode 0 \
    -outdir ${name} \
    -outformat HTML
    """
}

process MarkDuplicates {
    tag "$bam_sort"
    container 'broadinstitute/gatk:latest'

    input:
    set val(patientId), val(sampleId), val(status), val(name), file(bam_sort) from bam_sort

    output:
    set val(name), file("${name}.bam"), file("${name}.bai"), val(patientId), val(sampleId), val(status) into bam_markdup_baserecalibrator, bam_markdup_applybqsr
    file ("${name}.bam.metrics") into markDuplicatesReport

    """
    gatk MarkDuplicates \
    -I ${bam_sort} \
    --CREATE_INDEX true \
    --M ${name}.bam.metrics \
    -O ${name}.bam
    """
}

baserecalibrator_index = fasta_baserecalibrator.merge(fai_baserecalibrator, dict_baserecalibrator, dbsnp, dbsnp_idx, golden_indel, golden_indel_idx)
baserecalibrator = bam_markdup_baserecalibrator.combine(baserecalibrator_index)

process BaseRecalibrator {
    tag "$bam_markdup"
    container 'broadinstitute/gatk:latest'

    input:
    set val(name), file(bam_markdup), file(bai), val(patientId), val(sampleId), val(status), 
    file(fasta), file(fai), file(dict), file(dbsnp), file(dbsnp_idx), file(golden_indel), file(golden_indel_idx) from baserecalibrator

    output:
    set val(name), file("${name}_recal_data.table") into baserecalibrator_table

    """
    gatk BaseRecalibrator \
    -I $bam_markdup \
    --known-sites $dbsnp \
    --known-sites $golden_indel \
    -O ${name}_recal_data.table \
    -R $fasta
    """
}

applybqsr = baserecalibrator_table.join(bam_markdup_applybqsr)

process ApplyBQSR {
    tag "$baserecalibrator_table"
    container 'broadinstitute/gatk:latest'

    input:
    set val(name), file(baserecalibrator_table), file(bam), file(bai), val(patientId), val(sampleId), val(status) from applybqsr

    output:
    set val(patientId), val(sampleId), val(status), val(name), file("${name}_bqsr.bam"), file("${name}_bqsr.bai") into bam_mutect

    script:
    """
    gatk ApplyBQSR -I $bam -bqsr $baserecalibrator_table -OBI -O ${name}_bqsr.bam
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

process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'
    container 'ewels/multiqc:v1.7'

    when:
    !params.skip_multiqc

    input:
    file (fastqc:'fastqc/*') from fastqc_results.collect().ifEmpty([])
    file (bam_metrics) from markDuplicatesReport.collect().ifEmpty([])
    file (bamQC) from bamQCmappedReport.collect().ifEmpty([])

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    """
    multiqc . -m fastqc -m picard -m qualimap
    """
}