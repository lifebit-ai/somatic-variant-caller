#!/usr/bin/env nextflow

log.info "===================================================================="
log.info "                          Mutect                                    "
log.info "===================================================================="

// set threadmem equal to total memory divided by number of threads
int threads = Runtime.getRuntime().availableProcessors()
threadmem = (((Runtime.getRuntime().maxMemory() * 4) / threads) as nextflow.util.MemoryUnit)


Channel.fromPath(params.samples)
    .ifEmpty { exit 1, "samples file not found: ${params.samples}" }
    .splitCsv(sep: '\t')
    .map{ patientId, sampleId, status, fastq1, fastq2 -> [patientId, sampleId, status, file(fastq1).baseName, [file(fastq1), file(fastq2)]] }
    .into { samples; reads; reads_kraken }

params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
if (params.fasta) {
    Channel.fromPath(params.fasta)
           .ifEmpty { exit 1, "fasta annotation file not found: ${params.fasta}" }
           .into { fasta_bwa; fasta_baserecalibrator; fasta_haplotypecaller; fasta_mutect; fasta_variant_eval }
}
params.fai = params.genome ? params.genomes[ params.genome ].fai ?: false : false
if (params.fai) {
    Channel.fromPath(params.fai)
           .ifEmpty { exit 1, "fasta index file not found: ${params.fai}" }
           .into { fai_mutect; fai_baserecalibrator; fai_haplotypecaller; fai_variant_eval }
}
params.dict = params.genome ? params.genomes[ params.genome ].dict ?: false : false
if (params.dict) {
    Channel.fromPath(params.dict)
           .ifEmpty { exit 1, "dict annotation file not found: ${params.dict}" }
           .into { dict_interval; dict_mutect; dict_baserecalibrator; dict_haplotypecaller; dict_variant_eval }
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
if (params.bed) {
    Channel.fromPath(params.bed)
           .ifEmpty { exit 1, "BED file for to define regions not found: ${params.bed}" }
           .into { bed; bed_basename }

    bed_basename.map { file -> tuple(file.baseName, file) }.set{ bed_interval }
}

kraken_db = params.kraken_db

process BedToIntervalList {
    tag "$bed"
    container 'broadinstitute/gatk:latest'

    input:
    set val(name), file(bed) from bed_interval
    file dict from dict_interval

    output:
    file("${name}.interval_list") into interval_list

    when: params.bed

    script:
    """
    gatk BedToIntervalList \
    -I ${bed} \
    -O ${name}.interval_list \
    -SD ${dict}

    # remove header, columns 3 onwards & reformat
    sed -i '/^@/d' ${name}.interval_list
    cut -f 1-3 ${name}.interval_list > tmp.interval_list
    awk 'BEGIN { OFS = "" }{ print \$1,":",\$2,"-",\$3 }' tmp.interval_list > ${name}.interval_list
    """
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
    container 'flowcraft/fastqc:0.11.7-1'

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
    -M ${name}.bam.metrics \
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
    file("*data.table") into baseRecalibratorReport

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
    set val(patientId), val(sampleId), val(status), val(name), file("${name}_bqsr.bam"), file("${name}_bqsr.bai") into bam_for_qc, bam_haplotypecaller, bam_mutect

    script:
    """
    gatk ApplyBQSR -I $bam -bqsr $baserecalibrator_table -OBI -O ${name}_bqsr.bam
    """
}

process RunBamQCrecalibrated {
    tag "$bam"

    container 'maxulysse/sarek:latest'

    input:
    set val(patientId), val(sampleId), val(status), val(name), file(bam), file(bai) from bam_for_qc

    output:
    file("${name}_recalibrated") into bamQCrecalibratedReport

    when: !params.skip_multiqc

    script:
    // TODO: add --java-mem-size=${task.memory.toGiga()}G \
    """
    qualimap \
    bamqc \
    -bam ${bam} \
    --paint-chromosome-limits \
    --genome-gc-distr HUMAN \
    -nt ${task.cpus} \
    -skip-duplicated \
    --skip-dup-mode 0 \
    -outdir ${name}_recalibrated \
    -outformat HTML
    """
}


interval_list
    .splitText()
    .map { it -> it.trim() }
    .set { intervals }

haplotypecaller_index = fasta_haplotypecaller.merge(fai_haplotypecaller, dict_haplotypecaller, bam_haplotypecaller)
haplotypecaller = intervals.combine(haplotypecaller_index)

process HaplotypeCaller {
    tag "$intervals"
    container 'broadinstitute/gatk:latest'

    memory threadmem

    input:
    set val(intervals), file(fasta), file(fai), file(dict),
    val(patientId), val(sampleId), val(status), val(name), file(bam), file(bai) from haplotypecaller

    output:
    file("${name}.g.vcf") into haplotypecaller_gvcf
    file("${name}.g.vcf.idx") into index
    val(name) into name_mergevcfs

    when: !params.skip_haplotypecaller

    script:
    """
    gatk HaplotypeCaller \
    --java-options -Xmx${task.memory.toMega()}M \
    -R $fasta \
    -O ${name}.g.vcf \
    -I $bam \
    -ERC GVCF \
    -L $intervals
    """
}



process MergeVCFs {
    tag "${name[0]}.g.vcf"
    publishDir "${params.outdir}/GermlineVariantCalling", mode: 'copy'
    container 'broadinstitute/gatk:latest'

    input:
    file ('*.g.vcf') from haplotypecaller_gvcf.collect()
    file ('*.g.vcf.idx') from index.collect()
    val name from name_mergevcfs.collect()

    output:
    set file("${name[0]}.g.vcf"), file("${name[0]}.g.vcf.idx") into mergevcfs

    script:
    """
    ## make list of input variant files
    for vcf in \$(ls *vcf); do
    echo \$vcf >> input_variant_files.list
    done
    gatk MergeVcfs \
    --INPUT= input_variant_files.list \
    --OUTPUT= ${name[0]}.g.vcf
    """
}

bamsNormal = Channel.create()
bamsTumour = Channel.create()

bam_mutect.choice( bamsTumour, bamsNormal ) { it[2] =~ 1 ? 0 : 1 }

combined_bam = bamsNormal.combine(bamsTumour, by: 0)

ref_mutect = fasta_mutect.merge(fai_mutect, dict_mutect)
variant_calling = combined_bam.combine(ref_mutect)
variant_calling.into{ mutect; manta_no_bed}


process Mutect2 {
    tag "$bam"
    container 'broadinstitute/gatk:latest'
    publishDir "${params.outdir}/Somatic", mode: 'copy'

    input:
    set val(patientId), val(sampleId), val(status), val(name), file(bam), file(bai),
    val(tumourSampleId), val(tumourStatus), val(tumourName), file(tumourBam), file(tumourBai),
    file(fasta), file(fai), file(dict) from mutect

    output:
    set val("${tumourSampleId}_vs_${sampleId}"), file("${tumourSampleId}_vs_${sampleId}.vcf") into vcf_variant_eval

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

if (params.bed) {
    manta = manta_no_bed.merge(bed)
} else {
    manta = manta_no_bed
}

process Manta {
    tag {tumourSampleId + "_vs_" + sampleId}
    container 'maxulysse/sarek:latest'
    publishDir "${params.outdir}/StructuralVariance/${patientId}", mode: 'copy'

    input:
    set val(patientId), val(sampleId), val(status), val(name), file(bam), file(bai),
    val(tumourSampleId), val(tumourStatus), val(tumourName), file(tumourBam), file(tumourBai),
    file(fasta), file(fai), file(dict), file(bed) from manta

    output:
    set val(patientId), val(sampleId), val(tumourSampleId), file("*.vcf.gz"), file("*.vcf.gz.tbi") into manta_results
    
    script:
    beforeScript = params.bed ? "bgzip --threads ${task.cpus} -c *.bed > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
    options = params.bed ? "--exome --callRegions call_targets.bed.gz" : ""
    """
    ${beforeScript}
    configManta.py \
    --normalBam ${bam} \
    --tumorBam ${tumourBam} \
    --reference ${fasta} \
    ${options} \
    --runDir Manta

    python Manta/runWorkflow.py -m local -j ${task.cpus}
    mv Manta/results/variants/candidateSmallIndels.vcf.gz \
        Manta_${tumourSampleId}_vs_${sampleId}.candidateSmallIndels.vcf.gz
    mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
        Manta_${tumourSampleId}_vs_${sampleId}.candidateSmallIndels.vcf.gz.tbi
    mv Manta/results/variants/candidateSV.vcf.gz \
        Manta_${tumourSampleId}_vs_${sampleId}.candidateSV.vcf.gz
    mv Manta/results/variants/candidateSV.vcf.gz.tbi \
        Manta_${tumourSampleId}_vs_${sampleId}.candidateSV.vcf.gz.tbi
    mv Manta/results/variants/diploidSV.vcf.gz \
        Manta_${tumourSampleId}_vs_${sampleId}.diploidSV.vcf.gz
    mv Manta/results/variants/diploidSV.vcf.gz.tbi \
        Manta_${tumourSampleId}_vs_${sampleId}.diploidSV.vcf.gz.tbi
    mv Manta/results/variants/somaticSV.vcf.gz \
        Manta_${tumourSampleId}_vs_${sampleId}.somaticSV.vcf.gz
    mv Manta/results/variants/somaticSV.vcf.gz.tbi \
        Manta_${tumourSampleId}_vs_${sampleId}.somaticSV.vcf.gz.tbi
    """
}

variant_eval = vcf_variant_eval.merge(fasta_variant_eval, fai_variant_eval, dict_variant_eval)

process VariantEval {
    tag "$vcf"
    container 'broadinstitute/gatk:latest'

    input:
    set val(name), file(vcf), file(fasta), file(fai), file(dict) from variant_eval

    output:
    file("${name}.eval.grp") into variantEvalReport

    when: !params.skip_multiqc

    script:
    // TODO: add dbsnp & gold standard
    """
    touch ${name}.eval.grp

    gatk VariantEval \
    -R ${fasta} \
    --eval:${name} $vcf \
    -O ${name}.eval.grp
    """
}

process kraken2 {
    tag "$name"
    publishDir "${params.outdir}/kraken/${patientId}", mode: 'copy'
    container 'flowcraft/kraken2:2.0.7-1'

    input:
    set val(patientId), val(sampleId), val(status), val(name), file(reads) from reads_kraken

    output:
    set file("${sampleId}_kraken.out"), file("${sampleId}_kraken.report") into kraken_results

    script:
    """
    kraken2 \
    --threads ${task.cpus} \
    --paired \
    --db ${kraken_db} \
    --fastq-input $reads \
    --output ${sampleId}_kraken.out \
    --report ${sampleId}_kraken.report
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
    file (bamQCrecalibrated) from bamQCrecalibratedReport.collect().ifEmpty([])
    file (baseRecalibrator) from baseRecalibratorReport.collect().ifEmpty([])
    file (variantEval) from variantEvalReport.collect().ifEmpty([])
    
    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    """
    multiqc . -m fastqc -m picard -m qualimap -m gatk
    """
}
