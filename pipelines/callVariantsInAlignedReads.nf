#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

	Channel
		.fromPath(params.outputDir + 'aligned/' + params.cohortId + "*.bam")
		.map { file -> tuple(file.baseName, file) } \
		| indexInputBamFile \
		| callVariantsForEachSample \
		| map { path -> "-V ${path} " } \
		| collect \
		| map { x -> x.join(" ") } \
		| combineSampleGvcfFiles \
		| genotypeCombinedGvcfFile

}

process indexInputBamFile {
    beforeScript "source ${params.processConfigFile}"
    container params.samtoolsImage
    clusterOptions = params.clusterOptions
    queue = params.serverOptions['queue']
    time =  '10m'
    jobName = 'indexInputBamFile'

	input:
	tuple val(bamId), path(bamFile)

	output:
	tuple val(bamId), path(bamFile), path("${bamFile}.bai")

	script:
	"""
	samtools index \
		-b \
		-@ ${params.threads} \
		${bamFile}
	"""
}

process callVariantsForEachSample {
    beforeScript "source ${params.processConfigFile}"
    container params.gatk4Image
    clusterOptions = params.clusterOptions
    queue = params.serverOptions['queue']
    time =  '10m'
    jobName = 'callVariantsForEachSample'

	input:
	tuple val(bamId), path(bamFile), path(bamIndex)

	output:
	path "${bamId}.g.vcf"

	script:
	"""
	gatk HaplotypeCaller \
		-R ${params.referenceSequence['path']} \
		-I ${bamFile} \
		-O ${bamId}.g.vcf \
		--lenient true \
		-ERC GVCF
	"""
}


process combineSampleGvcfFiles {
    beforeScript "source ${params.processConfigFile}"
    container params.gatk4Image
    clusterOptions = params.clusterOptions
    queue = params.serverOptions['queue']
    time =  '10m'
    jobName = 'combineSampleGvcfFiles'

	input:
	val gvcfList

	output:
	path "${params.cohortId}.g.vcf"

	script:
	"""
	gatk CombineGVCFs \
		-R ${params.referenceSequence['path']} \
		${gvcfList} \
		-O ${params.cohortId}.g.vcf
	"""
}

process genotypeCombinedGvcfFile {
    beforeScript "source ${params.processConfigFile}"
    container params.gatk4Image
    clusterOptions = params.clusterOptions
    queue = params.serverOptions['queue']
    time =  '10m'
    jobName = 'genotypeCombinedGvcfFile'

	input:
	path combinedGvcfFile

	script:
	"""
	gatk GenotypeGVCFs \
		-R ${params.referenceSequence['path']} \
		-V ${params.cohortId}.g.vcf \
		-O ${params.outputDir}${params.cohortId}.genotyped.g.vcf
	"""
}
