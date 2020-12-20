#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

	Channel
		.fromPath(params.alignedReadsDir + params.cohortId + "*.bam")
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
    container params.samtoolsImage

    label 'parallel'

	input:
	tuple val(bamId), path(bamFile)

	output:
	tuple val(bamId), path(bamFile), path("${bamFile}.bai")

	script:
	"""
	samtools index \
		-b \
		-@ ${params.nThreadsPerProcess} \
		${bamFile}
	"""
}

process callVariantsForEachSample {
    container params.gatk4Image
    label 'bigMemory'
    label 'bigDuration'
    label 'parallel'

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
    container params.gatk4Image

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
    container params.gatk4Image

	input:
	path combinedGvcfFile

	script:
	"""
	gatk GenotypeGVCFs \
		-R ${params.referenceSequence['path']} \
		-V ${params.cohortId}.g.vcf \
		-O ${params.variantSetsDir}${params.cohortId}.genotyped.g.vcf
	"""
}
