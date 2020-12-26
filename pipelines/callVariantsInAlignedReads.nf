#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

	Channel
		.fromPath(params.alignedReadsDir + params.cohortId + "*.bam")
		.map { file -> tuple(file.baseName, file) } \
		| indexInputBamFile \
		| callVariantsForEachSample
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

    label 'reallyBigDuration'
    label 'parallel'

	input:
	tuple val(bamId), path(bamFile), path(bamIndex)

	script:
	"""
    mkdir -p ${params.variantSetsDir}
	gatk HaplotypeCaller \
		-R ${params.referenceSequence['path']} \
		-I ${bamFile} \
		-O ${params.variantSetsDir}${bamId}.g.vcf \
		--lenient true \
		-ERC GVCF
	"""
}
