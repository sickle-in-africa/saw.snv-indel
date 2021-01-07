#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

	Channel
		.fromPath("${params.alignedReadsDir}/${params.cohortId}*.bam")
		.ifEmpty { error "Cannot find any read file pairs matching: ${params.alignedReadsDir}/${params.cohortId}*.bam" }
		.map { file -> tuple(file.baseName, file) } \
		| indexInputBamFile \
		| callVariantsForEachSample
}


process indexInputBamFile {
	label 'withMaxMemory'
	label 'withMaxCpus'
	label 'withMaxTime'
	container params.samtoolsImage

	input:
	tuple val(bamId), path(bamFile)

	output:
	tuple val(bamId), path(bamFile), path("${bamFile}.bai")

	script:
	"""
	samtools index \
		-b \
		-@ ${task.cpus} \
		${bamFile}
	"""
}

process callVariantsForEachSample {
    label 'withMaxMemory'
    label 'withMaxCpus'
    label 'withMaxTime'
    container params.gatk4Image

	input:
	tuple val(bamId), path(bamFile), path(bamIndex)

	script:
	"""
    mkdir -p ${params.variantSetsDir}
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
		HaplotypeCaller \
		-R ${params.referenceSequence['path']} \
		-I ${bamFile} \
		-O ${params.variantSetsDir}/${bamId}.g.vcf \
		--lenient true \
		-ERC GVCF
	"""
}
