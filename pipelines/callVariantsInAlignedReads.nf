#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

	inputAlignedReads = channel
		.fromPath(params.outputDir + 's1.aligned.bam')
		.ifEmpty { error "Cannot find any aligned read file in ${params.outputDir}" }
		.view()

	indexInputBamFile(inputAlignedReads) \
		| callVariantsPerSample \
		| saveVCFToOutputDir
}

process indexInputBamFile {
	container params.samtoolsImage

	input:
	path bamFile

	output:
	tuple val('s1'), path("s1.aligned.bam"), path("s1.aligned.bam.bai")

	script:
	"""
	samtools index \
		-b \
		-@ ${params.threads} \
		${bamFile}
	"""
}

process callVariantsPerSample {
	container params.gatk4Image

	input:
	tuple val(sampleID), path(alignedBamFile), path(alignedBamIndex)

	output:
	tuple val("${sampleID}"), path("${sampleID}.raw.g.vcf")

	script:
	"""
	gatk HaplotypeCaller \
		-R ${params.referenceSequence['path']} \
		-I ${alignedBamFile} \
		-O ${sampleID}.raw.g.vcf
	"""
}

process saveVCFToOutputDir {

	input:
	tuple val(sampleID), path(vcfFile)

	script:
	"""
	cat ${vcfFile} > ${params.outputDir}${sampleID}.raw.g.vcf
	"""
}