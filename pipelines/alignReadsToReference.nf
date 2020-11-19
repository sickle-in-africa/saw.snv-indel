#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

	inputPairReads = channel
		.fromFilePairs(params.outputDir + 's1.trimmed_P{1,2}.fq')
		.ifEmpty { error "Cannot find any read file pairs in ${params.rawReadsDir}" }
		.view()

	alignReadsToReference(inputPairReads) \
		| convertAlignedReadsFromSamToBamAndSort \
		| addReadGroupInfo \
		| markDuplicateReads \
		| checkBamFile \
		| saveBamFileToOutputDir

}

process alignReadsToReference {
	container params.bwaImage

	input:
	tuple val(name), path(readsFilePair)

	output:
	tuple val(name), path("${name}.sam")

	script:
	"""
	bwa mem \
		-t ${params.threads} \
		${params.referenceSequence['dir']}bwa.${params.referenceSequence['label']} \
		${readsFilePair[0]} ${readsFilePair[1]} > ${name}.sam
	"""

}

process convertAlignedReadsFromSamToBamAndSort {
	container params.samtoolsImage

	input:
	tuple val(name), path(samFile)

	output:
	tuple val(name), path("${name}.bam")

	script:
	"""
	samtools view -b ${samFile} | samtools sort -@ ${params.threads} -o ${name}.bam
	"""
}

process addReadGroupInfo {
	container params.gatk4Image

	input:
	tuple val(name), path(bamFile)

	output:
	tuple val(name), path("${name}.withRGs.bam")

	script:
	"""
	gatk AddOrReplaceReadGroups \
		-I ${bamFile} \
		-O ${name}.withRGs.bam \
		-SM not_sure \
		-PL ILLUMINA \
		-PU not_sure \
		-ID ${name} \
		-LB not_sure
	"""
}

process markDuplicateReads {
	container params.gatk4Image

	input:
	tuple val(name), path(bamFile)

	output:
	tuple val(name), path("${name}.marked.bam")

	script:
	"""
	gatk MarkDuplicates \
		-I ${bamFile} \
		-O "${name}.marked.bam" \
		-M "${name}.marked_dup_metrics.txt"
	"""

}

process checkBamFile {
	echo true
	container params.gatk4Image

	input:
	tuple val(name), path(markedBamFile)

	output:
	tuple val(name), path("${markedBamFile}")

	script:
	"""
	gatk ValidateSamFile \
		-I ${markedBamFile} \
		-R ${params.referenceSequence['path']} \
		--TMP_DIR ${params.tempDir} \
		-M SUMMARY
	"""
}

process recalibrateBaseQualityScores {
	container params.gatk4Image

	input:
	tuple val(name), path(markedBamFile)

	script:
	"""
	echo "we will write this process when we are working with human samples. 
	"""
}

process saveBamFileToOutputDir {
	echo true
	container params.samtoolsImage

	input:
	tuple val(name), path(bamFile)

	script:
	"""
	samtools view \
		-b \
		-o ${params.outputDir}s1.aligned.bam \
		${bamFile}
	"""
}

