#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

	Channel
		.fromFilePairs(params.outputDir + 'trimmedReads/' + params.readFilePairGlob)
		.ifEmpty { error "Cannot find any read trimmed file pairs matching: ${params.readFilePairGlob}" } \
		| alignReadsToReference \
		| addReadGroupInfo \
		| markDuplicateReads \
		| checkBamFile \
		| saveBamFileToOutputDir

}

process alignReadsToReference {
    beforeScript "source ${params.processConfigFile}"
    container params.bwaImage
    clusterOptions = params.clusterOptions
    queue = params.serverOptions['queue']
    time =  '10m'
    jobName = 'alignReadsToReference'

	input:
	tuple val(name), path(readsFilePair)

	output:
	tuple val(name), path("${name}.bam")

	script:
	"""
	bwa mem \
		-t ${params.threads} \
		${params.referenceSequence['dir']}bwa.${params.referenceSequence['label']} \
		${readsFilePair[0]} ${readsFilePair[1]} | \
			samtools view \
				-b - | \
				samtools sort \
					-@ ${params.threads} \
					-o ${name}.bam
	"""

}

process addReadGroupInfo {
    beforeScript "source ${params.processConfigFile}"
    container params.gatk4Image
    clusterOptions = params.clusterOptions
    queue = params.serverOptions['queue']
    time =  '10m'
    jobName = 'addReadGroupInfo'

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
    beforeScript "source ${params.processConfigFile}"
    container params.gatk4Image
    clusterOptions = params.clusterOptions
    queue = params.serverOptions['queue']
    time =  '10m'
    jobName = 'markDuplicateReads'

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
    beforeScript "source ${params.processConfigFile}"
    container params.gatk4Image
    clusterOptions = params.clusterOptions
    queue = params.serverOptions['queue']
    time =  '10m'
    jobName = 'checkBamFile'

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
    beforeScript "source ${params.processConfigFile}"
    container params.samtoolsImage
    clusterOptions = params.clusterOptions
    queue = params.serverOptions['queue']
    time =  '10m'
    jobName = 'saveBamFileToOutputDir'

	input:
	tuple val(name), path(bamFile)

	script:
	"""
	mkdir -p ${params.outputDir}aligned
	samtools view \
		-b \
		-o ${params.outputDir}aligned/${name}.bam \
		${bamFile}
	"""
}

