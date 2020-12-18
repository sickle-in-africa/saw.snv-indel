#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

	Channel
		.fromFilePairs(params.outputDir + 'trimmedReads/' + params.readFilePairGlob)
		.ifEmpty { error "Cannot find any read trimmed file pairs matching: ${params.readFilePairGlob}" } \
		| alignReadsToReference \
		| addReadGroupInfo \
		| markDuplicateReads \
        | addTagInfo \
		| checkBamFile \
		| saveBamFileToOutputDir
}

process alignReadsToReference {
    beforeScript "source ${params.processConfigFile}"
    container params.bwaImage
    clusterOptions = params.clusterOptions
    queue = params.serverOptions['queue']
    time = '09:00:00' // params.serverOptions['time']

	input:
	tuple val(name), path(readsFilePair)

	output:
	tuple val(name), path("${name}.bam")

	script:
	"""
	bwa mem \
		-t ${params.nThreadsPerProcess} \
		${params.referenceSequence['dir']}bwa.${params.referenceSequence['label']} \
		${readsFilePair[0]} ${readsFilePair[1]} | \
			samtools view \
				-b - | \
				samtools sort \
					-@ ${params.nThreadsPerProcess} \
					-o ${name}.bam
	"""

}

process addReadGroupInfo {
    beforeScript "source ${params.processConfigFile}"
    container params.gatk4Image
    clusterOptions = params.clusterOptions
    queue = params.serverOptions['queue']
    time = '01:00:00' // params.serverOptions['time']

	input:
	tuple val(name), path(bamFile)

	output:
	tuple val(name), path("${name}.withRGs.bam")

	script:
	"""
	gatk AddOrReplaceReadGroups \
		-I ${bamFile} \
		-O ${name}.withRGs.bam \
		-SM ${name} \
		-PL ILLUMINA \
		-PU abcd \
		-LB 1234
	"""
}

process markDuplicateReads {
    beforeScript "source ${params.processConfigFile}"
    container params.gatk4Image
    clusterOptions = params.clusterOptions
    queue = params.serverOptions['queue']
    time = '01:30:00' //params.serverOptions['time']

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

process addTagInfo {
    beforeScript "source ${params.processConfigFile}"
    container params.gatk4Image
    clusterOptions = params.clusterOptions
    queue = params.serverOptions['queue']
    time = '01:30:00' //params.serverOptions['time']

    input:
    tuple val(name), path(bamFile)

    output:
    tuple val(name), path("${name}.fixed.bam")

    script:
    """
    gatk SetNmMdAndUqTags \
        -R ${params.referenceSequence['path']} \
        -I ${bamFile} \
        -O "${name}.fixed.bam"
    """

}

process checkBamFile {
    beforeScript "source ${params.processConfigFile}"
    container params.gatk4Image
    clusterOptions = params.clusterOptions
    queue = params.serverOptions['queue']
    time = '01:30:00' //params.serverOptions['time']

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
		-M VERBOSE
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
    time = '02:00:00' //params.serverOptions['time']

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

