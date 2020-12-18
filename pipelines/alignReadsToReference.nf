#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

	Channel
		.fromFilePairs(params.trimmedReadsDir + params.readFilePairGlob)
		.ifEmpty { error "Cannot find any read trimmed file pairs matching: ${params.readFilePairGlob}" } \
		| alignReadsToReference \
		| addReadGroupInfo \
		| markDuplicateReads \
        | addTagInfo \
        | recalibrateBaseQualityScores \
		| checkBamFile \
		| saveBamFileToOutputDir
}

process alignReadsToReference {
    container params.bwaImage
    label 'bigMemory'
    label 'bigDuration'
    label 'parallel'

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
		-SM ${name} \
		-PL ILLUMINA \
		-PU abcd \
		-LB 1234
	"""
}

process markDuplicateReads {
    container params.gatk4Image
    label 'bigMemory'

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
    container params.gatk4Image

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

process recalibrateBaseQualityScores {
	container params.gatk4Image

	input:
	tuple val(name), path(bamFile)

    output:
    tuple val(name), path("${name}.recal.bam")

	script:
	if ( file("${params.baseQualityRecalibrationTable}").exists() )
		"""
		gatk ApplyBQSR \
			-R reference.fasta \
			-I ${bamFile} \
			--bqsr-recal-file ${params.baseQualityRecalibrationTable} \
			-O ${name}.recal.bam
	   	"""
	else
		"""
		ln -s ${bamFile} ${name}.recal.bam
		"""
}

process checkBamFile {
    container params.gatk4Image

	input:
	tuple val(name), path(bamFile)

	output:
	tuple val(name), path("${bamFile}")

	script:
	"""
	gatk ValidateSamFile \
		-I ${bamFile} \
		-R ${params.referenceSequence['path']} \
		--TMP_DIR ${params.tempDir} \
		-M VERBOSE
	"""
}

process saveBamFileToOutputDir {
    container params.samtoolsImage

	input:
	tuple val(name), path(bamFile)

	script:
	"""
	mkdir -p ${params.alignedReadsDir}
	samtools view \
		-b \
		-o ${params.alignedReadsDir}${name}.bam \
		${bamFile}
	"""
}

