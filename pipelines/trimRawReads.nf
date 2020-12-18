#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

	Channel
		.fromFilePairs(params.rawReadsDir + params.readFilePairGlob)
		.ifEmpty { error "Cannot find any read file pairs in ${params.rawReadsDir}" } \
		| trimRawReads
}

process trimRawReads {
    container params.trimmomaticImage

	input:
	tuple val(name), path(reads)

	script:
	"""
	mkdir -p ${params.trimmedReadsDir}
	java -jar ${params.trimmomaticJar} PE \
		-threads ${params.nThreadsPerProcess} \
		${reads[0]} ${reads[1]} \
		${params.trimmedReadsDir}${reads[0]} ${params.trimmedReadsDir}unpaired_${reads[0]} \
		${params.trimmedReadsDir}${reads[1]} ${params.trimmedReadsDir}unpaired_${reads[1]} \
		LEADING:${params.trimLeadX} \
		TRAILING:${params.trimTrailX} \
		SLIDINGWINDOW:5:${params.trimMinAverageQuality} \
		MINLEN:${params.trimMinReadLength}
	"""
}


