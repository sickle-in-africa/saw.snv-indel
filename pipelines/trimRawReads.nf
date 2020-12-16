#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

	Channel
		.fromFilePairs(params.rawReadsDir + params.readFilePairGlob)
		.ifEmpty { error "Cannot find any read file pairs in ${params.rawReadsDir}" } \
		| trimRawReads
}

process trimRawReads {
    beforeScript "source ${params.processConfigFile}"
    container params.trimmomaticImage
    clusterOptions = params.clusterOptions
    queue = params.serverOptions['queue']
    time =  params.serverOptions['time']

	input:
	tuple val(name), path(reads)

	script:
	"""
	mkdir -p ${params.outputDir}trimmedReads/
	java -jar ${params.trimmomaticJar} PE \
		-threads ${params.nThreadsPerProcess} \
		${reads[0]} ${reads[1]} \
		${params.outputDir}trimmedReads/${reads[0]} ${params.outputDir}trimmedReads/unpaired_${reads[0]} \
		${params.outputDir}trimmedReads/${reads[1]} ${params.outputDir}trimmedReads/unpaired_${reads[1]} \
		LEADING:${params.trimLeadX} \
		TRAILING:${params.trimTrailX} \
		SLIDINGWINDOW:5:${params.trimMinAverageQuality} \
		MINLEN:${params.trimMinReadLength}
	"""
}


