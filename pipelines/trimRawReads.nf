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
    time =  '15m'
    jobName = 'trimRawReads'

	input:
	tuple val(name), path(reads)

	script:
	"""
	mkdir -p ${params.outputDir}trimmedReads/
	java -jar ${params.trimmomaticJar} PE \
		-threads ${params.threads} \
		${reads[0]} ${reads[1]} \
		${params.outputDir}trimmedReads/${name}_R1.fq ${params.outputDir}trimmedReads/${name}_UR1.fq \
		${params.outputDir}trimmedReads/${name}_R2.fq ${params.outputDir}trimmedReads/${name}_UR2.fq \
		LEADING:${params.trimLeadX} \
		TRAILING:${params.trimTrailX} \
		SLIDINGWINDOW:5:${params.trimMinAverageQuality} \
		MINLEN:${params.trimMinReadLength}
	"""
}


