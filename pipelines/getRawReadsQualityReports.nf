#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

	Channel
		.fromFilePairs(params.rawReadsDir + params.readFilePairGlob)
		.ifEmpty { error "Cannot find any read file pairs that match: ${params.readFilePairGlob}" } \
		| getFastaQualityReport
}


process getFastaQualityReport {
    beforeScript "source ${params.processConfigFile}"
    container params.fastqcImage
    clusterOptions = params.clusterOptions
    queue = params.serverOptions['queue']
    time =  '15m'
    jobName = 'getFastaQualityReport'

	input:
	tuple val(name), path(reads)

	script:
	"""
	mkdir -p ${params.outputDir}fastqc
	fastqc \
		-t ${params.threads} \
		-o ${params.outputDir}fastqc/ \
		${reads[0]} ${reads[1]}
	"""
}
