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
    time =  params.serverOptions['time']

	input:
	tuple val(name), path(reads)

	script:
	"""
	mkdir -p ${params.outputDir}quality-reports
	fastqc \
		-t ${params.nThreadsPerProcess} \
		-o ${params.outputDir}quality-reports/ \
		${reads[0]} ${reads[1]}
	"""
}
