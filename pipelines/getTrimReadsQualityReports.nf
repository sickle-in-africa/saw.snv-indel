#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

	Channel
		.fromFilePairs(params.outputDir + 'trimmedReads/' + params.readFilePairGlob)
		.ifEmpty { error "Cannot find any read trimmed file pairs matching: ${params.readFilePairGlob}" } \
		| getFastaQualityReport
}

process getFastaQualityReport {
    beforeScript "source ${params.processConfigFile}"
    container params.fastqcImage
    clusterOptions = params.clusterOptions
    queue = params.serverOptions['queue']
    time =  params.serverOptions['time']

	input:
	tuple val(name), file(reads)

    script:
    """
    mkdir -p ${params.outputDir}quality-reports/trimmed/
    fastqc \
    	-t ${params.nThreadsPerProcess} \
    	-o ${params.outputDir}quality-reports/trimmed/ \
    	${reads[0]} ${reads[1]}

    """
}
