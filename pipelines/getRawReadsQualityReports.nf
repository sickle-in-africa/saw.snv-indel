#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

	Channel
		.fromFilePairs(params.rawReadsDir + params.readFilePairGlob)
		.ifEmpty { error "Cannot find any read file pairs that match: ${params.readFilePairGlob}" } \
		| getFastaQualityReport
}


process getFastaQualityReport {
    container params.fastqcImage

	input:
	tuple val(name), path(reads)

	script:
	"""
	mkdir -p ${params.qualityReportsDir}
	fastqc \
		-t ${params.nThreadsPerProcess} \
		-o ${params.qualityReportsDir} \
		${reads[0]} ${reads[1]}
	"""
}
