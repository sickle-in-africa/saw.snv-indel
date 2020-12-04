#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

	Channel
		.fromFilePairs(params.outputDir + 'trimmedReads/' + params.readFilePairGlob)
		.ifEmpty { error "Cannot find any read trimmed file pairs matching: ${params.readFilePairGlob}" } \
		| getFastaQualityReport
	
}

process getFastaQualityReport {
	container params.fastqcImage

	input:
	tuple val(name), file(reads)

    script:
    """
    mkdir -p ${params.outputDir}fastqc/trimmed/
    fastqc \
    	-t ${params.threads} \
    	-o ${params.outputDir}fastqc/trimmed/ \
    	${reads[0]} ${reads[1]}

    """
}