#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

	Channel
		.fromFilePairs("${params.trimmedReadsDir}/${params.readFilePairGlob}")
		.ifEmpty { error "Cannot find any read trimmed file pairs matching: ${params.trimmedReadsDir}/${params.readFilePairGlob}" } \
		| getFastaQualityReport
}

process getFastaQualityReport {
    label 'withMaxMemory'
    label 'withMaxCpus'
    label 'withMaxTime'
    container params.fastqcImage

	input:
	tuple val(name), file(reads)

    script:
    """
    mkdir -p ${params.qualityReportsDir}/trimmed
    fastqc \
    	-t ${task.cpus} \
    	-o ${params.qualityReportsDir}/trimmed \
    	${reads[0]} ${reads[1]}
    """
}
