#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

	inputPairReads = channel
		.fromFilePairs(params.outputDir + 's1.trimmed_P{1,2}.fq')
		.ifEmpty { error "Cannot find any read file pairs in ${params.rawReadsDir}" }
		.view()

	getFastaQualityReport(inputPairReads)
	
}

process getFastaQualityReport {
	echo true
	container params.fastqcImage

	input:
	set val(name), file(reads)

    script:
    """
    fastqc \
    	-t ${params.threads} \
    	-o ${params.outputDir} \
    	${reads[0]} ${reads[1]}

    """
}