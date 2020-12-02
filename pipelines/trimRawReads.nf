#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

	inputPairReads = channel
		.fromFilePairs(params.rawReadsDir + 'c1.raw_{1,2}P.fq')
		.ifEmpty { error "Cannot find any read file pairs in ${params.rawReadsDir}" }
		.view()

	trimRawReads(inputPairReads)

}

process trimRawReads {
	echo true
	container params.trimmomaticImage

	input:
	tuple val(name), file(reads)

	script:
	"""
	java -jar ${params.trimmomaticJar} PE \
		-threads ${params.threads} \
		${reads[0]} ${reads[1]} \
		${params.outputDir}${name}.trimmed_P1.fq ${params.outputDir}${name}.trimmed_UP1.fq \
		${params.outputDir}${name}.trimmed_P2.fq ${params.outputDir}${name}.trimmed_UP2.fq \
		LEADING:${params.trimLeadX} \
		TRAILING:${params.trimTrailX} \
		SLIDINGWINDOW:5:${params.trimMinAverageQuality} \
		MINLEN:${params.trimMinReadLength}
	"""
}


