#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

	inputPairReads = channel
		.fromFilePairs(params.outputDir + 's1.trimmed_P{1,2}.fq')
		.ifEmpty { error "Cannot find any read file pairs in ${params.rawReadsDir}" }
		.view()

	convertFastqToUnmappedBam(inputPairReads)

	alignReadsToReference(convertFastqToUnmappedBam.out)

}


process convertFastqToUnmappedBam {
	echo true
	container params.gatk4Image 

	input:
	tuple val(name), file(reads)

	output:
	path "s1.unmapped.bam"

	script:
	"""
	gatk FastqToSam \
		-F1 ${reads[0]} \
		-F2 ${reads[1]} \
		-SM s_1 \
		-PL ILLUMINA \
		-RG s_1 \
		-O s1.unmapped.bam
	"""

}

process alignReadsToReference {
	echo true
	container params.bwaImage

	input:
	path "s1.unmapped.bam"

	output:
	stdout

	script:
	"""
	bwa mem
	"""

}


process viewUnmappedBam {
	echo true
	container params.samtoolsImage

	input:
	path unmappedBam

	script:
	"""
	samtools view ${unmappedBam} | head -5
	"""
}
