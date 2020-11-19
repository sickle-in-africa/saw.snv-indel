#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

	buildBwaIndex()

	buildSamtoolsIndex()

	buildGATKDictionary()

}

process buildBwaIndex {
	container params.bwaImage

	when:
	!file(params.referenceSequence['dir'] + 'bwa.' + params.referenceSequence['label'] + '.amb').exists()

	script:
	"""
	bwa index \
			-p ${params.referenceSequence['dir']}bwa.${params.referenceSequence['label']} \
			-a is \
			${params.referenceSequence['path']}
	"""
}

process buildSamtoolsIndex {
	container params.samtoolsImage

	when:
	!file(params.referenceSequence['path'] + '.fai').exists()

	script:
	"""
	samtools faidx ${params.referenceSequence['path']}
	"""
}

process buildGATKDictionary {
	container params.gatk4Image

	when:
	!file(params.referenceSequence['dir'] + params.referenceSequence['label'] + '.dict').exists()

	script:
	"""
	gatk CreateSequenceDictionary \
			-R ${params.referenceSequence['path']} \
			-O ${params.referenceSequence['dir']}${params.referenceSequence['label']}'.dict'
	"""
}