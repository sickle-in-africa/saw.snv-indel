#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

	buildBwaIndex()

	buildSamtoolsIndex()

	buildGATKDictionary()

}

process buildBwaIndex {
    beforeScript "source ${params.processConfigFile}"
    container params.bwaImage
    clusterOptions = params.clusterOptions
    queue = params.serverOptions['queue']
    time =  params.serverOptions['time']

	when:
	!file(params.referenceSequence['dir'] + 'bwa.' + params.referenceSequence['label'] + '.amb').exists()

	script:
	"""
	bwa index \
			-p ${params.referenceSequence['dir']}bwa.${params.referenceSequence['label']} \
			${params.referenceSequence['path']}
	"""
}

process buildSamtoolsIndex {
    beforeScript "source ${params.processConfigFile}"
    container params.samtoolsImage
    clusterOptions = params.clusterOptions
    queue = params.serverOptions['queue']
    time =  params.serverOptions['time']

	when:
	!file(params.referenceSequence['path'] + '.fai').exists()

	script:
	"""
	samtools faidx ${params.referenceSequence['path']}
	"""
}

process buildGATKDictionary {
    beforeScript "source ${params.processConfigFile}"
    container params.gatk4Image
    clusterOptions = params.clusterOptions
    queue = params.serverOptions['queue']
    time =  params.serverOptions['time']

	when:
	!file(params.referenceSequence['dir'] + params.referenceSequence['label'] + '.dict').exists()

	script:
	"""
	gatk CreateSequenceDictionary \
			-R ${params.referenceSequence['path']} \
			-O ${params.referenceSequence['dir']}${params.referenceSequence['label']}'.dict'
	"""
}
