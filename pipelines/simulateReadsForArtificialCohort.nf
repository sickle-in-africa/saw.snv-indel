#!/usr/bin/env nextflow
nextflow.enable.dsl=2


workflow {

	Channel
		.fromPath(params.simulationInputs) \
		| mutateReference \
		| sortTruthVcfFile

	Channel
		.of(1..params.nSamples)
		.map { x -> "${params.cohortId}_S${x}" }
		.map { x -> "${params.rawReadsDir}${x}" }
		.combine(sortTruthVcfFile.out) \
		| generateReads
}


process mutateReference {
    beforeScript "source ${params.processConfigFile}"
    container params.simulatorImage
    clusterOptions = params.clusterOptions
    queue = params.serverOptions['queue']
    time =  '15m'
    jobName = 'mutateReference'

	container params.simulatorImage

	input:
	path simulationInputs

	output:
	tuple path(simulationInputs), path("${params.referenceSequence['label']}.mutated.fa"), path("${params.cohortId}.truth.g.vcf")

	script:
	"""
	simulate MutateReference \
		--input_ref ${params.referenceSequence['path']} \
		--output_ref ${params.referenceSequence['label']}.mutated.fa \
		--input_json ${simulationInputs} \
		--output_vcf ${params.cohortId}.truth.g.vcf
	"""
}


process sortTruthVcfFile {
    beforeScript "source ${params.processConfigFile}"
    container params.bcftoolsImage
    clusterOptions = params.clusterOptions
    queue = params.serverOptions['queue']
    time =  '5m'
    jobName = 'sortTruthVcfFile'

	input:
	tuple path(simulationInputs), path(mutatedReference), path(truthVcfFile)

	output:
	tuple path(simulationInputs), path(mutatedReference)

	script:
	"""
	bcftools sort \
		${truthVcfFile} \
		-o ${params.outputDir}/${params.cohortId}.truth.g.vcf
	"""
}


process generateReads {
    beforeScript "source ${params.processConfigFile}"
    container params.simulatorImage
    clusterOptions = params.clusterOptions
    queue = params.serverOptions['queue']
    time =  '20m'
    jobName = 'generateReads'

	input:
	tuple val(readsPrefix),  path(simulationInputs), path(mutatedReference)

	script:
	"""
	simulate GenerateReads \
		--input_ref ${mutatedReference} \
		--reads_prefix ${readsPrefix} \
		--input_json ${simulationInputs}
	"""
}
