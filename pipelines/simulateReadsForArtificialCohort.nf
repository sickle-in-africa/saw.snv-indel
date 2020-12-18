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
	container params.bcftoolsImage

	input:
	tuple path(simulationInputs), path(mutatedReference), path(truthVcfFile)

	output:
	tuple path(simulationInputs), path(mutatedReference)

	script:
	"""
	mkdir -p ${params.variantSetsDir}
	bcftools sort \
		${truthVcfFile} \
		-o ${params.variantSetsDir}/${params.cohortId}.truth.g.vcf
	"""
}


process generateReads {
	container params.simulatorImage

	input:
	tuple val(readsPrefix),  path(simulationInputs), path(mutatedReference)

	script:
	"""
	mkdir -p ${params.rawReadsDir}
	simulate GenerateReads \
		--input_ref ${mutatedReference} \
		--reads_prefix ${readsPrefix} \
		--input_json ${simulationInputs}
	"""
}
