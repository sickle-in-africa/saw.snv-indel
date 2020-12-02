#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

	simulationInputs = channel
		.fromPath(params.simulationInputs)
		.ifEmpty { error "Cannot find any input file in ${params.projectDir}" }
		.view()

	writeSampleIdList(simulationInputs) | view

	//mutateReference(simulationInputs)
	
	//generateReads(mutateReference.out)

	//sortTruthVcfFile(mutateReference.out)
}

process mutateReference {
	container params.simulatorImage

	input:
	path simulationInputs

	output:
	tuple path(simulationInputs), path("${params.referenceSequence['label']}.mutated.fa"), path("c_1.truth.g.vcf")

	script:
	"""
	simulate MutateReference \
		--input_ref ${params.referenceSequence['path']} \
		--output_ref ${params.referenceSequence['label']}.mutated.fa \
		--input_json ${simulationInputs} \
		--output_vcf c_1.truth.g.vcf
	"""
}

process generateReads {
	container params.simulatorImage

	input:
	tuple path(simulationInputs), path(mutatedReference), path(truthVcfFile)

	script:
	"""
	simulate GenerateReads \
		--input_ref ${mutatedReference} \
		--reads_prefix ${params.rawReadsDir}c1.raw \
		--input_json ${simulationInputs}
	"""
}

process sortTruthVcfFile {
	container params.bcftoolsImage

	input:
	tuple path(simulationInputs), path(mutatedReference), path(truthVcfFile)

	script:
	"""
	bcftools sort \
		${truthVcfFile} \
		-o ${params.outputDir}/c_1.truth.g.vcf
	"""
}

process writeSampleIdList {

	input:
	path simulationInputs

	output:
	stdout

	script:
	"""
	python3 ${params.toolsDir}writeSimulationSampleList.py \
		--inputJSON ${simulationInputs} \
		--readsDir ${params.rawReadsDir}
	"""
}