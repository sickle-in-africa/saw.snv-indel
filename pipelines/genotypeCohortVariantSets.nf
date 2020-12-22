#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

	Channel
		.fromPath(params.variantSetsDir + params.cohortId + "*.g.vcf.gz")
		.map { path -> "-V ${path} " } \
		| collect \
		| map { x -> x.join(" ") } \
		| combineSampleGvcfFiles \
		| genotypeCombinedGvcfFile

}

process combineSampleGvcfFiles {
    container params.gatk4Image

	input:
	val gvcfList

	output:
	path "${params.cohortId}.g.vcf"

	script:
	"""
	gatk CombineGVCFs \
		-R ${params.referenceSequence['path']} \
		${gvcfList} \
		-O ${params.cohortId}.g.vcf
	"""
}

process genotypeCombinedGvcfFile {
    container params.gatk4Image

	input:
	path combinedGvcfFile

	script:
	"""
	gatk GenotypeGVCFs \
		-R ${params.referenceSequence['path']} \
		-V ${params.cohortId}.g.vcf \
		-O ${params.variantSetsDir}${params.cohortId}.genotyped.g.vcf
	"""
} 