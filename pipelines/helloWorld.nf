#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

	Channel
		.of(1..5) \
		| printHelloWorldMessage \
		| view
}

process printHelloWorldMessage {

	input:
	val x

	output:
	stdout

	script:
	"""
	echo "${x}: hello world"
	echo "${launchDir}"
	echo "${params.rawReadsDir}"
	"""
}