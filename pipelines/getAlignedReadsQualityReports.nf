#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

    Channel
        .fromPath(params.outputDir + 'aligned/' + params.cohortId + "*.bam")
        .map { file -> tuple(file.baseName, file) }
        .ifEmpty { error 'Cannot find any bam files in output/aligned/' } \
        | getBamStatsReport

/*
    Channel
        .fromPath(params.outputDir + 'aligned/' + params.cohortId + "*.bam")
        .map { file -> tuple(file.baseName, file) }
        .ifEmpty { error 'Cannot find any bam files in output/aligned/' } \
        | getDepthOfCoverageReport
*/
}

process getBamStatsReport {
    beforeScript "source ${params.processConfigFile}"
    container params.bamstatsImage
    clusterOptions = params.clusterOptions
    queue = params.serverOptions['queue']
    time =  params.serverOptions['time']

    input:
    tuple val(bamBase), path(bamFile)

    script:
    """
    mkdir -p ${params.outputDir}quality-reports/aligned/
    java -Xmx8g -jar ${params.bamstatsJar} \
        -i ${bamFile} \
        -v html \
        -o ${params.outputDir}quality-reports/aligned/${bamBase}.bamstats.html 
    """
}

/*
 * THE GATK DEPTHOFCOVERAGE FUNCTION
 * IS NOT YET READY FOR PRODUCTION
 * SO LEAVE THIS PROCESS COMMENTED OUT
 * FOR NOW.
 */
process getDepthOfCoverageReport {
    beforeScript "source ${params.processConfigFile}"
    container params.gatk4Image
    clusterOptions = params.clusterOptions
    queue = params.serverOptions['queue']
    time =  params.serverOptions['time']

    input:
    tuple val(bamBase), path(bamFile)

    script:
    """
    mkdir -p ${params.outputDir}quality-reports/aligned/
    gatk \
        DepthOfCoverage \
        -R ${params.referenceSequence['path']} \
        -O ${params.outputDir}quality-reports/alined/${bamBase}.depthofcov \
        -I ${bamFile} \
        -L []
    """
}
