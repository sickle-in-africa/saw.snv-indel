#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

    Channel
        .fromFilePairs("${params.trimmedReadsDir}/${params.readFilePairGlob}")
        .ifEmpty { error "Cannot find any read trimmed file pairs matching: ${params.trimmedReadsDir}/${params.readFilePairGlob}" } \
        | alignReadsToReference \
        | markDuplicateReads \
        | addTagInfo \
        | recalibrateBaseQualityScores \
        | checkBamFile \
        | saveBamFileToOutputDir
}

process alignReadsToReference {
    label 'withMaxMemory'
    label 'withMaxCpus'
    label 'withMaxTime'
    container params.bwaImage

    input:
    tuple val(name), path(readsFilePair)

    output:
    tuple val(name), path("${name}.bam")

    script:
    """
    bwa mem \
        -t ${task.cpus} \
        -R \"@RG\\tID:${params.runId}\\tPU:${params.runId}\\tSM:${name}\\tLB:${name}\\tPL:${params.sequencingPlatform}\" \
        ${params.referenceSequence['dir']}/bwa.${params.referenceSequence['label']} \
        ${readsFilePair[0]} ${readsFilePair[1]} | \
            samtools view \
                -b - | \
                samtools sort \
                    -@ ${task.cpus} \
                    -o ${name}.bam
    """
}

process markDuplicateReads {
    label 'withMaxMemory'
    label 'withMaxCpus'
    label 'withMaxTime'
    container params.gatk4Image

    input:
    tuple val(name), path(bamFile)

    output:
    tuple val(name), path("${name}.marked.bam")

    script:
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        MarkDuplicates \
        -I ${bamFile} \
        -O "${name}.marked.bam" \
        -M "${name}.marked_dup_metrics.txt"
    """
}

process addTagInfo {
    label 'withMaxMemory'
    label 'withMaxCpus'
    label 'withMaxTime'
    container params.gatk4Image

    input:
    tuple val(name), path(bamFile)

    output:
    tuple val(name), path("${name}.fixed.bam")

    script:
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        SetNmMdAndUqTags \
        -R ${params.referenceSequence['path']} \
        -I ${bamFile} \
        -O "${name}.fixed.bam"
    """

}

process recalibrateBaseQualityScores {
    label 'withMaxMemory'
    label 'withMaxCpus'
    label 'withMaxTime'
    container params.gatk4Image

    input:
    tuple val(name), path(bamFile)

    output:
    tuple val(name), path("${name}.recal.bam")

    script:
    if ( file("${params.baseQualityRecalibrationTable}").exists() )
        """
        gatk ApplyBQSR \
            -R reference.fasta \
            -I ${bamFile} \
            --bqsr-recal-file ${params.baseQualityRecalibrationTable} \
            -O ${name}.recal.bam
        """
    else
        """
        ln -s ${bamFile} ${name}.recal.bam
        """
}

process checkBamFile {
    label 'withMaxMemory'
    label 'withMaxCpus'
    label 'withMaxTime'
    container params.gatk4Image

    input:
    tuple val(name), path(bamFile)

    output:
    tuple val(name), path("${bamFile}")

    script:
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        ValidateSamFile \
        -I ${bamFile} \
        -R ${params.referenceSequence['path']} \
        --TMP_DIR ${params.tempDir} \
        -M VERBOSE
    """
}

process saveBamFileToOutputDir {
    label 'withMaxMemory'
    label 'withMaxCpus'
    label 'withMaxTime'
    container params.samtoolsImage

    input:
    tuple val(name), path(bamFile)

    script:
    """
    mkdir -p ${params.alignedReadsDir}
    samtools view \
        -b \
        -o ${params.alignedReadsDir}/${name}.bam \
        ${bamFile}
    """
}

