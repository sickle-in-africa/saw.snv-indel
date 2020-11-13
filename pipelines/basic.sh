#
#  BASIC
#
#	simple variant caller
#
#	Jack Morrice
#
###############################
#!/usr/bin/env bash

source ./base_mod.sh

pipeline() {
	local argv=("@")

	custom_call "testing..." test || exit 1

	custom_call "aligning reads to reference..." align_reads || exit 1

}

#
#  tasks
#
test() {
	echo 'hi'
}

align_reads() {
	${bwa} mem \
		-M \
		-t ${threads} \
		-R "@RG\tID:$sample_id\tSM:$sample_id\tPL:Illumina" \
		${!idx} ${fastq1} ${fastq2} | \
			${samtools} view -bS - | \
				${samtools} sort \
					-@ ${threads} \
					-m ${maxmem} \
					-o ${work_dir}/${prefix}/${prefix}.bam \
					-T ${tmp_dir}/${prefix}
}

#
#  run pipeline
#
pipeline "$@"