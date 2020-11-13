#
#  LOCATIONS
#
################

pro_dir=/home/jackmo/computer/genemap/saw.snv-indel

tls_dir=${pro_dir}/tools
dat_dir=${pro_dir}/data
rds_dir=${dat_dir}/reads
ref_dir=${dat_dir}/references
tmp_dir=${dat_dir}/temp
log_dir=${dat_dir}/logs
med_dir=${dat_dir}/media

pip_dir=${pro_dir}/pipelines

## alignment & variant calling tools
bowtie2=${tls_dir}/bowtie2-2.4.2-linux-x86_64
samtools=${tls_dir}/samtools-1.11/samtools
bwa=${tls_dir}/bwa-0.7.17/bwa
gatk=${tls_dir}/gatk-4.1.7.0/gatk
fastqc=${tls_dir}/FastQC/fastqc
trimmomatic=${tls_dir}/Trimmomatic-0.39/trimmomatic-0.39.jar