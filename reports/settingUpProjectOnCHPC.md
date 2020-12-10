Setting up a newsnv-indel project on CHPC
=========================================

To use the snv-indel workflow on the CHPC cluster, set up a working directory with the following steps:

1. navigate to your personal workspace and create a new project directory:
```
$ cd /mnt/lustre/users/<username>
$ mkdir <new-project>
```
2. create the folders `output`, `temp`, `reads`:
```
$ cd /path/to/new-project
$ mkdir output temp reads
```
3. create a soft link to the snv-indel pipelines directory:
```
$ cd /path/to/new-project
$ ln -s /mnt/lustre/groups/<project-shortname>/SADaCC/sequence-analysis-workflows/snv-indel/pipelines ./
```
4. copy across the config and parameter files:
```
$ cd /path/to/new-project
$ cp /mnt/lustre/groups/<project-shortname>/SADaCC/sequence-analysis-workflows/snv-indel/example-chpc.config ./nextflow.config
$ cp /mnt/lustre/groups/<project-shortname>/SADaCC/sequence-analysis-workflows/snv-indel/config-chpc.sh ./config.sh
$ cp /mnt/lustre/groups/<project-shortname>/SADaCC/sequence-analysis-workflows/snv-indel/example-simulationInputs.json ./simulationInputs.json
$ cp /mnt/lustre/groups/<project-shortname>/SADaCC/sequence-analysis-workflows/snv-indel/ processConfig-chpc.sh ./processConfig.sh
```
5. open the `nextflow.config` file you have just created and make sure all the locations point to the right places, and that the parameter values are correct, for example the number of treads, of the desired CHPC queue.
6. run the workflow scripts in the following order:
```
$ cd /path/to/new-project
$ source config.sh
$ nextflow pipelines/simulateReadsForArtificialCohort.nf
$ nextflow pipelines/getRawReadsQualityReports.nf
$ nextflow pipelines/trimRawReads.nf
$ nextflow pipelines/getTrimReads/QualityReports.nf
$ nextflow pipelines/alignReadsToReference.nf
$ nextflow pipeline/callVariantsInAlignedReads.nf
```
the `config.sh` script loads the nextflow module needed to run the workflow scipts. You only need to run this script once per session.
7. verify that the `.vcf` files overlap:
```
$ cd /path/to/new-project/output
$ less <cohort-id>.genotyped.g.vcf
$ less <cohort-id>.truth.g.vcf
```
for example I find the following variant is present in both `.vcf` files:
```
gi|9626243|ref|NC_001416.1|     585     .       T       G       199.05  ...
```
Your will find different variants for a different simulation run.
8. your project directory is now verified and set up for analysis!
