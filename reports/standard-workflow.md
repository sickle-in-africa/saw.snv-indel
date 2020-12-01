the standard workflow
=====================

the different sub-workflow scripts (pipelines) are saved in
```
saw.snv-indel/pipelines
```
and each script is a step in the standard workflow. Between the execution of each script (sub-workflow) human interaction is usually needed, for example after running `getRawReadsQualityReports.nf` you should check the quality reports before deciding on suitable trimming quality control parameters to use. 

Workflow outline
-----------------

the full workflow is run as:

0. create input parameter file `nextflow.config`
1. run `getRawReadsQualityReports.nf`
2. check output quality reports and chose trimming parameters
3. run `trimRawReads.nf`
4. run `getTrimReadsQualityReports.nf`
5. check the trimmed-reads quality reports; go back to step 2 if necessary
6. run `alignReadsToReference.nf`
7. run `callVariantsInAlignedReads.nf`