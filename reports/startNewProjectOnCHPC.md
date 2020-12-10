Start a new project on CHPC
===========================

1. create new project directory
```
$ touch /my/project/directory 
```
2. create sort link to pipelines directory
```
$ cd /my/project/directory
$ ln -s /path/to/sequence-analysis-workflows/snv-indel/pipelines/ ./pipelines
```
3. copy template parameter file to project directory
```
$ cd /my/project/directory
$ cp /path/to/sequence-analysis-workflows/snv-indel/example.config ./nextflow.config
```
