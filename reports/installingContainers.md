installing containers
=====================

we use singularity. The [Biocontainers project]() and [galaxy project depot](https://depot.galaxyproject.org/singularity/) are particularly helpful.

samtools
--------
```
$ cd /path/to/containers
$ singularity pull https://depot.galaxyproject.org/singularity/samtools:1.9--h91753b0_8
```
sambamba
--------
```
$ singularity pull https://depot.galaxyproject.org/singularity/sambamba:0.7.1--h984e79f_3 
```
I haven't had succes with sambamba so far. 

trimmomatic
-----------
```
$ cd /path/to/containers
$ singularity pull https://depot.galaxyproject.org/singularity/trimmomatic:0.39--1
```
This container actually contains a jar file, `trimmomatic.jar`, which is *inside the container*, at the following path: `/usr/local/share/trimmomatic-0.39-1/trimmomatic.jar`. 

To check this, run:
```
$ cd /path/to/containers
$ singularity run trimmomatic:039--1
Singularity> ls /usr/local/share
trimmomatic-0.39-1/
```
Then, if you search the contents of `/usr/local/share` on your own machine outside of the container, you should not find this folder. 

fastqc
------
```
cd /path/to/containers
$ singularity pull https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0
```

bwa
---
```
$ cd /path/to/containers
$ singularity pull https://depot.galaxyproject.org/singularity/bwa:0.7.8--hed695b0_5
```
We need to make our own container here that contains both `bwa` and `samtools`. They can be used together in the same process (tha read alignment one) but this means we need a single container. Need to write a Dockerfile and make the image. 