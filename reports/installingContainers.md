Installing containers
=====================

We use singularity to manage software used by the workflow. Sickle In Africa has its own Docker Hub organisation page, which you can view [here](https://hub.docker.com/orgs/sickleinafrica/repositories).

We have followed the do-one-thing rule from functional programming here, and chosen to create a container image for each tool.

samtools
--------
```
$ cd /path/to/containers
$ singularity pull docker://sickleinafrica/samtools:1.11
```

trimmomatic
-----------
```
$ cd /path/to/containers
$ singularity pull docker://sickleinafrica/trimmomatic:0.39
This container actually contains a jar file, `trimmomatic-0.39.jar`, which is *inside the container*, at the following path: `/usr/local/share/applications/Trimmomatic-0.39/trimmomatic-0.39.jar`.  

fastqc
------
```
cd /path/to/containers
$ singularity pull docker://sickleinafrica/fastqc:0.11.9
```

bwa
---
```
$ cd /path/to/containers
$ singularity pull docker://sickleinafrica/bwa:0.7.17
```
We have made our own container here that contains both `bwa 0.7.17` and `samtools 1.11`. They are used together in the same process (alignReadsToReference.nf::alignReadsToReference) which means we need a single container.

sequence-simulator
------------------
```
$ cd /path/to/containers
$ singularity pull docker://sickleinafrica/sequence-simulator:0.1
```

gatk 4
------
```
$ cd /path/to/containers
$ singularity pull docker://broadinstitute/gatk
```

bcftools
--------
```
$ cd /path/to/containers
$ singularity pull docker://sickleinafrica/bcftools:1.11
```

bamstats
--------

Taken from [here](http://bamstats.sourceforge.net/).
```
$ cd /path/to/containers
$ singularity pull docker://sickleinafrica/bamstats:1.25
```

