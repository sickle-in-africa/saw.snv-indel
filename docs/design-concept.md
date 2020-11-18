Design concept
==============

We will design the workflow around the following three tools:
1. Nextflow
2. Singularity
3. Conda

Nextflow will be the *workflow description language*. We will specify our pipelines as nextflow scripts. A pipeline is a automated chain of processes. Processes in workflow DSLs have the same place as functions do in other scipting languages. So we will follow the "do-one-thing" rule for processes. 

Singularity will be used to run containers to handle all dependencies. For portbility reasons, *all processes will be run in containters*. Following best practices guides published [here](https://biocontainers-edu.readthedocs.io/en/latest/best_practices.html) by the [Biocontainers](https://biocontainers.pro/#/) project, we will enforce the "one container per process" rule. The focus will be on containers rather than 

Conda will be used in some cases for certain dependencies. However, following the "one container per process" rule, each process will have a container, and if this process requires any packages from conda then *these packages will be installed within a container.*
