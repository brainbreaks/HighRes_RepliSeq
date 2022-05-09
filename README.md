 GRO-seq pipeline package (Wei, Pei-Chi group; DKFZ)
====================================================
Global Run-On Sequencing (GRO-Seq) pipeline for analyzing transcription activity of genes from engaged RNA polymerase. 


Table of Contents
=================

  * Using GRO-seq docker image
    * [Singularity](#singularity) 
      * [Pull Repli-seq image from DockerHUB](#singularity-pull)
      * [Run container](#singularity-run)
      * [Inspect container](#singularity-inspect)
    * [Docker](#docker) 
      * [Pull Repli-seq image from DockerHUB](#docker-pull)
      * [Run container](#docker-run)
      * [Inspect container](#docker-inspect)
    * [Build Repliseq-seq docker image](#build) 
      * [Download required packages and files](#build-download)
      * [Build docker image](#build-build)
      * [Push docker image to Docker HUB](#build-build)
      * [Convert cached docker image to singularity (for local testing)](#build-convert)

<a name="singularity">Singularity</a>
====================================================

<a name="singularity-pull">Pull Repli-seq image from data server</a>
----------------------------------------------------
```console
singularity pull docker://sandrejev/repliseq:latest
```

<a name="singularity-run">Run container</a>
----------------------------------------------------
Before running Repli-seq pipeline you will need to obtain genome(fasta). For mm9, mm10 and hg19 these can be downloaded automatically with `download` command. 
```console
singularity exec -B `pwd` repliseq.sif download mm10
```

Run Repli-seq pipeline
```console
 
singularity exec -B `pwd` repliseq.sif align -m ~/Workspace/Datasets/Repliseq/raw/B400_RS_001_23911/23911_meta.tsv -g ~/Workspace/genomes/mm10/mm10 -o ~/Workspace/Datasets/Repliseq -t 32 --fastq-dir ~/Workspace/Datasets/Repliseq/raw/B400_RS_001_23911
singularity exec -B `pwd` repliseq.sif coverage -m ~/Workspace/Datasets/Repliseq/raw/B400_RS_001_23911/23911_meta.tsv -o ~/Workspace/Datasets/Repliseq --binsizes 30000,40000,60000
```


<a name="singularity-inspect">Inspect container</a>
----------------------------------------------------
```console
singularity shell repliseq.sif
```

<a name="docker">Docker</a>
====================================================

<a name="docker-pull">Pull Repli-seq image from DockerHUB</a>
----------------------------------------------------
```console
docker pull sandrejev:repliseq
```

<a name="docker-run">Run container</a>
----------------------------------------------------
Before running Repli-seq pipeline you will need to obtain genome(fasta). For mm9, mm10 and hg19 these can be downloaded automatically with `download` command. 
```console
docker run -v ${PWD}:/mount -u $(id -g ${USER}):$(id -g ${USER}) -it --entrypoint download repliseq mm10
```

Run groseq pipeline. Keep in mind that annotation file (`-a` flag) is created automatically for you from geneRef.gtf
```console
docker run -v ${PWD}:/mount -u $(id -g ${USER}):$(id -g ${USER}) -it --entrypoint align repliseq -m ~/Workspace/Datasets/Repliseq/raw/B400_RS_001_23911/23911_meta.tsv -g ~/Workspace/genomes/mm10/mm10 -o ~/Workspace/Datasets/Repliseq -t 32 --fastq-dir ~/Workspace/Datasets/Repliseq/raw/B400_RS_001_23911
docker run -v ${PWD}:/mount -u $(id -g ${USER}):$(id -g ${USER}) -it --entrypoint coverage repliseq -m ~/Workspace/Datasets/Repliseq/raw/B400_RS_001_23911/23911_meta.tsv -o ~/Workspace/Datasets/Repliseq --binsizes 30000,40000,60000
```

<a name="docker-inspect">Inspect container</a>
----------------------------------------------------
```console
docker run -v ${PWD}:/mount -u $(id -g ${USER}):$(id -g ${USER}) -it --entrypoint bash repliseq
```

<a name="build">Building Repli-seq docker image</a>
====================================================


<a name="build-build">Build docker image</a>
----------------------------------------------------
To build Docker image you need to execute
```console
docker build --squash --build-arg http_proxy="http://www-int2.inet.dkfz-heidelberg.de:80" --build-arg https_proxy="http://www-int2.inet.dkfz-heidelberg.de:80" --progress=plain --rm -t sandrejev/repliseq:latest .
```

<a name="build-push">Push docker image to Docker HUB</a>
----------------------------------------------------
```console
docker login
docker push sandrejev/repliseq:latest
```

<a name="build-convert">Convert cached docker image to singularity (for local testing)</a>
----------------------------------------------------
```console
singularity pull docker-daemon:sandrejev/repliseq:latest
```