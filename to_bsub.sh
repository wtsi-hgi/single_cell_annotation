#!/usr/bin/env bash

# temporary cache dir used by singularity when pulling images from dockerhub
# this is different than the nextflow singularity cache dir when it caches images (defined in NF conf)
export SINGULARITY_CACHEDIR="$PWD/singularity_cache"
mkdir -p $SINGULARITY_CACHEDIR

# tmp dir used by singularity when pulling images from dockerhub
export TMPDIR="$PWD/tmpdir"
mkdir -p $TMPDIR

# load conda env where nextflow is installed
eval "$(conda shell.bash hook)"
conda activate nextflow

nextflow run workflows/main.nf \
	 -c inputs.nf \
         --utilise_gpu \
         --singularity_use_pre_cached_images \
	 -profile lsf \
	 --nf_ci_loc $PWD \
	 -resume

