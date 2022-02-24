#!/usr/bin/env bash

# pull into singularity.cacheDir - which is defined in Nextflow config (lsf.conf)
cd /lustre/scratch118/humgen/resources/containers/

mkdir -p cache_dir
export SINGULARITY_CACHEDIR=$PWD/cache_dir

mkdir -p tmp_dir
export TMPDIR=$PWD/tmp_dir

# pull all singularity images used by pipeline - must match all images defined in Nextflow config lsf_tasks.conf
singularity pull wtsihgi-nf_cellbender_container-3cc9983.img docker://wtsihgi/nf_cellbender_container:3cc9983
singularity pull wtsihgi-nf_scrna_qc-6bb6af5.img docker://wtsihgi/nf_scrna_qc:6bb6af5
