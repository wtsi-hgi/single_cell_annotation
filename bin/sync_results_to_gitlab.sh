#!/usr/bin/env bash
set -e
RESULTS_DIR=$1
cd $RESULTS_DIR

# source private gitlab token to clone repo, see instructions below
source /nfs/users/nfs_m/mercury/secrets/gitlab_scrna_cellranger_tokens.sh
export REPO_HTTPS=gitlab.internal.sanger.ac.uk/hgi-projects/scrna_cellranger.git
export REPO_DIR=scrna_cellranger
if [ ! -d "${REPO_DIR}" ]; then
  git clone https://${CI_DEPLOY_USER}:${CI_DEPLOY_TOKEN}@${REPO_HTTPS} ${REPO_DIR} > /dev/null 2>&1
fi

cd ${REPO_DIR}
git pull
mkdir -p ./sync_status/keras_celltypes

FULL_PATH_RESULTS=$(dirname $PWD)
echo \`"${FULL_PATH_RESULTS}/\`" >| ./sync_status/keras_celltypes/README.md
git add ./sync_status/keras_celltypes/README.md

# list input samples
mkdir -p ./sync_status/keras_celltypes/inputs
cp ../*.csv ./sync_status/keras_celltypes/inputs/
git add ./sync_status/keras_celltypes/inputs/*.csv

# add Nextflow html report and tasks trace
cp ../../reports/trace.txt ./sync_status/keras_celltypes/
git add ./sync_status/keras_celltypes/trace.txt
cp ../../reports/timeline.html ./sync_status/keras_celltypes/
git add ./sync_status/keras_celltypes/timeline.html

# cp nf_scrna_qc inputs
cp ../*.nf_scrna_qc_cellbender_inputs.tsv ./sync_status/keras_celltypes/
git add ./sync_status/keras_celltypes/*.nf_scrna_qc_cellbender_inputs.tsv

# list all samples processed
mkdir -p ./sync_status/keras_celltypes/cellbender
mkdir -p ./sync_status/keras_celltypes/multiplets
mkdir -p ./sync_status/keras_celltypes/celltypes
mkdir -p ./sync_status/keras_celltypes/merge
find ../cellbender/3_preprocess_output -maxdepth 2 -mindepth 2 | sort -o ./sync_status/keras_celltypes/cellbender/processed_samples.txt
find ../multiplet/1_scrublet -maxdepth 2 -mindepth 2 | sort -o ./sync_status/keras_celltypes/multiplets/processed_samples.txt
find ../celltype_prediction/1_keras -maxdepth 2 -mindepth 2 | sort -o ./sync_status/keras_celltypes/celltypes/processed_samples.txt
find ../merge/1_merge -maxdepth 1 -mindepth 1 | sort -o ./sync_status/keras_celltypes/merge/processed_biopsy_types.txt
git add ./sync_status/keras_celltypes/cellbender/processed_samples.txt
git add ./sync_status/keras_celltypes/multiplets/processed_samples.txt
git add ./sync_status/keras_celltypes/celltypes/processed_samples.txt
git add ./sync_status/keras_celltypes/merge/processed_biopsy_types.txt

git commit -m "keras pipeline run complete"
git push


## uses secret script:
#!/usr/bin/env bash

# cf. token generated from https://gitlab.internal.sanger.ac.uk/hgi-projects/scrna_cellranger/-/settings/access_tokens
# cf. doc at https://docs.gitlab.com/ee/user/project/deploy_tokens/
#export CI_DEPLOY_USER=gn5
#export CI_DEPLOY_TOKEN=****secret gitlab token****

# to be used with

# source /path_to/gitlab_scrna_cellranger_tokens.sh # this script
# export REPO_HTTPS=gitlab.internal.sanger.ac.uk/hgi-projects/scrna_cellranger.git > /dev/null 2>&1
# export REPO_DIR=repo
# git clone https://${CI_DEPLOY_USER}:${CI_DEPLOY_TOKEN}@${REPO_HTTPS} ${REPO_DIR} > /dev/null 2>&1
