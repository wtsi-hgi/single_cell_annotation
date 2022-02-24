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
cp ../*.csv ./sync_status/keras_celltypes/
git add ./sync_status/keras_celltypes/*.csv
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
