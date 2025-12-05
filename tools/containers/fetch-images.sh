#!/bin/bash

module load apptainer

export PIPELINE_IMAGES=${PIPELINE_IMAGES:-/misc/bioinf-ops/shared/singularity}
mkdir -p $PIPELINE_IMAGES

for uri in $(cat required-images.txt)
do
    mangled=$(echo $uri | sed 's/docker:\/\///; s/:/./; s/$/.sif/')
    dest="${PIPELINE_IMAGES}/${mangled}"
    mkdir -p $(dirname ${dest})
    echo "${uri}    ${dest}"
    singularity pull --force --disable-cache ${dest} ${uri}
done
