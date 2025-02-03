#!/bin/bash

study='oriStudy'

for hemi in lh rh; do
  for smoothness in 10; do
    for meas in volume thickness; do
      mri_glmfit \
        --y "${hemi}.${meas}.${study}.${smoothness}.mgh" \
        --fsgd "FSGD/${study}.fsgd" \
        --C Contrasts/contrast.mtx  \
        --surf fsaverage "${hemi}"  \
        --cortex  \
        --glmdir "${hemi}.${meas}.${study}.${smoothness}.glmdir"
    done
  done
done
