#!/bin/bash

study='oriStudy'

for hemi in lh rh; do
  for smoothing in 10; do
    for meas in volume thickness; do
      mris_preproc --fsgd FSGD/"$study".fsgd \
        --cache-in "$meas".fwhm"$smoothing".fsaverage \
        --target fsaverage \
        --hemi "$hemi" \
        --out "$hemi"."$meas"."$study"."$smoothing".mgh
    done
  done
done
