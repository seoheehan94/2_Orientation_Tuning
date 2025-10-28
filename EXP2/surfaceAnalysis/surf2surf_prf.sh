#!/bin/bash

for cond in old ori control; do
  for sub in sub{1..8}; do
    # Loop through hemispheres
    for hemi in lh rh; do
        # Convert surface data to fsaverage surface
        mri_surf2surf \
            --srcsubject "$sub" \
            --sval "${sub}/prf_eccentricity_${hemi}.mgh" \
            --trgsubject fsaverage \
            --tval "${sub}/prf_eccentricity_${hemi}_fsaverage.mgh" \
            --hemi "$hemi"
        mri_surf2surf \
            --srcsubject "$sub" \
            --sval "${sub}/prf_size_${hemi}.mgh" \
            --trgsubject fsaverage \
            --tval "${sub}/prf_size_${hemi}_fsaverage.mgh" \
            --hemi "$hemi"
    done
  done
done
