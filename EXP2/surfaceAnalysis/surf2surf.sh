#!/bin/bash


#for cond in old ori control; do
for cond in old control; do
  for sub in sub{1..8}; do
    # Loop through hemispheres
    for hemi in lh rh; do
        # Convert surface data to fsaverage surface
        mri_surf2surf \
            --srcsubject "$sub" \
            --sval "${sub}/${cond}Brain_sfmean_${sub}_${hemi}.mgh" \
            --trgsubject fsaverage \
            --mapmethod nnf \
            --tval "${sub}/${cond}Brain_sfmean_${sub}_${hemi}_fsaverage.mgh" \
            --hemi "$hemi"
    done
  done
done

#for cond in old ori control; do
for cond in old control; do
  for sub in sub{1..8}; do
    # Loop through hemispheres
    for hemi in lh rh; do
        # Convert surface data to fsaverage surface
        mri_surf2surf \
            --srcsubject "$sub" \
            --sval "${sub}/${cond}BrainR2_sfmean_${sub}_${hemi}.mgh" \
            --trgsubject fsaverage \
            --tval "${sub}/${cond}BrainR2_sfmean_${sub}_${hemi}_fsaverage.mgh" \
            --hemi "$hemi"
    done
  done
done
