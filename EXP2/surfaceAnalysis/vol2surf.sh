#!/bin/bash


#for cond in old ori control; do
for cond in old control; do
  for sub in sub{1..8}; do
    for hemi in lh rh; do
        mri_vol2surf \
            --mov "../brainVolume_regress/${cond}Brain_sfmean_${sub}.nii" \
            --regheader "$sub" \
            --projfrac 0.5 \
            --interp nearest \
            --hemi "$hemi" \
            --o "${sub}/${cond}Brain_sfmean_${sub}_${hemi}.mgh"
    done
  done
done

#for cond in old ori control; do
for cond in old control; do
  for sub in sub{1..8}; do
    for hemi in lh rh; do
        mri_vol2surf \
            --mov "../brainVolume_regress/${cond}BrainR2_sfmean_${sub}.nii" \
            --regheader "$sub" \
            --projfrac 0.5 \
            --interp nearest \
            --hemi "$hemi" \
            --o "${sub}/${cond}BrainR2_sfmean_${sub}_${hemi}.mgh"
    done
  done
done

