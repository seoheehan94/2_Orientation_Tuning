#!/bin/bash


for cond in photoSP ldSP contour; do
  for sub in sub{1..8}; do
    for hemi in lh rh; do
        mri_vol2surf \
            --mov "../brainVolume/${cond}Brain_${sub}.nii" \
            --regheader "$sub" \
            --projfrac 0.5 \
            --interp nearest \
            --hemi "$hemi" \
            --o "${sub}/${cond}Brain_${sub}_${hemi}.mgh"
    done
  done
done

for cond in photoSP ldSP contour; do
  for sub in sub{1..8}; do
    for hemi in lh rh; do
        mri_vol2surf \
            --mov "../brainVolume/${cond}BrainR2_${sub}.nii" \
            --regheader "$sub" \
            --projfrac 0.5 \
            --interp nearest \
            --hemi "$hemi" \
            --o "${sub}/${cond}BrainR2_${sub}_${hemi}.mgh"
    done
  done
done

