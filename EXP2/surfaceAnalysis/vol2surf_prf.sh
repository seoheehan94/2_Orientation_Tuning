#!/bin/bash


 for sub in sub{1..8}; do
    sub_padded=$(printf "subj%02d" "${sub#sub}")

    for hemi in lh rh; do
        mri_vol2surf \
            --mov "/bwdata/NSDData/nsddata/ppdata/${sub_padded}/func1pt8mm/prf_eccentricity.nii.gz" \
            --regheader "$sub" \
            --projfrac 0.5 \
            --interp nearest \
            --hemi "$hemi" \
            --o "${sub}/prf_eccentricity_${hemi}.mgh"
        mri_vol2surf \
            --mov "/bwdata/NSDData/nsddata/ppdata/${sub_padded}/func1pt8mm/prf_size.nii.gz" \
            --regheader "$sub" \
            --projfrac 0.5 \
            --interp nearest \
            --hemi "$hemi" \
            --o "${sub}/prf_size_${hemi}.mgh"
    done
done



