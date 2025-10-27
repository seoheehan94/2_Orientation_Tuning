# Orientation Tuning in humans
This study investigates how humans perceive orientation in complex images by comparing texture-based and contour-based computational methods. Through behavioural experiments and statistical modelling, I show that people rely more on shapes and edges than textures when judging orientation. Using neural data and an image-computable model, I demonstrate that contour-based orientation computations better align with human brain activity. These findings challenge filter-based models, highlighting the importance of edge-based processing in vision science.

[Ppt link](https://drive.google.com/file/d/1S12H1B69XabA3PpYL-Oc38ubC4eHQlrM/view?usp=sharing)

---------
**EXP 1**
Orientation Judgement Experiment
 - getImages
   - getPatches.m: generate max/min patches used for experiment stimuli
     - Orientation computed by Steerable pyramid filter from Roth et al. (2022) (https://github.com/elimerriam/nsdOtopy)
     - Orientation computed by Contour from Walther et al. (2023) (https://github.com/bwlabToronto/MLV_toolbox)
   - makeGratings.m: generate grating images used for experiment stimuli
   - saveImages.m: saves figures for paper
 - runExperiment: psychophysics experiment files for exp1a and exp1b
 - dataAnalysis: experiment data and analysis for exp1a and exp1b
  
---------
**EXP 2**
Orientation Selectivity in Visual Cortex
Three models
1. Photo-Steerable Pyramid (photoSP) - replication of Roth et al. (2022)
   - https://github.com/elimerriam/nsdOtopy 
   - Roth, Z. N., Kay, K., & Merriam, E. P. (2022). Natural scene sampling reveals reliable coarse-scale orientation tuning in human V1. Nature communications, 13(1), 6469.
2. Line Drawing-Steerable Pyramid (ldSP)
3. Contour
   
 - model_computation: image-computable model with three different types of orientation computation
   - nsd_stim.m: run natural scene stimuli through the models, get energy responses
   - prfSampleModel.m: use voxel pRFs to sample model response outputs
   - regressPrfSplit.m: perform linear regression on the response amplitudes for each voxel with orientation response output values as predictors
   - getVoxProf.m: extract preferred orientation from regression weights for each voxel
   - runAllFiles.m: run prfSampleModel, regressPrfSplit, getVoxProf for all three models
 - R2_analysis: explained variance of neural data
 - surfaceAnalysis: surface mapping of neural data
 - ori_histogram: orientation distribution in images 

