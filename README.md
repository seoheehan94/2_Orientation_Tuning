# Orientation Tuning in humans
This study investigates how humans perceive orientation in complex images by comparing different computational methods. Through behavioural experiments and statistical modelling, I show that people rely more on shapes and edges when judging orientation. Using neural data and an image-computable model, I demonstrate that contour-based orientation computations better align with human brain activity. These findings challenge filter-based models, highlighting the importance of edge-based processing in vision science.

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
   - saveAnalysisFiles.m: save R2/AIC/BIC and create brainVolumes
   - saveEccentricity.m: analyze, plot, and save R2 by eccentricity
   - controlAnalysis.m: control analysis where the second model is fit to the residuals of the first model
   - fig_scatterplots.m: scatterplots of R2/prf R2/prf eccentricity for V1
   - fig_preferredOri.m: plot preferred orientations for V1
   - fig_quantPlot.m: Quantitative plots of the deviation of preferred orientation from three idealized orientation preference maps
   - fig_modelOriResponses.m: sample orientation responses for model figure
   - fig_controlAnalysis.m: plot preferred orientations from controlAnalysis for V1 
 - ori_histogram: orientation distribution in images 
 - R2_analysis: explained variance of neural data
 - surfaceAnalysis: surface mapping of neural data
   - vol2surf.sh: convert preferred orientation and R2 from volume to surface
   - surf2surf.sh: convert preferred orientation and R2 from subject surface to fsaverage surface
   - getGroupMean.m: get the average of subjects
   - colorScale: colorScale for orientation preference in freeview
   - vol2surf_prf.sh: convert prf_eccentricity from volume to surface
   - surf2surf_prf.sh: convert prf_eccentricity from subject surface to fsaverage surface
   - script_visualRegion: get visual regions
 

