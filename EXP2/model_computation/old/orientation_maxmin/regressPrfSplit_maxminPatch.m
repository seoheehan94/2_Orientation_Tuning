
% regressPrfSplit.m
%
% associated with the following publication: Roth, ZN, Kay, K, and Merriam, EP (2022).
% Massive natural scene sampling reveals reliable coarse-scale orientation tuning in human V1
% DOI:
%
%   usage: regressPrfSplit(1,1)
%   by: zvi roth
%   date: 7/29/2022
%   purpose: Perform linear regression on the response amplitudes for each voxel with filter output values as predictors
%   uses files created by: prfSampleModel.m, prfSampleModel_synth.m
%   creates files used by: getVoxPref.m

function regressPrfSplit_maxminPatch(isub,visualRegions,pairType, imgType)
addpath '/home/hanseohe/Documents/GitHub/nsdOtopy'
tic

%%%%%%%%%%%%%%%%%%%%%%%%
nsessionsSub = [40 40 32 30 40 32 40 30];
nsessions=nsessionsSub(isub);
nsplits=2;
bandpass = 1; bandMin = 1; bandMax = 7;

boxfolder = '/bwdata/NSDData/Seohee/Orientation/prfsample/';
savefolder = '/bwdata/NSDData/Seohee/Orientation/prfsample_maxmin/';
indicefolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/MaxMinPatch/';

betasfolder = ['/bwdata/NSDData/nsddata_betas/ppdata/subj0' num2str(isub) '/func1pt8mm/betas_fithrf_GLMdenoise_RR/'];
% stimfilename = fullfile(folder,'nsdsynthetic_colorstimuli_subj01.hdf5');
nsdfolder = '/bwdata/NSDData/nsddata/experiments/nsd/';
visualRoisfolder = ['/bwdata/NSDData/nsddata/ppdata/subj0' num2str(isub) '/func1pt8mm/'];

visualRoisFile = fullfile(visualRoisfolder,'roi/prf-visualrois.nii.gz');%V1v, V1d, V2v, V2d, V3v, V3d, and hV4
visRoiData = niftiread(visualRoisFile);
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4'};
visRoiData = visRoiData(:);

if strcmp(pairType, 'old_control')
    if strcmp(imgType, 'Top')
        variableName = 'indicesTop_old_control';
    elseif strcmp(imgType, 'Bottom')
        variableName = 'indicesBottom_old_control';
    else
        error('Invalid imgType. Must be "Top" or "Bottom".');
    end
elseif strcmp(pairType, 'old_ori')
    if strcmp(imgType, 'Top')
        variableName = 'indicesTop_old_ori';
    elseif strcmp(imgType, 'Bottom')
        variableName = 'indicesBottom_old_ori';
    else
        error('Invalid imgType. Must be "Top" or "Bottom".');
    end
elseif strcmp(pairType, 'control_ori')
    if strcmp(imgType, 'Top')
        variableName = 'indicesTop_control_ori';
    elseif strcmp(imgType, 'Bottom')
        variableName = 'indicesBottom_control_ori';
    else
        error('Invalid imgType. Must be "Top" or "Bottom".');
    end
else
    error('Invalid pairType. Must be "old_control" or "old_ori" or "control_ori".');
end

for visualRegion=visualRegions
    loadedData = load(fullfile(indicefolder,['prfSampleIndices_v', num2str(visualRegion), '_sub', num2str(isub), '.mat']), ...
        variableName);
    curIndices = loadedData.(variableName);
    load(fullfile(boxfolder,['prfSampleStim_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']),'prfSampleLevOri','prfSampleLev',...
        'rois','allImgs','numLevels','numOrientations','interpImgSize','backgroundSize','pixPerDeg',...
        'roiPrf');
    %if prf sampling was done with the nonlinear CSS prf, then we want to
    %define the weights for the constrained model as a sum of the
    %orientation model across orientations:
    % for roinum=1:length(rois)
    %     prfSampleLev{roinum} = squeeze(sum(prfSampleLevOri{roinum},4));
    % end
    
    if bandpass
        for roinum=1:length(rois)
            prfSampleLevOri{roinum} = prfSampleLevOri{roinum}(:,:,bandMin:bandMax,:);
            % prfSampleLev{roinum} = prfSampleLev{roinum}(:,:,bandMin:bandMax);
        end
        numLevels = bandMax-bandMin+1;
    end
    
    for roinum=1:length(rois); iroi = rois(roinum); roiBetas{roinum}=[]; end
    for isession=1:nsessions
        betasfilename = fullfile(betasfolder,['betas_session' num2str(isession,'%02.f') '.nii.gz']);
        betas = niftiread(betasfilename);
        betas = cast(betas,'double');
        betas = betas/300;
        betas=reshape(betas,[],size(betas,4));
        for roinum=1:length(rois)
            iroi = rois(roinum);
            roiBetas{roinum} = [roiBetas{roinum} betas(visRoiData==iroi,:)];
            roiInd{roinum} = find(visRoiData==iroi);
        end
    end
    
    nsdDesignFilename = fullfile(nsdfolder, 'nsd_expdesign.mat');
    nsdDesign = load(nsdDesignFilename);
    nsdDesign.masterordering;%for each of 30000 trials, what is corresponding image (out of 10000 images)
    subDesign = nsdDesign.subjectim(isub,nsdDesign.masterordering);%for each of 30000 trials, what is corresponding image (out of 73000 images)
    [imgTrials, imgNum] = ismember(subDesign, allImgs);%logical array
    
    %if less than 40 sessions, only use image trials that were actually presented
    imgTrials = imgTrials(1:size(roiBetas{roinum},2));
    imgNum = imgNum(1:size(roiBetas{roinum},2));

    for roinum=1:length(rois)
        for ivox = 1:size(curIndices{roinum},2)
            matchingIndices = find(curIndices{roinum}(:,ivox) == 1);
            logicalArray{roinum}(:,ivox) = ismember(imgNum, matchingIndices);
            totalMidImg{roinum}(ivox) = ceil(median(imgNum(logicalArray{roinum}(:,ivox)'>0)));
            totalNumTrials{roinum}(ivox) = size(imgNum(logicalArray{roinum}(:,ivox)'>0),2);
        end
         maxNumTrials{roinum} = ceil(max(totalNumTrials{roinum})/2);

    end
    
    % imgTrialSum = cumsum(imgTrials);
  
    % maxNumTrials = max(sum(splitImgTrials,2));

    r2 = cell(length(rois),1);
    r2ori = cell(length(rois),1);
    r2oriSplit = cell(length(rois),1);
    r2split = cell(length(rois),1);
    voxOriResidOriR2 = cell(length(rois),1);
    voxOriPredOriR2 = cell(length(rois),1);
    voxResidOriR2 = cell(length(rois),1);
    voxPredOriR2 = cell(length(rois),1);
    aicOriSplit = cell(length(rois),1);
    bicOriSplit = cell(length(rois),1);
    
    pearsonRori = cell(length(rois),1);
    pearsonR = cell(length(rois),1);
    
    for roinum=1:length(rois)
        nvox(roinum) = size(roiBetas{roinum},1);
        voxOriResidual{roinum} = NaN(nsplits, nvox(roinum),maxNumTrials{roinum});
        voxResidual{roinum} = NaN(nsplits, nvox(roinum),maxNumTrials{roinum});
        voxOriResidualSplit{roinum} = NaN(nsplits, nvox(roinum),maxNumTrials{roinum});
        voxResidualSplit{roinum} = NaN(nsplits, nvox(roinum),maxNumTrials{roinum});
        voxOriCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels*numOrientations+1);
        voxCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels+1);
        voxPredOriCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels*numOrientations+1);
        voxPredCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels+1);
        voxOriPredOriCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels*numOrientations+1);
        voxResidOriCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels*numOrientations+1);
        voxOriResidOriCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels*numOrientations+1);
        voxPrfOriSampleNormalized{roinum} = NaN(nsplits, nvox(roinum),maxNumTrials{roinum},numOrientations);

        %get model coefficients for each voxel, within each split
        for ivox=1:nvox(roinum)
            splitImgTrials = repmat(logicalArray{roinum}(:,ivox)',2,1);
            splitImgTrials(1,imgNum<totalMidImg{roinum}(ivox)) = zeros;
            splitImgTrials(2,imgNum>=totalMidImg{roinum}(ivox)) = zeros;
            for isplit=1:nsplits
                imgTrials = splitImgTrials(isplit,:);
                numTrials = sum(imgTrials);


                voxBetas = roiBetas{roinum}(ivox,imgTrials>0)';
                voxPrfSample = squeeze(prfSampleLev{roinum}(imgNum(imgTrials>0),ivox,:));
                %add constant predictor
                voxPrfSample(:,end+1) = ones;
                voxCoef{roinum}(isplit,ivox,:) = voxPrfSample\voxBetas;
                
                
                voxPrfOriSample = squeeze(prfSampleLevOri{roinum}(imgNum(imgTrials>0),ivox,:,:));
                
                voxPrfOriSample2 = squeeze(mean(voxPrfOriSample,2));
                rowSums = sum(voxPrfOriSample2, 2);
                rowSums(rowSums == 0) = NaN; % Replace zero sums with NaN to avoid division by zero
                voxPrfOriSampleNormalized{roinum}(isplit,ivox,1:numTrials,:) = voxPrfOriSample2 ./ rowSums;
                
                voxPrfOriSample = reshape(voxPrfOriSample,[],numLevels*numOrientations);

                %add constant predictor
                voxPrfOriSample(:,end+1) = ones;
                voxOriCoef{roinum}(isplit,ivox,:) = voxPrfOriSample\voxBetas;%check vox 144 in first ROI
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %regress vignetting predicted timecourse with orientation model
                voxPred = squeeze(voxCoef{roinum}(isplit,ivox,:))'*voxPrfSample';
                voxPredOriCoef{roinum}(isplit,ivox,:) = regress(voxPred',voxPrfOriSample);
                
                voxOriPred = squeeze(voxOriCoef{roinum}(isplit,ivox,:))'*voxPrfOriSample';
                voxOriPredOriCoef{roinum}(isplit,ivox,:) = regress(voxOriPred',voxPrfOriSample);
                
                %compute within-split residuals
                voxOriResidual{roinum}(isplit,ivox,1:numTrials) = voxBetas' - squeeze(voxOriCoef{roinum}(isplit,ivox,:))'*voxPrfOriSample';
                voxResidual{roinum}(isplit,ivox,1:numTrials) = voxBetas' - squeeze(voxCoef{roinum}(isplit,ivox,:))'*voxPrfSample';
                
                %regress residuals of vignetting model with full model
                voxResidOriCoef{roinum}(isplit,ivox,:) = regress(squeeze(voxResidual{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:)))),voxPrfOriSample);
                %regress residuals of orientation model with orientation model
                voxOriResidOriCoef{roinum}(isplit,ivox,:) = regress(squeeze(voxOriResidual{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:)))),voxPrfOriSample);
                
                %r2 within split
                r2{roinum}(isplit,ivox) = rsquared(voxResidual{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                r2ori{roinum}(isplit,ivox) = rsquared(voxOriResidual{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                
                %resid r2 within split
                %residuals of full orientation model:
                oriResidOriResidual = voxBetas' - squeeze(voxOriResidOriCoef{roinum}(isplit,ivox,:))'*voxPrfOriSample';
                voxOriResidOriR2{roinum}(isplit,ivox) = rsquared(oriResidOriResidual(1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                
                %prediction of full orientation model:
                oriPredOriResidual = voxBetas' - squeeze(voxOriPredOriCoef{roinum}(isplit,ivox,:))'*voxPrfOriSample';
                voxOriPredOriR2{roinum}(isplit,ivox) = rsquared(oriPredOriResidual(1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                
                
                %residuals of constrained model:
                residOriResidual = voxBetas' - squeeze(voxResidOriCoef{roinum}(isplit,ivox,:))'*voxPrfOriSample';
                voxResidOriR2{roinum}(isplit,ivox) = rsquared(residOriResidual(1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                
                
                %prediction of constrained model:
                predOriResidual = voxBetas' - squeeze(voxPredOriCoef{roinum}(isplit,ivox,:))'*voxPrfOriSample';
                voxPredOriR2{roinum}(isplit,ivox) = rsquared(predOriResidual(1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
            end
        end
        
        for ivox=1:nvox(roinum)
            splitImgTrials = repmat(logicalArray{roinum}(:,ivox)',2,1);
            splitImgTrials(1,imgNum<totalMidImg{roinum}(ivox)) = zeros;
            splitImgTrials(2,imgNum>=totalMidImg{roinum}(ivox)) = zeros;
            
            for isplit=1:nsplits
                imgTrials = splitImgTrials(isplit,:);
                numTrials = sum(imgTrials);
                
                voxBetas = roiBetas{roinum}(ivox,imgTrials>0)';
                voxPrfSample = squeeze(prfSampleLev{roinum}(imgNum(imgTrials>0),ivox,:));
                %add constant predictor
                voxPrfSample(:,end+1) = ones;
                voxPrfOriSample = squeeze(prfSampleLevOri{roinum}(imgNum(imgTrials>0),ivox,:,:));
                voxPrfOriSample = reshape(voxPrfOriSample,[],numLevels*numOrientations);
                
                %add constant predictor
                voxPrfOriSample(:,end+1) = ones;
                
                voxOriResidualSplit{roinum}(isplit,ivox,1:numTrials) = voxBetas' - squeeze(voxOriCoef{roinum}(nsplits-isplit+1,ivox,:))'*voxPrfOriSample';
                voxResidualSplit{roinum}(isplit,ivox,1:numTrials) = voxBetas' - squeeze(voxCoef{roinum}(nsplits-isplit+1,ivox,:))'*voxPrfSample';
                
                %r2 between splits
                r2split{roinum}(isplit,ivox) = rsquared(voxResidualSplit{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                r2oriSplit{roinum}(isplit,ivox) = rsquared(voxOriResidualSplit{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                aicOriSplit{roinum}(isplit,ivox) = computeAIC(voxOriResidualSplit{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:))), numTrials, size(voxOriCoef{roinum},3));
                bicOriSplit{roinum}(isplit,ivox) = computeBIC(voxOriResidualSplit{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:))), numTrials, size(voxOriCoef{roinum},3));
                
                %corr between splits
                % pearsonRori{roinum}(isplit,ivox) = corr(voxBetas,(squeeze(voxOriCoef{roinum}(nsplits-isplit+1,ivox,:))'*voxPrfOriSample')');
                % pearsonR{roinum}(isplit,ivox) = corr(voxBetas,(squeeze(voxCoef{roinum}(nsplits-isplit+1,ivox,:))'*voxPrfSample')');
                % 
            end
        end
        
        
    end
    
    nsd.voxResidual = voxResidual;
    nsd.voxOriResidual = voxOriResidual;
    nsd.voxResidualSplit = voxResidualSplit;
    nsd.voxOriResidualSplit = voxOriResidualSplit;
    nsd.r2 = r2;
    nsd.r2ori = r2ori;
    nsd.r2split = r2split;
    nsd.r2oriSplit = r2oriSplit;
    nsd.aicOriSplit = aicOriSplit;
    nsd.bicOriSplit = bicOriSplit;
    % nsd.pearsonRori = pearsonRori;
    % nsd.pearsonR = pearsonR;
    nsd.imgNum = imgNum;
    nsd.totalNumTrials = totalNumTrials;
    nsd.splitImgTrials = splitImgTrials;
    nsd.voxCoef = voxCoef;
    nsd.voxOriCoef = voxOriCoef;
    nsd.voxPredOriCoef = voxPredOriCoef;
    nsd.voxOriPredOriCoef = voxOriPredOriCoef;
    nsd.voxResidOriCoef = voxResidOriCoef;
    nsd.voxOriResidOriCoef = voxOriResidOriCoef;
    nsd.voxPredOriR2 = voxPredOriR2;
    nsd.voxOriPredOriR2 = voxOriPredOriR2;
    nsd.voxResidOriR2 = voxResidOriR2;
    nsd.voxOriResidOriR2 = voxOriResidOriR2;
    nsd.roiInd = roiInd;
    nsd.voxPrfOriSampleNormalized = voxPrfOriSampleNormalized;
    
    
    %% SYNTHETIC STIMULI
%     load(fullfile(prffolder,['prfSampleSynth_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']),'prfSampleLevOri','prfSampleLev',...
%         'rois','allImgs','numLevels','numOrientations','interpImgSize','backgroundSize','pixPerDeg','imgScaling',...
%         'roiPrf');
%     if bandpass
%         for roinum=1:length(rois)
%             prfSampleLevOri{roinum} = prfSampleLevOri{roinum}(:,:,bandMin:bandMax,:);
%             prfSampleLev{roinum} = prfSampleLev{roinum}(:,:,bandMin:bandMax);
%         end
%         numLevels = bandMax-bandMin+1;
%     end
%     
%     allImgs = 105:216;%1:220; ONLY GRATINGS AND SPIRALS
%     synthfolder = ['~/NSD/sub' num2str(isub) '_synth_func1pt8mm/'];
%     betasfilename = fullfile(synthfolder,'betas_nsdsynthetic.nii');
%     betas = niftiread(betasfilename);
%     betas = cast(betas,'double');
%     betas = betas/300;
%     
%     synthDesignFilename = fullfile(nsdfolder, 'nsdsynthetic_expdesign.mat');
%     synthDesign = load(synthDesignFilename);
%     synthDesign.masterordering;
%     
%     [imgTrials, imgNum] = ismember(synthDesign.masterordering, allImgs);%logical array
%     betas=reshape(betas,[],size(betas,4));
%     nimgs=length(allImgs);
%     imgTrialSum = cumsum(imgTrials);
%     totalTrials = max(imgTrialSum);%246
%     
%     clear roiBetas pearsonRori pearsonR voxOriResidual voxResidual
%     for roinum=1:length(rois)
%         iroi = rois(roinum);
%         roiBetas{roinum} = betas(visRoiData==iroi,:);
%     end
%     r2 = cell(length(rois),1);
%     r2ori = cell(length(rois),1);
%     voxOriResidual=cell(length(rois),1);
%     voxResidual=cell(length(rois),1);
%     voxOriCoef = cell(length(rois),1);
%     voxCoef = cell(length(rois),1);
%     for roinum=1:length(rois)
%         nvox(roinum) = size(roiBetas{roinum},1);
%         voxOriResidual{roinum} = zeros(nsplits, nvox(roinum),nimgs);
%         voxResidual{roinum} = zeros(nsplits, nvox(roinum),nimgs);
%         %compute R^2 for both models, within splits, and between splits
%         for isplit=1:nsplits
%             
%             %average across trials of each image
%             clear imgBetas imgPrfSample imgPrfOriSample
%             for iimg=1:nimgs
%                 imgBetas(:,iimg) = mean(roiBetas{roinum}(:,imgNum==iimg),2);
%             end
%             
%             for ivox=1:nvox(roinum)
%                 voxBetas = imgBetas(ivox,:)';
%                 voxPrfSample = squeeze(prfSampleLev{roinum}(:,ivox,:));
%                 voxPrfSample(:,end+1) = ones;
%                 voxCoef{roinum}(ivox,:) = voxPrfSample\voxBetas;
%                 
%                 voxPrfOriSample = squeeze(prfSampleLevOri{roinum}(:,ivox,:,:));
%                 voxPrfOriSample = reshape(voxPrfOriSample,[],numLevels*numOrientations);
%                 voxPrfOriSample(:,end+1) = ones;
%                 voxOriCoef{roinum}(ivox,:) = voxPrfOriSample\voxBetas;%check vox 144 in first ROI
%                 voxOriResidual{roinum}(isplit,ivox,:) = voxBetas' - squeeze(nsd.voxOriCoef{roinum}(isplit,ivox,:))'*voxPrfOriSample';
%                 voxResidual{roinum}(isplit,ivox,:) = voxBetas' - squeeze(nsd.voxCoef{roinum}(isplit,ivox,:))'*voxPrfSample';
%                 pearsonRori{roinum}(isplit,ivox) = corr(voxBetas,(squeeze(nsd.voxOriCoef{roinum}(isplit,ivox,:))'*voxPrfOriSample')');
%                 pearsonR{roinum}(isplit,ivox) = corr(voxBetas,(squeeze(nsd.voxCoef{roinum}(isplit,ivox,:))'*voxPrfSample')');
%             end
%         end
%     end
%     
%     synth.roiBetas = roiBetas;
%     synth.prfSampleLevOri = prfSampleLevOri;
%     synth.prfSampleLev = prfSampleLev;
%     synth.imgTrials = imgTrials;
%     synth.imgNum = imgNum;
%     synth.pearsonR = pearsonR;
%     synth.pearsonRori = pearsonRori;
%     synth.voxOriCoef = voxOriCoef;
%     synth.voxCoef = voxCoef;
%     synth.voxResidual = voxResidual;
%     synth.voxOriResidual = voxOriResidual;
    
    %% SAVE RESULTS
    bandpassStr = ['_bandpass' num2str(bandMin) 'to' num2str(bandMax)];
    save(fullfile(savefolder,[variableName,'regressPrfSplit' bandpassStr '_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']), ...
        'nsd',...
        'numLevels', 'numOrientations','rois','nvox','roiPrf','nsplits');
    toc
end

end

