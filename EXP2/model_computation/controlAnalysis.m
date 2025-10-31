function controlAnalysis(isub,numregions, con)

nsdfolder = '/bwdata/NSDData/nsddata/experiments/nsd/';
photoSPfolder = '/bwdata/NSDData/Seohee/Orientation/prfsample_photoSP/';
ldSPfolder = '/bwdata/NSDData/Seohee/Orientation/prfsample_ldSP';
contourfolder = '/bwdata/NSDData/Seohee/Orientation/prfsample_contour';

if con == 1
    condition = 'ldSPresidual_photoSP';
    fittingMfolder = photoSPfolder;
    residualMfolder = ldSPfolder;
    bandpassStr = '_bandpass1to7';
elseif con == 2
    condition = 'photoSPresidual_ldSP';
    fittingMfolder = ldSPfolder;
    residualMfolder = photoSPfolder;
    bandpassStr = '_bandpass1to7';
elseif con == 3
    condition = 'contourresidual_ldSP';
    fittingMfolder = ldSPfolder;
    residualMfolder = contourfolder;
    bandpassStr = '_bandpass1to1';
end



nsdDesignFilename = fullfile(nsdfolder, 'nsd_expdesign.mat');
nsdDesign = load(nsdDesignFilename);
nsdDesign.masterordering;%for each of 30000 trials, what is corresponding image (out of 10000 images)


bandpass = 1; bandMin = 1; bandMax = 7;
numOrientations = 8;
nsessionsSub = [40 40 32 30 40 32 40 30];

for iregion=1:numregions
    visualRegion = iregion;%V1,V2,V3,V4
    nsessions=nsessionsSub(isub);


    fittingM = load(fullfile(fittingMfolder,['prfSampleStim_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']),'prfSampleLevOri','prfSampleLev',...
        'rois','allImgs','numLevels','numOrientations','interpImgSize','backgroundSize','pixPerDeg',...
        'roiPrf');  % fitting model - prfSampleLevOri
    residualM = load(fullfile(residualMfolder, ['regressPrfSplit' bandpassStr '_v' num2str(visualRegion) '_sub' num2str(isub)  '.mat']), ...
        'nsd','numLevels', 'numOrientations','rois','nvox','roiPrf','nsplits'); % residual of the model - nsd.voxOriResidual
    rois = fittingM.rois;

    if bandpass
        for roinum=1:length(fittingM.rois)
            fittingM.prfSampleLevOri{roinum} = fittingM.prfSampleLevOri{roinum}(:,:,bandMin:bandMax,:);
            fittingM.prfSampleLevOri{roinum} = sum(fittingM.prfSampleLevOri{roinum},3);
        end
        % numLevels = bandMax-bandMin+1;
        numLevels = 1;
    end
    nsplits = 2;

    subDesign = nsdDesign.subjectim(isub,nsdDesign.masterordering);%for each of 30000 trials, what is corresponding image (out of 73000 images)
    [imgTrials, imgNum] = ismember(subDesign, fittingM.allImgs);%logical array

    %if less than 40 sessions, only use image trials that were actually presented
    imgTrials = imgTrials(1:size(residualM.nsd.imgNum,2));
    imgNum = imgNum(1:size(residualM.nsd.imgNum,2));
    splitImgTrials = repmat(imgTrials,2,1);
    midImg = ceil(median(imgNum));
    splitImgTrials(1,imgNum<midImg) = zeros;
    splitImgTrials(2,imgNum>=midImg) = zeros;





    %% regress residuals of ldSP model with photoSP
    % voxOriResidual from ldSP
    % voxPrfOriSample from photoSp
    for roinum=1:length(rois)
        nvox(roinum) = size(fittingM.prfSampleLevOri{roinum},2);
        voxResidOriCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels*numOrientations+1);

        for isplit=1:nsplits
            imgTrials = splitImgTrials(isplit,:);
            numTrials = sum(imgTrials);
            for ivox=1:nvox(roinum)

                voxPrfOriSample = squeeze(fittingM.prfSampleLevOri{roinum}(imgNum(imgTrials>0),ivox,:,:));
                % voxPrfOriSample = squeeze(voxPrfOriSample(:,1,:));

                voxPrfOriSample = reshape(voxPrfOriSample,[],numLevels*numOrientations);

                %add constant predictor
                voxPrfOriSample(:,end+1) = ones;


                voxResidOriCoef{roinum}(isplit,ivox,:) = regress(squeeze(residualM.nsd.voxOriResidual{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:)))),voxPrfOriSample);

            end
        end
    end

    nsd.voxResidOriCoef = voxResidOriCoef;
    nsd.r2oriSplit = residualM.nsd.r2oriSplit;
    if length(rois)>1 %combine across ventral and dorsal ROIs
        oldNsd = nsd;
        nsd.r2oriSplit{1} = [];
        nsd.voxResidOriCoef{1} = [];
        for iroi=1:length(rois)
            nsd.r2oriSplit{1} = cat(2,nsd.r2oriSplit{1},oldNsd.r2oriSplit{iroi});
            nsd.voxResidOriCoef{1} = cat(2,nsd.voxResidOriCoef{1},oldNsd.voxResidOriCoef{iroi});
        end

        oldPrf = residualM.roiPrf; clear roiPrf;
        roiPrf{1}.ecc=[];
        roiPrf{1}.ang=[];
        roiPrf{1}.sz=[];
        %         roiPrf{1}.exponent=[];
        %         roiPrf{1}.gain=[];
        roiPrf{1}.r2=[];
        roiPrf{1}.x=[];
        roiPrf{1}.y=[];
        for iroi=1:length(rois)
            roiPrf{1}.ecc = cat(1,roiPrf{1}.ecc,oldPrf{iroi}.ecc);
            roiPrf{1}.ang = cat(1,roiPrf{1}.ang,oldPrf{iroi}.ang);
            roiPrf{1}.sz = cat(1,roiPrf{1}.sz,oldPrf{iroi}.sz);
            %             roiPrf{1}.exponent = cat(1,roiPrf{1}.exponent,oldPrf{iroi}.exponent);
            %             roiPrf{1}.gain = cat(1,roiPrf{1}.gain,oldPrf{iroi}.gain);
            roiPrf{1}.r2 = cat(1,roiPrf{1}.r2,oldPrf{iroi}.r2);
            roiPrf{1}.x = cat(1,roiPrf{1}.x,oldPrf{iroi}.x);
            roiPrf{1}.y = cat(1,roiPrf{1}.y,oldPrf{iroi}.y);
        end
        rois = 1;
    end


    nsd.voxResidOriCoef{1}(nsplits+1,:,:) = mean(nsd.voxResidOriCoef{1},1);
    nsd.r2oriSplit{1}(nsplits+1,:) = mean(nsd.r2oriSplit{1},1);
    nsplits = nsplits+1;

    %% get preferred ori
    %fullCoef = squeeze(nsd.voxResidOriCoef{iroi}(isplit,:,1:end-1));
    %residPrefOri{iroi}(isplit,:) = gratingPrefOri(fullCoef,gratings.modelOriEnergy);

    % get regress angle
    clear fullPrefOri oriDeviation vertDeviation cardDeviation fullCoef
    for  iroi=1:length(rois)%rois=1
        for isplit=1:nsplits
            prefAngle =[];
            fullCoef = squeeze(nsd.voxResidOriCoef{iroi}(isplit,:,1:end-1));
            numvox = size(fullCoef,1);
            % if con == 1 || con ==3
            coefMat = reshape(fullCoef,numvox,numLevels,numOrientations);%vox x levels x orientations
            coefMat_meansf = squeeze(mean(coefMat, 2));
            % elseif con ==2
            %     coefMat_meansf = fullCoef;
            % end
            coefMat_meansf = coefMat_meansf - min(coefMat_meansf,[],2);

            theta = linspace(0,2*pi,numOrientations+1);%for circular calculation
            theta = theta(1:end-1);
            for ivox=1:numvox
                prefAngle(ivox) = circ_mean(theta',coefMat_meansf(ivox,:)');
            end
            prefAngle = mod(prefAngle,2*pi);%from [-pi, pi] to [0 2pi]
            prefAngle = prefAngle./2;%range 0 to pi.
            % prefAngle = prefAngle / pi *180;

            % allMeanCoef.(thisfield) = prefAngle;
            fullPrefOri{iroi}(isplit,:) = prefAngle;
        end
    end



    allRoiPrf{iregion} = roiPrf{iroi};%iroi=1
    residOri{iregion} = fullPrefOri{iroi};
    roiNsdOriR2{iregion} = nsd.r2oriSplit{iroi};

    %% save

    save(['/bwdata/NSDData/Seohee/Orientation/' condition '_sub' num2str(isub) '.mat'],'allRoiPrf',...
        'residOri','roiNsdOriR2', 'nsplits');

end