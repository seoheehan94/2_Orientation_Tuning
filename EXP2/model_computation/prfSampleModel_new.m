% prfSampleModel_new.m
%
% associated with the following publication: Roth, ZN, Kay, K, and Merriam, EP (2022).
% Massive natural scene sampling reveals reliable coarse-scale orientation tuning in human V1
% DOI:
%
%   usage: prfSampleModel(1,1)
%   by: zvi roth
%   date: 7/29/2022
%   purpose: use voxel pRFs to sample model energy outputs for natural
%   images stimuli
%   uses files created by: nsdStim.m
%   creates files used by: regressPrfSplit.m

function prfSampleModel_new(isub,visualRegion)
cd '/home/hanseohe/Documents/GitHub/nsdOtopy';
delete(gcp('nocreate'));
g=gcp
distcomp.feature( 'LocalUseMpiexec', false ); % https://www.mathworks.com/matlabcentral/answers/447051-starting-matlab-pool-hangs-in-2018b

nsdfolder = '/bwdata/NSDData/nsddata/experiments/nsd/';
nsdDesignFilename = fullfile(nsdfolder, 'nsd_expdesign.mat');
nsdDesign = load(nsdDesignFilename);
allImgs = nsdDesign.subjectim(isub,nsdDesign.masterordering);%indices of all 10000 images used for this subject
allImgs = unique(allImgs);


colors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250], [0.4940, 0.1840, 0.5560], ...
    [0.4660, 0.6740, 0.1880],[0.3010, 0.7450, 0.9330], [0.6350, 0.0780, 0.1840]	};

interpImgSize = 714;
backgroundSize = 1024;
imgScaling = 0.5;

bandwidth = 1;
numOrientations = 8;
dims = [512 512];
%numLevels = maxLevel(dims,bandwidth);
numLevels = 1;

switch visualRegion
    case 1
        rois=1:2;%ventral and dorsal V1
    case 2
        rois=3:4;%ventral and dorsal V2
    case 3
        rois=5:6;%ventral and dorsal V3
    case 4
        rois=7;%V4
end

pixPerDeg = imgScaling*714/8.4;%=85
x = -(backgroundSize*imgScaling)/2+0.5:(backgroundSize*imgScaling)/2-0.5;
y = -(backgroundSize*imgScaling)/2+0.5:(backgroundSize*imgScaling)/2-0.5;
[X,Y] = meshgrid(x,-y);%flip up-down

orifolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli/orientationfilter/';%to save model outputs
betasfolder = ['/bwdata/NSDData/nsddata/ppdata/subj0' num2str(isub) '/func1pt8mm/'];
angFile = fullfile(betasfolder,'prf_angle.nii.gz');
eccFile = fullfile(betasfolder,'prf_eccentricity.nii.gz');
sizeFile = fullfile(betasfolder,'prf_size.nii.gz');
r2file = fullfile(betasfolder,'prf_R2.nii.gz');

angData = niftiread(angFile);
eccData = niftiread(eccFile);
sizeData = niftiread(sizeFile);
r2Data = niftiread(r2file);

visualRoisFile = fullfile(betasfolder,'roi/prf-visualrois.nii.gz');%V1v, V1d, V2v, V2d, V3v, V3d, and hV4
visRoiData = niftiread(visualRoisFile);

tic
%%

roiData = visRoiData;

roiPrf = cell(length(rois),1);
prfSampleLevOri = cell(length(rois),1);


for roinum=1:length(rois)
    iroi = rois(roinum);
    ecc = eccData(roiData(:)==iroi);
    ang = angData(roiData(:)==iroi);
    sz = sizeData(roiData(:)==iroi);% The definition of pRF size is one standard deviation of a Gaussian that describes the response of the model to point stimuli. Note that this definition of pRF size takes into account the nonlinear summation behavior of the pRF model and is not the same as the ?sigma? parameter used in the pRF model.
    nvox = length(ecc);
    r2 = r2Data(roiData(:)==iroi);
    
    roiPrf{roinum}.ecc = ecc;
    roiPrf{roinum}.ang = ang;
    roiPrf{roinum}.sz = sz;
    roiPrf{roinum}.r2 = r2;
    
    sz = sz*pixPerDeg;
    [x0, y0] = pol2cart(ang*2*pi/360,ecc*pixPerDeg);
    
    roiPrf{roinum}.x = x0;
    roiPrf{roinum}.y = y0;
    
    prfSampleLevOriRoi = zeros(length(allImgs),nvox,numLevels, numOrientations);
    
    %loop through  images
    parfor iimg=1:length(allImgs)
        ['sub: ' num2str(isub) ', roi: ' num2str(iroi) ', image: ' num2str(iimg)]
        imgNum = allImgs(iimg);
        orifilename = ['oriImg' num2str(imgNum) '.mat'];
        data = load(fullfile(orifolder, orifilename),'modelOri');
        imgPrfSampleLevOri = zeros(nvox,numLevels,numOrientations);
        %loop through voxels
        for ivox=1:nvox
            %sample for each ori level (summed across orientations!)
            voxPrfSampleLevOri = zeros(numLevels,numOrientations);
            %create the pRF for this voxel
            G = exp(-((X-x0(ivox)).^2+(Y-y0(ivox)).^2)/(2*sz(ivox)^2));
            
            for ilev = 1:numLevels
                voxPrfSampleOri = zeros(numOrientations,1);
                for iori=1:numOrientations
                    temp = data.modelOri(iori,:,:);
                    voxPrfSampleOri(iori) = temp(:)'*G(:);
                    keyboard;
                end
                voxPrfSampleLevOri(ilev,:) = voxPrfSampleOri;
            end
            imgPrfSampleLevOri(ivox,:,:) = voxPrfSampleLevOri;
        end
        prfSampleLevOriRoi(iimg,:,:,:) = imgPrfSampleLevOri;
    end
    prfSampleLevOri{roinum} = prfSampleLevOriRoi;

    prffolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/prfsample_Ori/';
    save(fullfile(prffolder,['prfSampleStim_ori_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']),'prfSampleLevOri',...
        'rois','allImgs','numLevels','numOrientations','interpImgSize','backgroundSize','pixPerDeg',...
        'roiPrf','-v7.3');
end

% prffolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/prfsample/';
% save(fullfile(prffolder,['prfSampleStim_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']),'prfSampleLevOri','prfSampleLev',...
%     'rois','allImgs','numLevels','numOrientations','interpImgSize','backgroundSize','pixPerDeg',...
%     'roiPrf','-v7.3');
delete(g);

%% CSS pRF from Kay et al., J Neurophysiology, 2013
    function r = prfResp(S,x,y,x0,y0,sigma, n, g)
        G = exp(-((x-x0).^2+(y-y0).^2)/(2*sigma^2));
        % r = g*(sum(S(:).*G(:))^n);
        r = g*((S(:)'*G(:))^n);
    end
end


%%
%%figure; hold on; for ori=1:8;plot(squeeze(prfSampleLevOri{1,1}(1:100,23,1,ori)),'-','Color',colors(ori,:));end
