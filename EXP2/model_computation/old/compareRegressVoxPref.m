clear all;

prffolder = ['/bwdata/NSDData/Seohee/Orientation/prfsample/'];
bandpass = 1; bandMin = 1; bandMax = 7;
bandpassStr = '';
if bandpass
    bandpassStr = ['_bandpass' num2str(bandMin) 'to' num2str(bandMax)];
end
%%
for isub =1:8
isub

    %roiOri
    load([prffolder 'voxModelPref_sub' num2str(isub) '.mat'],'roiOri');

    numregions = 4;
    totalDiff = cell(4,1);
    for iregion=1:numregions
        load(fullfile(prffolder,['regressPrfSplit' bandpassStr '_v' num2str(iregion) '_sub' num2str(isub)  '.mat']), ...
            'nsd',...
            'numLevels', 'numOrientations','rois','nvox','roiPrf','nsplits');
        if length(rois)>1 %combine across ventral and dorsal ROIs
            oldNsd = nsd;
            nsd.r2{1} = [];
            nsd.r2ori{1} = [];
            nsd.voxOriCoef{1} = [];
            nsd.roiInd{1} = [];

            for iroi=1:length(rois)
                nsd.r2{1} = cat(2,nsd.r2{1},oldNsd.r2{iroi});
                nsd.r2ori{1} = cat(2,nsd.r2ori{1},oldNsd.r2ori{iroi});
                nsd.voxOriCoef{1} = cat(2,nsd.voxOriCoef{1},oldNsd.voxOriCoef{iroi});
                nsd.roiInd{1} = cat(1,nsd.roiInd{1}, oldNsd.roiInd{iroi});

            end

            oldPrf = roiPrf; clear roiPrf;
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

        %% AVERAGE SPLITS
        nsd.voxOriCoef{1}(nsplits+1,:,:) = mean(nsd.voxOriCoef{1},1);
        nsd.r2{1}(nsplits+1,:) = mean(nsd.r2{1},1);
        nsd.r2ori{1}(nsplits+1,:) = mean(nsd.r2ori{1},1);

        %% get regress angle
        prefAngle =[];
        fullCoef = squeeze(nsd.voxOriCoef{1}(3,:,1:end-1));
        numvox = size(fullCoef,1);
        numLevels = 7;
        numOrientations = 8;
        coefMat = reshape(fullCoef,numvox,numLevels,numOrientations);%vox x levels x orientations
        coefMat_meansf = squeeze(mean(coefMat, 2));
        coefMat_meansf = coefMat_meansf - min(coefMat_meansf,[],2);
        theta = linspace(0,2*pi,numOrientations+1);%for circular calculation
        theta = theta(1:end-1);
        for ivox=1:numvox
            prefAngle(ivox) = circ_mean(theta',coefMat_meansf(ivox,:)');
        end
        prefAngle = mod(prefAngle,2*pi);%from [-pi, pi] to [0 2pi]
        prefAngle = prefAngle./2;%range 0 to pi.

        %% get grating angle
        gratingAngle = roiOri{iregion}(3,:);

        %% compute difference
        angleDiff = abs(prefAngle' - gratingAngle');
        angleDiff(angleDiff > pi/2) = pi - angleDiff(angleDiff > pi/2);
        totalAngle = [];
        totalAngle = [prefAngle', rad2deg(prefAngle'), gratingAngle', rad2deg(gratingAngle'), angleDiff, rad2deg(angleDiff)];
        totalAngle = array2table(totalAngle, "VariableNames",{'regress', 'grating', 'regress_deg', 'grating_deg', 'diff', 'diff_deg'});

        totalDiff{iregion} = totalAngle;
    end


    save_folder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/compareRegressVoxPref';

    save(fullfile(save_folder,['comparison_sub' num2str(isub) '.mat']), 'totalDiff');

end
%%
save_folder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/compareRegressVoxPref';
ratios = zeros(8, length(totalDiff)); % Assuming 8 subjects

for isub = 1:8
    % Load data for the current subject
    load(fullfile(save_folder, ['comparison_sub' num2str(isub) '.mat']), 'totalDiff');
    
    for iroi = 1:length(totalDiff)
        curRoi = totalDiff{iroi};
        
        % Calculate the ratio of rows where diff_deg > 10
        ratioAbove10 = sum(curRoi.diff_deg > 10) / height(curRoi);
        
        % Store the result in the ratios matrix
        ratios(isub, iroi) = ratioAbove10;
    end
end

    save(fullfile(save_folder,['above10Percentage.mat']), 'ratios');
