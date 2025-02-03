clear all;

oldfolder = '/bwdata/NSDData/Seohee/Orientation/prfsample/';
orifolder = '/bwdata/NSDData/Seohee/Orientation/prfsample_Ori/';
controlfolder = '/bwdata/NSDData/Seohee/Orientation/prfsample_Ori_control/';
saveFolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/MaxMin/';

for isub = 1:8
    isub
    for visualRegion = 1:4

        %% load prfSample
        prf_old = load(fullfile(oldfolder,['prfSampleStim_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']),'prfSampleLevOri', 'allImgs');
        prf_ori = load(fullfile(orifolder,['prfSampleStim_ori_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']),'prfSampleLevOri', 'allImgs');
        prf_control = load(fullfile(controlfolder,['prfSampleStim_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']),'prfSampleLevOri', 'allImgs');
        allImgs = prf_old.allImgs;
        indicesTop_old_control = []; 
        indicesBottom_old_control = [];
        indicesTop_old_ori = [];
        indicesBottom_old_ori = [];
        indicesTop_control_ori = [];
        indicesBottom_control_ori = [];

        for curRoi =1:size(prf_old.prfSampleLevOri,1)

            for curVox = 1: size(prf_old.prfSampleLevOri{curRoi},2)
                curprf_old = squeeze(prf_old.prfSampleLevOri{curRoi}(:,curVox,:,:));
                curprf_old = squeeze(mean(curprf_old, 2));
                curprf_control = squeeze(prf_control.prfSampleLevOri{curRoi}(:,curVox,:,:));
                curprf_control = squeeze(mean(curprf_control, 2));
                curprf_ori = squeeze(prf_ori.prfSampleLevOri{curRoi}(:,curVox,:,:));

                %% get circ mean
                theta = linspace(0,2*pi,8+1);%for circular calculation
                theta = theta(1:end-1);
                for iimg=1:length(curprf_old)
                    prefAngle_old(iimg,1) = circ_mean(theta',curprf_old(iimg,:)');
                end
                prefAngle_old = mod(prefAngle_old,2*pi);%from [-pi, pi] to [0 2pi]
                prefAngle_old = prefAngle_old./2;%range 0 to pi.

                for iimg=1:length(curprf_control)
                    prefAngle_control(iimg,1) = circ_mean(theta',curprf_control(iimg,:)');
                end
                prefAngle_control = mod(prefAngle_control,2*pi);%from [-pi, pi] to [0 2pi]
                prefAngle_control = prefAngle_control./2;%range 0 to pi.

                for iimg=1:length(curprf_ori)
                    prefAngle_ori(iimg,1) = circ_mean(theta',curprf_ori(iimg,:)');
                end
                prefAngle_ori = mod(prefAngle_ori,2*pi);%from [-pi, pi] to [0 2pi]
                prefAngle_ori = prefAngle_ori./2;%range 0 to pi.

                %% Compute the absolute difference
                angleDiff_old_control = zeros(size(prefAngle_old));
                angleDiff_old_ori = zeros(size(prefAngle_old));
                angleDiff_control_ori = zeros(size(prefAngle_old));

                % Compute the absolute differences and wrap to [0, pi/2]
                % 1. prefAngle_old vs prefAngle_control
                angleDiff_old_control = abs(prefAngle_old - prefAngle_control);
                angleDiff_old_control = mod(angleDiff_old_control, pi);
                angleDiff_old_control(angleDiff_old_control > pi/2) = pi - angleDiff_old_control(angleDiff_old_control > pi/2);

                % 2. prefAngle_old vs prefAngle_ori
                angleDiff_old_ori = abs(prefAngle_old - prefAngle_ori);
                angleDiff_old_ori = mod(angleDiff_old_ori, pi);
                angleDiff_old_ori(angleDiff_old_ori > pi/2) = pi - angleDiff_old_ori(angleDiff_old_ori > pi/2);

                % 3. prefAngle_control vs prefAngle_ori
                angleDiff_control_ori = abs(prefAngle_control - prefAngle_ori);
                angleDiff_control_ori = mod(angleDiff_control_ori, pi);
                angleDiff_control_ori(angleDiff_control_ori > pi/2) = pi - angleDiff_control_ori(angleDiff_control_ori > pi/2);

                % calculate the median of each angleDiff
                median_old_control = median(angleDiff_old_control);
                median_old_ori = median(angleDiff_old_ori);
                median_control_ori = median(angleDiff_control_ori);

                % top/bottom halves
                if median_old_control ~= 0
                    topHalf_old_control = angleDiff_old_control >= median_old_control;
                    bottomHalf_old_control = angleDiff_old_control < median_old_control;
                else
                    topHalf_old_control  = zeros(length(curprf_old),1);
                    bottomHalf_old_control = zeros(length(curprf_old),1);
                end

                if median_old_ori ~= 0
                    topHalf_old_ori = angleDiff_old_ori >= median_old_ori;
                    bottomHalf_old_ori = angleDiff_old_ori < median_old_ori;
                else
                    topHalf_old_ori  = zeros(length(curprf_old),1);
                    bottomHalf_old_ori = zeros(length(curprf_old),1);
                end

                if median_control_ori ~= 0
                    topHalf_control_ori = angleDiff_control_ori >= median_control_ori;
                    bottomHalf_control_ori = angleDiff_control_ori < median_control_ori;
                else
                    topHalf_control_ori  = zeros(length(curprf_old),1);
                    bottomHalf_control_ori = zeros(length(curprf_old),1);
                end

                % Save the indices as needed
                indicesTop_old_control{curRoi}(:, curVox) = topHalf_old_control;
                indicesBottom_old_control{curRoi}(:, curVox) = bottomHalf_old_control;
                indicesTop_old_ori{curRoi}(:, curVox) = topHalf_old_ori;
                indicesBottom_old_ori{curRoi}(:, curVox) = bottomHalf_old_ori;
                indicesTop_control_ori{curRoi}(:, curVox) = topHalf_control_ori;
                indicesBottom_control_ori{curRoi}(:, curVox) = bottomHalf_control_ori;


            end
        end


        saveName = [saveFolder,'prfSampleIndices_v', num2str(visualRegion), '_sub', num2str(isub), '.mat'];
        save(saveName, 'indicesTop_old_control', 'indicesBottom_old_control',...
            'indicesTop_old_ori', 'indicesBottom_old_ori', ...
            'indicesTop_control_ori', 'indicesBottom_control_ori', 'allImgs');

    end
end