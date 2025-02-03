% getGroupMean.m
% For the group map, we computed the circular mean across subjects for each vertex in V1–V4, 
% weighted by the full model R2 values.
clear all;

filedir = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/surfaceData/';
condition = {'old', 'control','ori'};

addpath(genpath('/usr/local/freesurfer/7.4.1/matlab'));
[~,M,mr] = load_mgh('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/surfaceData/sub1/oldBrain_sub1_lh_fsaverage.mgh');

% 50 percet cutoff for old R2 (participant averaged)
cutoff_lh = 0.014;
cutoff_rh = 0.014;

% number of voxels survived
% numVoxel_lh = [];
% numVoxel_rh = [];

for curcond = 3
    %% get all subjects volume in radians
    % what to do with -1 and 0?
    for isub=1:8
        vol_lh(:,isub) = load_mgh([filedir, 'sub', num2str(isub), '/', condition{curcond}, 'Brain_sub', num2str(isub), '_lh_fsaverage.mgh']);
        vol_rh(:,isub) = load_mgh([filedir, 'sub', num2str(isub), '/', condition{curcond}, 'Brain_sub', num2str(isub), '_rh_fsaverage.mgh']);
    end
    
    % multiply by 2
    vol_lh = vol_lh*2;
    vol_rh = vol_rh*2;
    vol_lh_r = deg2rad(vol_lh);
    vol_rh_r = deg2rad(vol_rh);

    %% get all subjects full model R2 values
    % what to do with -1 and 0?
    for isub=1:8
        R2_lh(:,isub) = load_mgh([filedir, 'sub', num2str(isub), '/', condition{curcond}, 'BrainR2_sub', num2str(isub), '_lh_fsaverage.mgh']);
        R2_rh(:,isub) = load_mgh([filedir, 'sub', num2str(isub), '/', condition{curcond}, 'BrainR2_sub', num2str(isub), '_rh_fsaverage.mgh']);
    end


    %% weighted circ mean
    prefAngle_lh = [];
    prefAngle_rh = [];
    for ivox=1:size(vol_lh,1)
        prefAngle_lh(ivox) = circ_mean(vol_lh_r(ivox,:)',R2_lh(ivox,:)', 'omitnan');
    end
    for ivox=1:size(vol_rh,1)
        prefAngle_rh(ivox) = circ_mean(vol_rh_r(ivox,:)',R2_rh(ivox,:)', 'omitnan');
    end
   
    prefAngle_lh = mod(prefAngle_lh, 2*pi);
    prefAngle_lh = prefAngle_lh./2;%range 0 to pi.
    prefAngle_lh = prefAngle_lh/pi*180;
    prefAngle_lh = prefAngle_lh';
    prefAngle_rh = mod(prefAngle_rh, 2*pi);
    prefAngle_rh = prefAngle_rh./2;%range 0 to pi.
    prefAngle_rh = prefAngle_rh/pi*180;
    prefAngle_rh = prefAngle_rh';

    %% only top 50% full model R2 of old
    meanR2_lh = mean(R2_lh,2, "omitnan");
    meanR2_rh = mean(R2_rh,2, "omitnan");
    % cutoff_lh = prctile(meanR2_lh, 50);
    % cutoff_rh = prctile(meanR2_rh, 50);
    top50_prefAngle_lh = prefAngle_lh;
    top50_prefAngle_lh(meanR2_lh < cutoff_lh) = 0;
    top50_prefAngle_rh = prefAngle_rh;
    top50_prefAngle_rh(meanR2_rh < cutoff_rh) = 0;
    
    % meanR2_lh_save=(meanR2_lh(~isnan(meanR2_lh)));
    % meanR2_rh_save=(meanR2_rh(~isnan(meanR2_rh)));

    % writematrix(meanR2_lh_save, ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/subMeanR2_lh', condition{curcond}, '.csv']);
    % writematrix(meanR2_rh_save, ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/subMeanR_rh', condition{curcond}, '.csv']);

    % numVoxel_lh(curcond) = size(nonzeros(top50_prefAngle_lh),1);
    % numVoxel_rh(curcond) = size(nonzeros(top50_prefAngle_rh),1);
    %% save
    % fileName_lh = [filedir, condition{curcond}, 'Brain_sfmean_groupmean_lh_fsaverage.mgh'];
    % fileName_rh = [filedir, condition{curcond}, 'Brain_sfmean_groupmean_rh_fsaverage.mgh'];
    % save_mgh(prefAngle_lh, fileName_lh, M,mr);
    % save_mgh(prefAngle_rh, fileName_rh, M,mr);

    fileName_lh = [filedir, condition{curcond}, 'Brain_sfmean_groupmean_top50_lh_fsaverage.mgh'];
    fileName_rh = [filedir, condition{curcond}, 'Brain_sfmean_groupmean_top50_rh_fsaverage.mgh'];
    save_mgh(top50_prefAngle_lh, fileName_lh, M,mr);
    save_mgh(top50_prefAngle_rh, fileName_rh, M,mr);
end




% plot number of voxels

% figure;
% plot(numVoxel_lh,'-o', 'Color','#0072BD','LineWidth',2);
% hold on;
% plot(numVoxel_rh,'-o', 'Color','#F35872','LineWidth',2);
% xticks(1:3);
% xticklabels({'photograph-filter','line drawing-filter','contour'});
% xlabel('Methods');
% legend('left','right','FontSize',15,'FontName','Helvetica');
% ax = gca;
% ax.YAxis.FontSize = 15;
% ax.YAxis.FontName = 'Helvetica';
% title('Number of voxels above R2 0.02', 'FontSize',20,'FontName','Helvetica');
% box off;
% legend boxoff



%% get R2 difference
clear all;

filedir = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/surfaceData/';
condition = {'old', 'control','ori'};

addpath(genpath('/usr/local/freesurfer/7.4.1/matlab'));
[~,M,mr] = load_mgh('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/surfaceData/sub1/oldBrain_sub1_lh_fsaverage.mgh');

for isub=1:8
    for curcond = 1:3
        R2_lh(:,curcond) = load_mgh([filedir, 'sub', num2str(isub), '/', condition{curcond}, 'BrainR2_sub', num2str(isub), '_lh_fsaverage.mgh']);
        R2_rh(:,curcond) = load_mgh([filedir, 'sub', num2str(isub), '/', condition{curcond}, 'BrainR2_sub', num2str(isub), '_rh_fsaverage.mgh']); 
    end

    %old_ori
    R2diff_lh(:,isub) = R2_lh(:,3)- R2_lh(:,1);
    R2diff_rh(:,isub) = R2_rh(:,3)- R2_rh(:,1);

end

meanR2diff_lh = mean(R2diff_lh,2, "omitnan");
meanR2diff_rh = mean(R2diff_rh,2, "omitnan");

fileName_lh = [filedir, 'R2diff_oriold_lh_fsaverage.mgh'];
fileName_rh = [filedir, 'R2diff_oriold_rh_fsaverage.mgh'];
  
% save_mgh(meanR2diff_lh, fileName_lh, M,mr);
% save_mgh(meanR2diff_rh, fileName_rh, M,mr);


%% save by ROI
% Path to your .label file
visualregions = {'V1', 'V2', 'V3', 'V4'};
hemis = {'lh', 'rh'};
numVertices = 163842;
for visualregion = 1:4
    for hemi = 1:2
        labelFilePath = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/visualRegion/labels/';
        curfile = [labelFilePath, hemis{hemi}, '.', visualregions{visualregion}, '.label'];
        
        % Read the .label file
        labelData = read_label(curfile);

        % Initialize the overlay with NaN (no value) for all vertices
        overlay = nan(numVertices, 1);

        % Assign a value (e.g., 1.0) to vertices in the label file
        overlay(labelData.vertexIndex) = labelData.values;

        if hemi == 1
            curmeanR2diff = meanR2diff_lh;
        else 
            curmeanR2diff = meanR2diff_rh;
        end

        curmeanR2diff(isnan(overlay))=NaN;

        fileName = [filedir, 'R2diff_oriold_', visualregions{visualregion}, '_', hemis{hemi}, '_fsaverage.mgh'];
        save_mgh(curmeanR2diff, fileName, M,mr);

    end
end