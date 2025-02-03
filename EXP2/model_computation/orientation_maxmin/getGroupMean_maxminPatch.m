% getGroupMean.m
% For the group map, we computed the circular mean across subjects for each vertex in V1–V4,
% weighted by the full model R2 values.
clear all;

filedir = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/surfaceData/';
conditions = {'old', 'ori'};
imgTypes = {'Top', 'Bottom'};
pairTypes = {'old_control', 'old_ori', 'control_ori'};


addpath(genpath('/usr/local/freesurfer/7.4.1/matlab'));
[~,M,mr] = load_mgh('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/surfaceData/sub1/Bottom_control_ori_controlBrain_sub1_lh_fsaverage.mgh');
%save_mgh(vol, 'orig.stripped.mgz', M,mr);
%result = MRIread('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/surfaceData/sub1/angleBrain_sub1_lh_fsaverage.mgh');
% 50 percet cutoff for old R2 (participant averaged)
% cutoff_lh = 0.02;
% cutoff_rh = 0.02;

%Bottom_control_ori_control
for curimgtype = 1:2
    for curcond = 1:2
        %% get all subjects volume in radians
        for isub=1:8
            vol_lh(:,isub) = load_mgh([filedir, 'sub', num2str(isub), '/', ...
                imgTypes{curimgtype}, '_old_ori_', conditions{curcond}, 'Brain_sub', num2str(isub), '_lh_fsaverage.mgh']);
            vol_rh(:,isub) = load_mgh([filedir, 'sub', num2str(isub), '/', ...
                imgTypes{curimgtype}, '_old_ori_', conditions{curcond}, 'Brain_sub', num2str(isub), '_rh_fsaverage.mgh']);
        end

        % multiply by 2
        vol_lh = vol_lh*2;
        vol_rh = vol_rh*2;


        %% get all subjects full model R2 values

        for isub=1:8
            R2_lh(:,isub) = load_mgh([filedir, 'sub', num2str(isub), '/', ...
                imgTypes{curimgtype}, '_old_ori_', conditions{curcond}, 'BrainR2_sub', num2str(isub), '_lh_fsaverage.mgh']);
            R2_rh(:,isub) = load_mgh([filedir, 'sub', num2str(isub), '/', ...
                imgTypes{curimgtype}, '_old_ori_', conditions{curcond}, 'BrainR2_sub', num2str(isub), '_rh_fsaverage.mgh']);
        end

        R2_lh = R2_lh - min(R2_lh,[],2);
        R2_rh = R2_rh - min(R2_rh,[],2);
        %% weighted circ mean
        prefAngle_lh = [];
        prefAngle_rh = [];

        for ivox=1:size(vol_lh,1)
            prefAngle_lh(ivox) = circ_mean(vol_lh(ivox,:)',R2_lh(ivox,:)', 'omitnan');
        end
        for ivox=1:size(vol_rh,1)
            prefAngle_rh(ivox) = circ_mean(vol_rh(ivox,:)',R2_rh(ivox,:)', 'omitnan');
        end

        prefAngle_lh = mod(prefAngle_lh, 2*pi);
        prefAngle_lh = prefAngle_lh./2;%range 0 to pi.
        prefAngle_lh = prefAngle_lh';
        prefAngle_rh = mod(prefAngle_rh, 2*pi);
        prefAngle_rh = prefAngle_rh./2;%range 0 to pi.
        prefAngle_rh = prefAngle_rh';

        allPrefAngle_lh.(imgTypes{curimgtype}).(conditions{curcond}) = prefAngle_lh;
        allPrefAngle_rh.(imgTypes{curimgtype}).(conditions{curcond}) = prefAngle_rh;
    end

    %% only top 50% full model R2 of old
    % meanR2_lh = mean(R2_lh,2, "omitnan");
    % meanR2_rh = mean(R2_rh,2, "omitnan");
    % % cutoff_lh = prctile(meanR2_lh, 50);
    % % cutoff_rh = prctile(meanR2_rh, 50);
    % top50_prefAngle_lh = prefAngle_lh;
    % top50_prefAngle_lh(meanR2_lh < cutoff_lh) = 0;
    % top50_prefAngle_rh = prefAngle_rh;
    % top50_prefAngle_rh(meanR2_rh < cutoff_rh) = 0;
    %% save
    % fileName_lh = [filedir, condition{curcond}, 'Brain_groupmean_lh_fsaverage.mgh'];
    % fileName_rh = [filedir, condition{curcond}, 'Brain_groupmean_rh_fsaverage.mgh'];
    % save_mgh(prefAngle_lh, fileName_lh, M,mr);
    % save_mgh(prefAngle_rh, fileName_rh, M,mr);
    %
    % fileName_lh = [filedir, condition{curcond}, 'Brain_groupmean_top50_lh_fsaverage.mgh'];
    % fileName_rh = [filedir, condition{curcond}, 'Brain_groupmean_top50_rh_fsaverage.mgh'];
    % save_mgh(top50_prefAngle_lh, fileName_lh, M,mr);
    % save_mgh(top50_prefAngle_rh, fileName_rh, M,mr);

end

% get difference between conditions
for curimgtype = 1:2
    fieldsCon = fieldnames(allPrefAngle_lh.(imgTypes{curimgtype}));
    cur1_lh = allPrefAngle_lh.(imgTypes{curimgtype}).(fieldsCon{1});
    cur2_lh = allPrefAngle_lh.(imgTypes{curimgtype}).(fieldsCon{2});
    cur1_rh = allPrefAngle_rh.(imgTypes{curimgtype}).(fieldsCon{1});
    cur2_rh = allPrefAngle_rh.(imgTypes{curimgtype}).(fieldsCon{2});

    diff_lh = abs(cur1_lh - cur2_lh);
    diff_lh = mod(diff_lh,pi);
    diff_lh(diff_lh > pi/2) = pi - diff_lh(diff_lh > pi/2);
    allDiff_lh.(imgTypes{curimgtype}) = rad2deg(diff_lh);

    diff_rh = abs(cur1_rh - cur2_rh);
    diff_rh = mod(diff_rh,pi);
    diff_rh(diff_rh > pi/2) = pi - diff_rh(diff_rh > pi/2);
    allDiff_rh.(imgTypes{curimgtype}) = rad2deg(diff_rh);

    fileName_lh = [filedir, (imgTypes{curimgtype}), '_old_oriBrain_groupmean_lh_fsaverage.mgh'];
    fileName_rh = [filedir, (imgTypes{curimgtype}), '_old_oriBrain_groupmean_rh_fsaverage.mgh'];
    save_mgh(allDiff_lh.(imgTypes{curimgtype}), fileName_lh, M,mr);
    save_mgh(allDiff_rh.(imgTypes{curimgtype}), fileName_rh, M,mr);
end

mean(nonzeros(allDiff_lh.Top))
mean(nonzeros(allDiff_lh.Bottom))
[h,p,ci,stats] = ttest(nonzeros(allDiff_lh.Top),nonzeros(allDiff_lh.Bottom));
fprintf('p:%d. tstat:%d.\n', p, stats.tstat);

