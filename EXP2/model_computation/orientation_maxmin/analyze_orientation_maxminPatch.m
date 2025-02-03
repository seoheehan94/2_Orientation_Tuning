%analyze_orientation.m
% get voxel preference from model weights and create brainVolume
%   uses files created by: regressPrfSplit.m
%   creates files used by:
clear all;
savefolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/brainVolume_regress/maxmin';
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4'};
combinedRoiNames = {'V1','V2','V3','hV4'};
conditions = {'old', 'ori', 'control'};
prffolders = {'prfsample_maxmin', 'prfsample_maxmin_Ori', 'prfsample_maxmin_control'};
imgTypes = {'Top', 'Bottom'};
pairTypes = {'old_control', 'old_ori', 'control_ori'};


%% 1.  mean R2, AIC, BIC %%
allroiR2OriSplit=struct;
V1R2OriSplit=struct;
allroiaicOriSplit=struct;
V1aicOriSplit=struct;
allroibicOriSplit=struct;
V1bicOriSplit=struct;
columnNames = imgTypes;
for curpairtype = 1:3
    totalR2OriSplit = struct;
    totalaicOriSplit = struct;
    totalbicOriSplit = struct;

    if curpairtype == 1
        condition = [conditions(1), conditions(3)];
        prffolder = [prffolders(1), prffolders(3)];
    elseif curpairtype == 2
        condition = conditions(1:2); 
        prffolder = prffolders(1:2); 
    elseif curpairtype == 3
        condition = conditions(2:3);
        prffolder = prffolders(2:3); 
    end
    rowNames = condition;
    allroiR2OriSplit.(pairTypes{curpairtype}) = array2table(nan(numel(rowNames), numel(columnNames)), 'RowNames', rowNames, 'VariableNames', columnNames);
    V1R2OriSplit.(pairTypes{curpairtype}) = array2table(nan(numel(rowNames), numel(columnNames)), 'RowNames', rowNames, 'VariableNames', columnNames);
    allroiaicOriSplit.(pairTypes{curpairtype}) = array2table(nan(numel(rowNames), numel(columnNames)), 'RowNames', rowNames, 'VariableNames', columnNames);
    V1aicOriSplit.(pairTypes{curpairtype}) = array2table(nan(numel(rowNames), numel(columnNames)), 'RowNames', rowNames, 'VariableNames', columnNames);
    allroibicOriSplit.(pairTypes{curpairtype}) = array2table(nan(numel(rowNames), numel(columnNames)), 'RowNames', rowNames, 'VariableNames', columnNames);
    V1bicOriSplit.(pairTypes{curpairtype}) = array2table(nan(numel(rowNames), numel(columnNames)), 'RowNames', rowNames, 'VariableNames', columnNames);
    

    for curimgtype =1:2
        for con = 1:2
            totalR2OriSplit.(condition{con}) = {};
            totalaicOriSplit.(condition{con}) = {};
            totalbicOriSplit.(condition{con}) = {};
            for isub = 1:8

                curPrf = ['/bwdata/NSDData/Seohee/Orientation/', prffolder{con}, '/'];
                fprintf('isub:%d. con:%d. pairType:%s. imgType:%s. ...\n',isub,con, pairTypes{curpairtype}, imgTypes{curimgtype});
                load([curPrf 'indices' imgTypes{curimgtype} '_' pairTypes{curpairtype} 'voxModelPref_regress_sub' num2str(isub) '.mat']);

                %% total values of R2, aic, bic
                totalR2OriSplit.(condition{con}){end+1} = roiNsdOriR2;
                totalaicOriSplit.(condition{con}){end+1} = allaicOriSplit;
                totalbicOriSplit.(condition{con}){end+1} = allbicOriSplit;

            end
        end

        %% mean R2
        
        fieldsCon = fieldnames(totalR2OriSplit);
        for i = 1:numel(fieldsCon)
            curRoiR2OriSplit = [];
            curV1R2OriSplit = [];
            for j = 1:numel(totalR2OriSplit.(fieldsCon{i}))
                for k = 1:size(roiNsdOriR2,2)
                    curRoiR2OriSplit = [curRoiR2OriSplit, totalR2OriSplit.(fieldsCon{i}){j}{k}(3,:)];
                end

                curV1R2OriSplit = [curV1R2OriSplit, totalR2OriSplit.(fieldsCon{i}){j}{1}(3,:)];
            end
            % writematrix(curRoiR2OriSplit', ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/MaxMin/allroiR2', imgType{curimgtype}, '_', (fieldsCon{i}), '.csv']);
            % writematrix(curV1R2OriSplit', ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/MaxMin/V1R2', imgType{curimgtype}, '_', (fieldsCon{i}), '.csv']);

            allroiR2OriSplit.(pairTypes{curpairtype}){fieldsCon{i}, columnNames{curimgtype}} = mean(curRoiR2OriSplit, 'omitnan');
            V1R2OriSplit.(pairTypes{curpairtype}){fieldsCon{i}, columnNames{curimgtype}} = mean(curV1R2OriSplit,'omitnan');
        end

         % save(fullfile(['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/MaxMinPatch/', 'meanR2.mat']), "allroiR2OriSplit", "V1R2OriSplit");


        %% AIC/BIC
        
        fieldsCon = fieldnames(totalaicOriSplit);
        for i = 1:numel(fieldsCon)
            curRoiaicOriSplit = [];
            curV1aicOriSplit = [];
            curRoibicOriSplit = [];
            curV1bicOriSplit = [];
            for j = 1:numel(totalR2OriSplit.(fieldsCon{i}))
                for k = 1:size(roiNsdOriR2,2)
                    curRoiaicOriSplit = [curRoiaicOriSplit, totalaicOriSplit.(fieldsCon{i}){j}{k}(3,:)];
                    curRoibicOriSplit = [curRoibicOriSplit, totalbicOriSplit.(fieldsCon{i}){j}{k}(3,:)];
                end

                curV1aicOriSplit = [curV1aicOriSplit, totalaicOriSplit.(fieldsCon{i}){j}{1}(3,:)];
                curV1bicOriSplit = [curV1bicOriSplit, totalbicOriSplit.(fieldsCon{i}){j}{1}(3,:)];
            end
            % writematrix(curRoiaicOriSplit, ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/MaxMin/allroiR2', imgType{curimgtype}, '_', (fieldsCon{i}), '.csv']);
            % writematrix(curRoiaicOriSplit, ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/MaxMin/V1R2', imgType{curimgtype}, '_', (fieldsCon{i}), '.csv']);

            allroiaicOriSplit.(pairTypes{curpairtype}){fieldsCon{i}, columnNames{curimgtype}} = mean(curRoiaicOriSplit,"omitnan");
            V1aicOriSplit.(pairTypes{curpairtype}){fieldsCon{i}, columnNames{curimgtype}} = mean(curV1aicOriSplit,"omitnan");
            allroibicOriSplit.(pairTypes{curpairtype}){fieldsCon{i}, columnNames{curimgtype}} = mean(curRoibicOriSplit,"omitnan");
            V1bicOriSplit.(pairTypes{curpairtype}){fieldsCon{i}, columnNames{curimgtype}} = mean(curV1bicOriSplit,"omitnan");
        end
        %row: condition column: top/bottom
        % allroiaicOriSplit
        % V1aicOriSplit
        % allroibicOriSplit
        % V1bicOriSplit

        % save(fullfile(['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/MaxMinPatch/', 'meanAICBIC.mat']), "allroiaicOriSplit", "V1aicOriSplit", "allroibicOriSplit", "V1bicOriSplit");
    end
end



%% 2.  check differences in preferred orientation %%
%old vs control
pairType = 'old_ori';
nrois=4;
totalDiff = struct;
prffolder = prffolders(1:2);

for curimgtype =1:2
    for isub = 1:8
        allOri = cell(length(imgTypes), nrois);
        for con = 1:2
            curPrf = ['/bwdata/NSDData/Seohee/Orientation/', prffolder{con}, '/'];
            fprintf('isub:%d. con:%d. pairType:%s. imgType:%s. ...\n',isub,con, pairType, imgTypes{curimgtype});
            load([curPrf 'indices' imgTypes{curimgtype} '_' pairType 'voxModelPref_regress_sub' num2str(isub) '.mat']);
            for iroi = 1:nrois
                allOri{con, iroi} = [allOri{con, iroi} roiOri{iroi}];
            end
        end

        for iroi = 1:nrois
            diff{isub,iroi} = abs(allOri{1, iroi}(3,:) - allOri{2, iroi}(3,:));
            diff{isub,iroi} = mod(diff{isub,iroi},pi);
            diff{isub,iroi}(diff{isub,iroi} > pi/2) = pi - diff{isub,iroi}(diff{isub,iroi} > pi/2);
        end
    end
    totalDiff.(imgTypes{curimgtype}) = diff;
end

columnMeans = nan(2, 4);
for curimgtype=1:2
    for col = 1:size(totalDiff.(imgTypes{curimgtype}), 2)
        allRows = [];
        for row= 1:size(totalDiff.(imgTypes{curimgtype}), 1)
            allRows = [allRows (totalDiff.(imgTypes{curimgtype}){row, col})];
        end
        totalDiffbyROI.(imgTypes{curimgtype}){col} = allRows;
        % Compute the mean across the 8 rows
        columnMeans(curimgtype, col) = mean(nonzeros(allRows), "omitnan");
    end
end
rad2deg(columnMeans(1,:)-columnMeans(2,:))

[h,p,ci,stats] = ttest(nonzeros(totalDiffbyROI.Top{4}),nonzeros(totalDiffbyROI.Bottom{4}));
fprintf('p:%d. tstat:%d.\n', p, stats.tstat);
%% 2. make a brain volume %%

%% make a brain volume
pairTypes = {'old_control', 'old_ori', 'control_ori'};
prffolders = {'prfsample_maxmin', 'prfsample_maxmin_Ori', 'prfsample_maxmin_control'};
imgTypes = {'Top', 'Bottom'};
savefolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/brainVolume_regress/maxminPatch';
condition = {'old', 'ori', 'control'};

for isub = 1:8
    roifolder = ['/bwdata/NSDData/nsddata/ppdata/subj0' num2str(isub) '/func1pt8mm/'];
    visualRoisFile = fullfile(roifolder,'roi/prf-visualrois.nii.gz');%V1v, V1d, V2v, V2d, V3v, V3d, and hV4
    visRoiData = niftiread(visualRoisFile);

    ourBrain = visRoiData;
    ourBrain(ourBrain == 2) = 1;
    ourBrain(ourBrain == 3 | ourBrain == 4) = 2;
    ourBrain(ourBrain == 5 | ourBrain == 6) = 3;
    ourBrain(ourBrain == 7) = 4;
    for con = 2:3
        if con == 1
            pairType = pairTypes(1:2);
        elseif con == 2
            pairType = pairTypes(2:3);
        elseif con == 3
            pairType = [pairTypes(1), pairTypes(3)];
        end
        for curpairtype = 1:2
            for curimgtype =1:2
                curPrf = ['/bwdata/NSDData/Seohee/Orientation/', prffolders{con}, '/'];
                fprintf('isub:%d. con:%d. pairType:%s. imgType:%s. ...\n',isub,con, pairType{curpairtype}, imgTypes{curimgtype});
                load([curPrf 'indices' imgTypes{curimgtype} '_' pairType{curpairtype} 'voxModelPref_regress_sub' num2str(isub) '.mat']);


                % make a brain volume
                newBrain = ourBrain;
                newBrain(newBrain <= 0) = NaN;
                newBrain(newBrain > 0) = 0;
                for visualRegion = 1:4
                    curOurBrain = ourBrain;
                    % if visualRegion == 2
                    %     curOurBrain(visRoiData == 3 | visRoiData == 4) = 2;
                    % elseif visualRegion == 3
                    %     curOurBrain(visRoiData == 5 | visRoiData == 6) = 3;
                    % end
                    newBrain(curOurBrain == visualRegion) = roiOri{visualRegion}(3,:);
                end


            
                % R2
                r2Brain = ourBrain;
                r2Brain(ourBrain <= 0) = NaN;
                r2Brain(ourBrain > 0) = 0;
                
                for visualRegion = 1:4
                    thisfield = combinedRoiNames{visualRegion};
                    r2Brain(ourBrain == visualRegion) = roiNsdOriR2{visualRegion}(3,:);
                end

               
                %% save nifti
                % cd('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature/brainVolume_regress');
                % load(['oriBrain_sub', num2str(isub), '.mat']);
                info_old = niftiinfo(['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/brainVolume/betas_session01_sub', num2str(isub),'.nii.gz']);

                load([curPrf 'indices' imgTypes{curimgtype} '_' pairType{curpairtype} 'voxModelPref_regress_sub' num2str(isub) '.mat']);

                niftiwrite(newBrain,[savefolder, '/', imgTypes{curimgtype}, '_', pairType{curpairtype}, '_', condition{con}, 'Brain_sub', num2str(isub),'.nii']);
                info_new = niftiinfo([savefolder, '/', imgTypes{curimgtype}, '_', pairType{curpairtype}, '_', condition{con}, 'Brain_sub', num2str(isub),'.nii']);
                info_new.PixelDimensions = [1.8, 1.8, 1.8];
                info_new.TransformName = info_old.TransformName;
                info_new.SpatialDimension = info_old.SpatialDimension;
                info_new.Transform = info_old.Transform;
                info_new.Qfactor = info_old.Qfactor;
                info_new.AuxiliaryFile = info_old.AuxiliaryFile;
                info_new.raw.pixdim = info_old.raw.pixdim;
                info_new.raw.aux_file = info_old.raw.aux_file;
                info_new.raw.sform_code = info_old.raw.sform_code;
                info_new.raw.srow_x = info_old.raw.srow_x;
                info_new.raw.srow_y = info_old.raw.srow_y;
                info_new.raw.srow_z = info_old.raw.srow_z;
                niftiwrite(newBrain,[savefolder, '/', imgTypes{curimgtype}, '_', pairType{curpairtype}, '_', condition{con}, 'Brain_sub', num2str(isub),'.nii'], info_new);

                
                niftiwrite(r2Brain,[savefolder, '/', imgTypes{curimgtype}, '_', pairType{curpairtype}, '_', condition{con}, 'BrainR2_sub', num2str(isub),'.nii']);
                info_new = niftiinfo([savefolder, '/', imgTypes{curimgtype}, '_', pairType{curpairtype}, '_', condition{con}, 'BrainR2_sub', num2str(isub),'.nii']);
                info_new.PixelDimensions = [1.8, 1.8, 1.8];
                info_new.TransformName = info_old.TransformName;
                info_new.SpatialDimension = info_old.SpatialDimension;
                info_new.Transform = info_old.Transform;
                info_new.Qfactor = info_old.Qfactor;
                info_new.AuxiliaryFile = info_old.AuxiliaryFile;
                info_new.raw.pixdim = info_old.raw.pixdim;
                info_new.raw.aux_file = info_old.raw.aux_file;
                info_new.raw.sform_code = info_old.raw.sform_code;
                info_new.raw.srow_x = info_old.raw.srow_x;
                info_new.raw.srow_y = info_old.raw.srow_y;
                info_new.raw.srow_z = info_old.raw.srow_z;
                niftiwrite(r2Brain,[savefolder, '/', imgTypes{curimgtype}, '_', pairType{curpairtype}, '_', condition{con}, 'BrainR2_sub', num2str(isub),'.nii'],info_new);

            end
        end
    end
end



%% distribution of orientation in top/bottom

savefolder = '/bwdata/NSDData/Seohee/Orientation/prfsample_maxmin/';
variablename = {'indicesTop_old_ori', 'indicesBottom_old_ori'};
for curimgtype = 1:2
    for isub = 1:8
        load(fullfile(savefolder,[variablename{curimgtype},'regressPrfSplit' '_bandpass1to7_v1_sub' num2str(isub) '.mat']),'nsd');
        allsubOld.(variablename{curimgtype}){isub} = nsd.voxPrfOriSampleNormalized;
    end
end

combinedOld = struct;
for curimgtype = 1:2
    for isub = 1:8
        % Extract the inner cell
        roiCell = allsubOld.(variablename{curimgtype}){isub}; % 1x2 cell
        curRoi = [];
        for icell = 1:size(roiCell,2)
            curSecCell = roiCell{icell};
            curSecCell = permute(curSecCell,[3 2 4 1]);
            newarray = [curSecCell(:,:,:,1);curSecCell(:,:,:,1)];
            newarray = permute(newarray, [2 1 3]);
            isAllZero = all(newarray == 0, [1, 3]);
            newarray(:, isAllZero, :) = [];

            newarray_mean = squeeze(mean(newarray, 2));
            curRoi = [curRoi; newarray_mean];
        end
        combinedOld.(variablename{curimgtype}){isub} = curRoi;
    end
end
clearvars allsubOld

savefolder = '/bwdata/NSDData/Seohee/Orientation/prfsample_maxmin_Ori/';
variablename = {'indicesTop_old_ori', 'indicesBottom_old_ori'};
for curimgtype = 1:2
    for isub = 1:8
        load(fullfile(savefolder,[variablename{curimgtype},'regressPrfSplit' '_bandpass1to1_v1_sub' num2str(isub) '.mat']),'nsd');
        allsubOri.(variablename{curimgtype}){isub} = nsd.voxPrfOriSampleNormalized;
    end
end

combinedOri = struct;
for curimgtype = 1:2
    for isub = 1:8
        % Extract the inner cell
        roiCell = allsubOri.(variablename{curimgtype}){isub}; % 1x2 cell
        curRoi = [];
        for icell = 1:size(roiCell,2)
            curSecCell = roiCell{icell};
            curSecCell = permute(curSecCell,[3 2 4 1]);
            newarray = [curSecCell(:,:,:,1);curSecCell(:,:,:,1)];
            newarray = permute(newarray, [2 1 3]);
            isAllZero = all(newarray == 0, [1, 3]);
            newarray(:, isAllZero, :) = [];

            newarray_mean = squeeze(mean(newarray, 2));
            curRoi = [curRoi; newarray_mean];
        end
        combinedOri.(variablename{curimgtype}){isub} = curRoi;
    end
end
clearvars allsubOri

totalMean = struct;
rowNames = {'Top', 'Bottom'};
thisMean = array2table(nan(size(rowNames,2), 8), "RowNames",rowNames);
for curimgtype = 1:2
    allsub = [];
    for isub = 1:8
        allsub = [allsub; combinedOld.(variablename{curimgtype}){isub}];
    end
    thisMean{rowNames{curimgtype},:} = mean(allsub, 'omitnan');
end
totalMean.old = thisMean;

thisMean = array2table(nan(size(rowNames,2), 8), "RowNames",rowNames);
for curimgtype = 1:2
    allsub = [];
    for isub = 1:8
        allsub = [allsub; combinedOri.(variablename{curimgtype}){isub}];
    end
    thisMean{rowNames{curimgtype},:} = mean(allsub, 'omitnan');
end
totalMean.ori = thisMean;



% Plot the bar graph
data_old = table2array(totalMean.old);
figure;
bar(data_old', 'grouped'); % Transpose the data to align columns with x-axis
xlabel('Bins (Columns)');
ylabel('Value');
ylim([0,0.35]);
title('totalMean.old');
legend(totalMean.old.Properties.RowNames, 'Location', 'best');
xticks(1:size(data_old, 2)); % Set x-ticks to match bin numbers
xticklabels(totalMean.old.Properties.VariableNames); % Optional: use column names if available

data_ori = table2array(totalMean.ori);
figure;
bar(data_ori', 'grouped'); % Transpose the data to align columns with x-axis
xlabel('Bins (Columns)');
ylabel('Value');
ylim([0,0.35]);
title('totalMean.ori');
legend(totalMean.ori.Properties.RowNames, 'Location', 'best');
xticks(1:size(data_ori, 2)); % Set x-ticks to match bin numbers
xticklabels(totalMean.ori.Properties.VariableNames); % Optional: use column names if available


theta = linspace(0,2*pi,8+1);%for circular calculation
theta = theta(1:end-1);
prefAngle_oldTop = circ_mean(theta',totalMean.old{'Top',:}');
prefAngle_oldTop = mod(prefAngle_oldTop,2*pi);%from [-pi, pi] to [0 2pi]
prefAngle_oldTop = prefAngle_oldTop./2;%range 0 to pi.
prefAngle_oldBottom = circ_mean(theta',totalMean.old{'Bottom',:}');
prefAngle_oldBottom = mod(prefAngle_oldBottom,2*pi);%from [-pi, pi] to [0 2pi]
prefAngle_oldBottom = prefAngle_oldBottom./2;%range 0 to pi.
prefAngle_oriTop = circ_mean(theta',totalMean.ori{'Top',:}');
prefAngle_oriTop = mod(prefAngle_oriTop,2*pi);%from [-pi, pi] to [0 2pi]
prefAngle_oriTop = prefAngle_oriTop./2;%range 0 to pi.
prefAngle_oriBottom = circ_mean(theta',totalMean.ori{'Bottom',:}');
prefAngle_oriBottom = mod(prefAngle_oriBottom,2*pi);%from [-pi, pi] to [0 2pi]
prefAngle_oriBottom = prefAngle_oriBottom./2;%range 0 to pi.
rad2deg(prefAngle_oldTop)
rad2deg(prefAngle_oldBottom)
rad2deg(prefAngle_oriTop)
rad2deg(prefAngle_oriBottom)
