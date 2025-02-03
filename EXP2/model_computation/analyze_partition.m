clear all;

savefolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/partition/';
%% 1.  mean R2, AIC, BIC %%
totalR2_old = cell(4,1);
totalR2_unique_ori = cell(4,1);
totalR2_combined = cell(4,1);
totalsharedVariance = cell(4,1);

for visualRegion =1:4
    for isub = 1:8
        fprintf('isub:%d. V:%d. ...\n',isub,visualRegion);
        load([savefolder 'regressPartitionCrossV_oldori_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']);
        oldNsd = nsd;
            nsd.R2_model1 = [];
            nsd.R2_model2_orth = [];
            nsd.R2_combined = [];
            nsd.sharedVariance = [];
        if length(rois)>1
            for iroi=1:length(rois)
                nsd.R2_model1 = cat(2,nsd.R2_model1,oldNsd.R2_model1{iroi});
                nsd.R2_model2_orth = cat(2,nsd.R2_model2_orth,oldNsd.R2_model2_orth{iroi});
                nsd.R2_combined = cat(2,nsd.R2_combined,oldNsd.R2_combined{iroi});
                nsd.sharedVariance = cat(2,nsd.sharedVariance,oldNsd.sharedVariance{iroi});
            end
        else
            nsd.R2_model1 = oldNsd.R2_model1{1};
            nsd.R2_model2_orth = oldNsd.R2_model2_orth{1};
            nsd.R2_combined = oldNsd.R2_combined{1};
            nsd.sharedVariance = oldNsd.sharedVariance{1};
        end
        nsd.R2_model1(3,:,:) = mean(nsd.R2_model1,1);
        nsd.R2_model2_orth(3,:,:) = mean(nsd.R2_model2_orth,1);
        nsd.R2_combined(3,:,:) = mean(nsd.R2_combined,1);
        nsd.sharedVariance(3,:,:) = mean(nsd.sharedVariance,1);



        %% total values
        totalR2_old{visualRegion}{end+1} = nsd.R2_model1(3,:,:);
        totalR2_unique_ori{visualRegion}{end+1} = nsd.R2_model2_orth(3,:,:);
        totalR2_combined{visualRegion}{end+1} = nsd.R2_combined(3,:,:);
        totalsharedVariance{visualRegion}{end+1} = nsd.sharedVariance(3,:,:);
    end

end


%% average participants
curR2_old=[];
curR2_unique_ori=[];
curR2_combined=[];
cur_sharedVariance=[];

for visualRegion =1:4
    
    for isub = 1:8
        curR2_old = [curR2_old, totalR2_old{visualRegion}{isub}];
        curR2_unique_ori = [curR2_unique_ori, totalR2_unique_ori{visualRegion}{isub}];
        curR2_combined = [curR2_combined, totalR2_combined{visualRegion}{isub}];
        cur_sharedVariance = [cur_sharedVariance, totalsharedVariance{visualRegion}{isub}];

    end
    % writematrix(curRoiR2OriSplit', ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/MaxMin/allroiR2', imgType{curimgtype}, '_', (fieldsCon{i}), '.csv']);
    % writematrix(curV1R2OriSplit', ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/MaxMin/V1R2', imgType{curimgtype}, '_', (fieldsCon{i}), '.csv']);
    if visualRegion == 1
        V1R2_model1 = mean(curR2_old, 'omitnan');
        V1R2_unique_model2 = mean(curR2_unique_ori,'omitnan');
        V1R2_combined = mean(curR2_combined, 'omitnan');
        V1sharedVariance = mean(cur_sharedVariance, 'omitnan');
    end


end

allroiR2_model1 = mean(curR2_old, 'omitnan');
allroiR2_unique_2 = mean(curR2_unique_ori,'omitnan');
allroiR2_combined = mean(curR2_combined, 'omitnan');
allroisharedVariance = mean(cur_sharedVariance, 'omitnan');