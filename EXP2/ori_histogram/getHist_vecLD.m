% photo
clear all;
rng(4228);

%% Add path
original_folder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli/orientationfilter'; % Original image file folder
addpath(original_folder);

save_folder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/ori_hist/vecLD';

allImgs = 1:73000;
imgSize = [512 512];

modelOriSum = struct;
ori_fields = {'ori1';'ori2';'ori3';'ori4';'ori5';'ori6';'ori7';'ori8'};

%% Read files
for imgIdx = 1:length(allImgs)
    fprintf('%d. ...\n', imgIdx);
    filename = ['oriImg' num2str(allImgs(imgIdx)) '.mat'];
    load(fullfile(original_folder, filename), 'modelOri');

    totalSum = sum(modelOri ,'all');
            for oriIdx = 1:8
                thisSum = sum(modelOri(oriIdx,:,:) ,'all');
                modelOriSum.(ori_fields{oriIdx})(imgIdx,1) = thisSum/totalSum;
            end
end

save(fullfile(save_folder, 'vecLD.mat'), 'modelOriSum');


