% photo
clear all;
rng(4228);

%% Add path
original_folder = '/bwdata/NSDData/Seohee/stimuli/pyramid'; % Original image file folder
addpath(original_folder);

save_folder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/ori_hist/photo';

allImgs = 1:73000;
imgSize = [512 512];

modelOriSum = struct;
sq_fields = {'sq1';'sq2';'sq3';'sq4';'sq5';'sq6';'sq7'};
ori_fields = {'ori1';'ori2';'ori3';'ori4';'ori5';'ori6';'ori7';'ori8'};

%% Read files
for imgIdx = 1:length(allImgs)
    fprintf('%d. ...\n', imgIdx);
    filename = ['pyrImg' num2str(allImgs(imgIdx)) '.mat'];
    load(fullfile(original_folder, filename), 'modelOri');

    totalSum = sum([modelOri{:}] ,'all');
        for sfIdx = 1:7
            for oriIdx = 1:8
                thisSum = sum(modelOri{sfIdx}(oriIdx,:,:) ,'all');
                modelOriSum.(sq_fields{sfIdx}).(ori_fields{oriIdx})(imgIdx,1) = thisSum/totalSum;
            end
        end
end

save(fullfile(save_folder, 'photo.mat'), 'modelOriSum');


