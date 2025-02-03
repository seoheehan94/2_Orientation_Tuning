% photo
clear all;
rng(4228);

%% Add path
original_folder = '/bwdata/NSDData/Seohee/stimuli/pyramid'; % Original image file folder
addpath(original_folder);
addpath(genpath('/home/hanseohe/Documents/GitHub/stimulusVignetting'))
save_folder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/experiment/getImages';

allImgs = 1:73000;
imgSize = [512 512];
topLeft = [77, 77];
width =360;
height = 360;
cropRect = [topLeft(1), topLeft(2), width-1, height-1];

% modelOriSum = struct;
% sq_fields = {'sq1';'sq2';'sq3';'sq4';'sq5';'sq6';'sq7'};
% ori_fields = {'ori1';'ori2';'ori3';'ori4';'ori5';'ori6';'ori7';'ori8'};
% Check if the file exists
save_file = fullfile(save_folder, 'mean_photo.mat');
if exist(save_file, 'file')
    % Load the file if it exists
    load(save_file, 'meanOri');
    
    % Set imgIdx to start from length of meanOri + 1
    startIdx = length(meanOri) + 1;
else
    % Initialize meanOri if the file does not exist
    meanOri = [];
    
    % Set imgIdx to start from 1
    startIdx = 1;
end 
% load([save_folder, '/mean_photo.mat']);
% bin_centers = linspace(0, 180 - 22.5, 8);
numAngles=8;
theta = linspace(0,2*pi,numAngles+1);%for circular calculation
theta = theta(1:end-1);
%% Read files
for imgIdx = startIdx:25000
    fprintf('%d. ...\n', imgIdx);
    filename = ['pyrImg' num2str(allImgs(imgIdx)) '.mat'];
    load(fullfile(original_folder, filename), 'modelOri');

    curModelOri = cat(4, modelOri{:});
    curM = squeeze(mean(curModelOri,4));
    curM_cropped = zeros(8,width, height);
    for j = 1:8
        img = squeeze(curM(j, :, :));
        img_cropped = imcrop(img, [topLeft(1), topLeft(2), width-1, height-1]);
        curM_cropped(j, :, :) = img_cropped;
    end
    curMM_cropped = squeeze(mean(curM_cropped,[2 3]));
    % curM_cropped = arrayfun(@(j) imcrop(squeeze(curM(j, :, :)), cropRect), 1:8, 'UniformOutput', false);
    % curM_cropped = cat(3, curM_cropped{:});
    % curMM_cropped = squeeze(mean(curM_cropped,[1 2]));

    thisImgM = circ_mean(theta',curMM_cropped);
    thisImgM = mod(thisImgM,2*pi);%from [-pi, pi] to [0 2pi]
    thisImgM = thisImgM./2;%range 0 to pi.
    % thisImgM_deg = rad2deg(thisImgM);

    
    % curMM_cropped(curMM_cropped < 100) = NaN;
    % 
    % [m maxInd] = max(curM,[],1,"omitnan");
    % allNaN = isnan(squeeze(m)); 
    % maxInd = squeeze(maxInd);
    % maxInd(allNaN) = NaN;         % Set maxInd to NaN where all values are NaN
    % maxInd_cropped = imcrop(maxInd, [topLeft(1), topLeft(2), width-1, height-1]);
    % curMaxInd = maxInd_cropped(:);
    % 
    % curMaxInd = curMaxInd(~isnan(curMaxInd));
    % for i = 1:8
    %     curMaxInd_deg(curMaxInd == i) = bin_centers(i);
    % end
    % curMaxInd_deg=curMaxInd_deg';
    % curMaxInd_rad = deg2rad(curMaxInd_deg);
    % curMaxInd_rad2 =curMaxInd_rad*2;
    % BB = circ_mean(curMaxInd_rad2);
    % BB = mod(BB,2*pi);%from [-pi, pi] to [0 2pi]
    % BB = BB./2;%range 0 to pi.
    % BBB = rad2deg(BB);

    % curMean = mean(maxInd_cropped(:),"omitnan");
    meanOri = [meanOri;thisImgM];
end

%%
save(fullfile(save_folder, 'mean_photo.mat'), 'meanOri');

% m=squeeze(m);
% maxInd=squeeze(maxInd);
% figure;imagesc(squeeze(curM(8,:,:)))
% maxInd_cropped(isnan(maxInd_cropped)) = 0;
% figure;imagesc(maxInd_cropped)
% 
% figure;imagesc(squeeze(curM(1,:,:)))
