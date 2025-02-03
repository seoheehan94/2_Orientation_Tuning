% vecLD
clear all;
rng(4228);

%% Add path
original_folder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli/orientationfilter'; % Original image file folder
addpath(original_folder);

save_folder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/experiment/mean';

allImgs = 1:73000;
imgSize = [512 512];
topLeft = [77, 77];
width =360;
height = 360;
% modelOriSum = struct;
% sq_fields = {'sq1';'sq2';'sq3';'sq4';'sq5';'sq6';'sq7'};
% ori_fields = {'ori1';'ori2';'ori3';'ori4';'ori5';'ori6';'ori7';'ori8'};

meanOri = [];
% numAngles=8;
% binWidth2 = 90/8;

%% Read files
for imgIdx = 1:length(allImgs)
    fprintf('%d. ...\n', imgIdx);

    filename = ['oriImg' num2str(allImgs(imgIdx)) '.mat'];
    load(fullfile(original_folder, filename), 'oriMap');
    oriMap_cropped = imcrop(oriMap, [topLeft(1), topLeft(2), width-1, height-1]);

    reverseIdx = (oriMap_cropped < 0); % Identify values modified in original operation
    oriMap_cropped(reverseIdx) = oriMap_cropped(reverseIdx) + 180; % Add 180 back to reverse the previous operation
   
    allValues = oriMap_cropped(~isnan(oriMap_cropped)); 
    allValues_rad = 2*deg2rad(allValues);
    thisImgM = circ_mean(allValues_rad);
    thisImgM = mod(thisImgM,2*pi);%from [-pi, pi] to [0 2pi]
    thisImgM = thisImgM./2;%range 0 to pi.
    % thisImgM_deg = rad2deg(thisImgM);


    % curSum = [];
    % for t = 1:size(croppedArray,1)
    %     curSum=[curSum, sum(croppedArray(t,:,:), 'all')];
    % end
    % 
    % curMean = (1*curSum(1)+2*curSum(2)+3*curSum(3)+4*curSum(4)+5*curSum(5)+6*curSum(6)+7*curSum(7)+8*curSum(8))/sum(curSum);
    meanOri = [meanOri;thisImgM];

end

%%
save(fullfile(save_folder, 'mean_vecLD.mat'), 'meanOri');
