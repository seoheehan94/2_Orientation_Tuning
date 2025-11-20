%%
clear all;

load('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/experiment/getImages/results_tableTotal.mat')
modelOriPhoto_folder = '/bwdata/NSDData/Seohee/stimuli/pyramid'; % Original image file folder
modelOriVecLD_folder = '/bwdata/NSDData/Seohee/stimuli/contour'; % Original image file folder
modelOriLD_folder = '/bwdata/NSDData/Seohee/stimuli/pyramid_LD'; % Original image file folder

%img size
topLeft = [77, 77];
imgWidth =360;
imgHeight = 360;

% theta
numAngles=8;
theta = linspace(0,2*pi,numAngles+1);%for circular calculation
theta = theta(1:end-1);

%colormap
numCols = 256;
base_cmap = turbo(numCols/2);  % half-length version
% Make it circular by mirroring
cmap = [base_cmap; base_cmap];

% Create a separate figure for the colorbar
% figCB = figure('Position',[100 100 150 400]); % width=100, height=400 (thicker)
% colormap(base_cmap);
% ax = axes('Position',[0.3 0.05 0.3 0.9]); % leave room for colorbar
% axis off; % no axes
% ax.CLim = [1 numCols/2];
% % Add colorbar with tick labels on the left
% cb = colorbar('Ticks', [1, numCols/4, numCols/2], ...
%               'TickLabels', {'0°','90°','180°'}, ...
%               'Location', 'West');  % left side
% pos = cb.Position;      % current position [x y width height]
% pos(3) = 0.1;           % increase width (thicker)
% cb.Position = pos;
% set(cb, 'Direction', 'reverse');
%exportgraphics(figCB, 'colourbar.png', 'BackgroundColor','none','Resolution',300);

% Create a circular mask
newRows = imgWidth/3;
newCols = imgHeight/3;
centerX = newCols / 2;
centerY = newRows / 2;
radius = newRows/2; % Set the radius of the circle
    
% Create a binary mask for the circle
[xGrid, yGrid] = meshgrid(1:imgWidth, 1:imgHeight);
circularMask = ((xGrid - centerX).^2 + (yGrid - centerY).^2) <= radius^2;

imgList = {'images01/img2114'};
imgNum = 1;
[~, mainFileName, ~] = fileparts(imgList{imgNum});
numStr = regexp(mainFileName, '\d+', 'match');  % returns a cell array
imgIdx = str2double(numStr{1});              % convert to number

%% Photograph
photoImg = imread(['/bwdata/NSDData/stimuli/',imgList{imgNum},'.png']);
photoImg = rgb2gray(photoImg);
photoImg = imresize(photoImg, [imgHeight, imgWidth]);
figure;imshow(photoImg)
imwrite(photoImg, [mainFileName, '_grey.png']); 

%% Line drawing 

load(['/bwdata/NSDData/stimuli/vecLD/',imgList{imgNum},'.mat']);
imgLD = renderLinedrawing(vecLD);
imgLD = squeeze(imgLD(:,:,1)); % use grayscale encoding
fig=figure; drawLinedrawing(vecLD); box on; set(gca, 'XTick', [],'YTick', []);      % (optional)
frame = getframe(gca);
LD_img = frame.cdata;
LD_img = imresize(LD_img, [imgHeight, imgWidth]);
imwrite(LD_img,strcat(mainFileName, '_LD.png'));

%% Contour - orientation map
filename = ['oriImg' num2str(imgIdx) '.mat'];
load(fullfile(modelOriVecLD_folder, filename), 'oriMap');
oriMap_cropped = imcrop(oriMap, [topLeft(1), topLeft(2), imgHeight-1, imgWidth-1]);
reverseIdx = (oriMap_cropped < 0); % Identify values modified in original operation
oriMap_cropped(reverseIdx) = oriMap_cropped(reverseIdx) + 180; % Add 180 back to reverse the previous operation

% colour scheme
cmap_with_nan = [cmap; 1 1 1];  % append white as last row
colorIdxImg = round((oriMap_cropped / 180) * (numCols/2)) + 1;
colorIdxImg(isnan(colorIdxImg)) = numCols + 1;

fig = figure; 
imagesc(colorIdxImg);
set(gca, 'XTick', [], 'YTick', []);
box on;colormap(cmap_with_nan);
caxis([1 numCols]); % ensures proper color scaling
axis image;
frame = getframe(gca);
LD_img = frame.cdata;
LD_img = imresize(LD_img, [imgHeight, imgWidth]);
% figure; imshow(LD_img)
imwrite(LD_img, [mainFileName, '_LD_coloured.png']);




%% Save Contour patch
rowIdx = strcmp(resultsTable.imgIdx_name, [mainFileName '.png']);

% Get the row(s)
matchingRow = resultsTable(rowIdx, :);

maxMask = circshift(circularMask, matchingRow.maxPos{1});
minMask = circshift(circularMask, matchingRow.minPos{1});

[row_max, col_max] = find(maxMask);
minRow_max = min(row_max);
maxRow_max = max(row_max);
minCol_max = min(col_max);
maxCol_max = max(col_max);
[row_min, col_min] = find(minMask);
minRow_min = min(row_min);
maxRow_min = max(row_min);
minCol_min = min(col_min);
maxCol_min = max(col_min);

% Create the alpha channel for transparency
alphaChannel_max = uint8(255 * maxMask); % Alpha channel (255 inside circle, 0 outside)
croppedAlphaChannel_max = alphaChannel_max(minRow_max:maxRow_max, minCol_max:maxCol_max, :);
alphaChannel_min = uint8(255 * minMask);
croppedAlphaChannel_min = alphaChannel_min(minRow_min:maxRow_min, minCol_min:maxCol_min, :);

for channel = 1:3
    currentChannel_max = LD_img(:,:,channel);  % Extract the channel
    currentChannel_max(~maxMask) = 255;  % Set background to white
    maskedImg_max(:,:,channel) = currentChannel_max; % Store it back in rgbImg
    currentChannel_min = LD_img(:,:,channel);  % Extract the channel
    currentChannel_min(~minMask) = 255;  % Set background to white
    maskedImg_min(:,:,channel) = currentChannel_min; % Store it back in rgbImg
end

maskedImg_max = uint8(maskedImg_max);
maskedImg_min = uint8(maskedImg_min);
croppedMaskedImg_max = maskedImg_max(minRow_max:maxRow_max, minCol_max:maxCol_max, :);
croppedMaskedImg_min = maskedImg_min(minRow_min:maxRow_min, minCol_min:maxCol_min, :);

figure; imagesc(croppedMaskedImg_max)
figure; imagesc(croppedMaskedImg_min)

outputFile_max = fullfile(['img',num2str(imgIdx),'_contour_max.png']);
imwrite(croppedMaskedImg_max, outputFile_max, 'Alpha', croppedAlphaChannel_max);
outputFile_min = fullfile(['img',num2str(imgIdx),'_contour_min.png']);
imwrite(croppedMaskedImg_min, outputFile_min, 'Alpha', croppedAlphaChannel_min);

%% photo steerable filter - orientation map

filename = ['pyrImg' num2str(imgIdx) '.mat'];
load(fullfile(modelOriPhoto_folder, filename), 'modelOri');
curModelOri = cat(4, modelOri{:});
curM = squeeze(mean(curModelOri,4));
curM_cropped = zeros(8,imgHeight, imgWidth);
for j = 1:8
    img = squeeze(curM(j, :, :));
    img_cropped = imcrop(img, [topLeft(1), topLeft(2), imgHeight-1, imgWidth-1]);
    curM_cropped(j, :, :) = img_cropped;
end
%curM_cropped=permute(curM_cropped, [2 3 1]);
curM_cropped = curM_cropped - min(curM_cropped(:));


photoRGB = cat(3, photoImg, photoImg, photoImg);
weight = 0.6;
photoRGB_d = im2double(photoRGB);

% Normalize curM_cropped to 0–1 for overlay (optional, if not already)
allVals = [];
for j = 1:8
    curImg = squeeze(curM_cropped(j,:,:));
    allVals = [allVals; curImg(:)];
end
clim = [min(allVals), max(allVals)];
imposeImg = zeros([8, size(photoRGB_d)]);

for j = 1:8
    curMap = squeeze(curM_cropped(j,:,:)); % pick one orientation map
    curMap = (curMap - min(allVals)) / (max(allVals) - min(allVals));

    % Blend red overlay
    imposeImg(j,:,:,1) = weight * photoRGB_d(:,:,1) + (weight) * curMap; % add red channel
    imposeImg(j,:,:,2) = weight * photoRGB_d(:,:,2); % keep grayscale G
    imposeImg(j,:,:,3) = weight * photoRGB_d(:,:,3); % keep grayscale B


end
%figure; imshow(curMap)
% figure; imshow(squeeze(imposeImg(3,:,:,:)))


% Plot all images with same color limits
fig = figure;
t = tiledlayout(2, 4,...
    'Padding', 'compact', ...   % minimal padding around edges
    'TileSpacing', 'tight');    % even smaller gaps between tiles
for j = 1:8
        nexttile
        imagesc(squeeze(imposeImg(j,:,:,:)), clim); % apply global color scale
        axis image off;  % remove ticks and preserve aspect ratio
        box on;
end
% Create colorbar for the entire layout
cb = colorbar;
cb.Layout.Tile = 'east';
redmap = [ones(256,1), linspace(1,0,256)', linspace(1,0,256)']; % white→red
colormap(fig, redmap);
cb.Ticks = [];
cb.TickLabels = {};
set(gcf, 'Position', [100 100 1000 550]);
exportgraphics(fig, [mainFileName, '_photoSP_meanSF_map_light.png'])
exportgraphics(fig, [mainFileName, '_photoSP_meanSF_map_light.pdf'])


%% Save photo steerable filter max/min patch

for j = 1:8
    for channel = 1:3
    currentChannel_max = squeeze(imposeImg(j, :,:,channel));  % Extract the channel
    currentChannel_max(~maxMask) = NaN;  % Set background to white
    maskedImg_max(j, :,:,channel) = currentChannel_max; % Store it back in rgbImg
    currentChannel_min = imposeImg(j, :,:,channel);  % Extract the channel
    currentChannel_min(~minMask) = NaN;  % Set background to white
    maskedImg_min(j,:,:,channel) = currentChannel_min; % Store it back in rgbImg
    end
end

croppedMaskedImg_max = maskedImg_max(:,minRow_max:maxRow_max, minCol_max:maxCol_max,:);
croppedMaskedImg_min = maskedImg_min(:,minRow_min:maxRow_min, minCol_min:maxCol_min,:);


% Plot all images with same color limits
fig = figure;
t = tiledlayout(2, 4,...
    'Padding', 'compact', ...   % minimal padding around edges
    'TileSpacing', 'tight');    % even smaller gaps between tiles
for j = 1:8
        nexttile
        imagesc(squeeze(croppedMaskedImg_max(j,:,:,:)), 'AlphaData', croppedAlphaChannel_max); % apply global color scale
        axis image off;  % remove ticks and preserve aspect ratio
        box on;
end
% cb = colorbar;
% cb.Layout.Tile = 'east';
set(gcf, 'Position', [100 100 1000 550]);
exportgraphics(fig, [mainFileName, '_photoSP_meanSF_max_light.png'])
exportgraphics(fig, [mainFileName, '_photoSP_meanSF_max_light.pdf'])


% Plot all images with same color limits
fig = figure;
t = tiledlayout(2, 4,...
    'Padding', 'compact', ...   % minimal padding around edges
    'TileSpacing', 'tight');    % even smaller gaps between tiles
for j = 1:8
        nexttile
        imagesc(squeeze(croppedMaskedImg_min(j,:,:,:)), 'AlphaData', croppedAlphaChannel_min); % apply global color scale
        axis image off;  % remove ticks and preserve aspect ratio
        box on;
end
% cb = colorbar;
% cb.Layout.Tile = 'east';
set(gcf, 'Position', [100 100 1000 550]);
exportgraphics(fig, [mainFileName, '_photoSP_meanSF_min_light.png'])
exportgraphics(fig, [mainFileName, '_photoSP_meanSF_min_light.pdf'])


