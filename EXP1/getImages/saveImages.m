%%
clear all;

load('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/experiment/getImages/results_tableTotal.mat')
modelOriPhoto_folder = '/bwdata/NSDData/Seohee/stimuli/pyramid'; % Original image file folder
modelOriVecLD_folder = '/bwdata/NSDData/Seohee/stimuli/orientationfilter'; % Original image file folder
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
[xGrid, yGrid] = meshgrid(1:imgHeight, 1:imgWidth);
circularMask = ((xGrid - centerX).^2 + (yGrid - centerY).^2) <= radius^2;

imgList = {'images01/img2114','images01/img4017'};
imgNum = 2;
[~, mainFileName, ~] = fileparts(imgList{imgNum});
numStr = regexp(mainFileName, '\d+', 'match');  % returns a cell array
imgIdx = str2double(numStr{1});              % convert to number

%% Line drawing 

load(['/bwdata/NSDData/stimuli/vecLD/',imgList{imgNum},'.mat']);
imgLD = renderLinedrawing(vecLD);
imgLD = squeeze(imgLD(:,:,1)); % use grayscale encoding
drawLinedrawing(vecLD);
%imwrite(imgLD,strcat(mainFileName, '_LD.png'));

%% Contour - orientation map
filename = ['oriImg' num2str(imgIdx) '.mat'];
load(fullfile(modelOriVecLD_folder, filename), 'oriMap');
oriMap_cropped = imcrop(oriMap, [topLeft(1), topLeft(2), imgWidth-1, imgHeight-1]);
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
LD_img = imresize(LD_img, [imgWidth,imgHeight]);
% figure; imshow(LD_img)
imwrite(LD_img, [mainFileName, '_contour_map.png']);




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
curM_cropped = zeros(8,imgWidth, imgHeight);
for j = 1:8
    img = squeeze(curM(j, :, :));
    img_cropped = imcrop(img, [topLeft(1), topLeft(2), imgWidth-1, imgHeight-1]);
    curM_cropped(j, :, :) = img_cropped;
end
%curM_cropped=permute(curM_cropped, [2 3 1]);
curM_cropped = curM_cropped - min(curM_cropped(:));

% reshape data so each pixel is a column
[nTheta, h, w] = size(curM_cropped);
curM_reshaped = reshape(curM_cropped, [nTheta, h*w]);  % -> 8 x 129600

% compute weighted circular mean for each pixel
thisImgM = [];
for curP = 1: size(curM_reshaped,2)
    thisImgM(curP) = circ_mean(theta', curM_reshaped(:,curP));
end
    
thisImgM = reshape(thisImgM, h, w);
thisImgM = mod(thisImgM,2*pi);%from [-pi, pi] to [0 2pi]
thisImgM = thisImgM./2;%range 0 to pi.

thisImgM_deg = rad2deg(thisImgM);
thisImgM_deg = mod(thisImgM_deg, 180);

colorIdxImg = round(thisImgM_deg / 180 * (numCols/2) + 1);

% colour scheme
fig = figure; imagesc(colorIdxImg);
set(gca, 'XTick', [], 'YTick', []);
box on;colormap(cmap);
caxis([1 numCols]); % ensures proper color scaling
frame = getframe(gca);
photoSP_img = frame.cdata;
photoSP_img = imresize(photoSP_img, [imgWidth,imgHeight]);
% figure; imshow(LD_img)
imwrite(photoSP_img, [mainFileName, '_photoSP_map.png']);



%% Save photo steerable filter patch
curM_cropped
for channel = 1:8
    currentChannel_max = curM_cropped(channel, :,:);  % Extract the channel
    currentChannel_max(~maxMask) = 255;  % Set background to white
    maskedImg_max(channel, :,:) = currentChannel_max; % Store it back in rgbImg
    currentChannel_min = curM_cropped(channel, :,:);  % Extract the channel
    currentChannel_min(~minMask) = 255;  % Set background to white
    maskedImg_min(channel,:,:) = currentChannel_min; % Store it back in rgbImg
end
maskedImg_max = uint8(maskedImg_max);
maskedImg_min = uint8(maskedImg_min);
croppedMaskedImg_max = maskedImg_max(:,minRow_max:maxRow_max, minCol_max:maxCol_max);
croppedMaskedImg_min = maskedImg_min(:,minRow_min:maxRow_min, minCol_min:maxCol_min);

%% sanity check
maxMask = double(maxMask); % Convert to double
maxMask(maxMask == 0) = NaN;
maskedImg_max = double(maskedImg_max);
% Apply the shifted mask to the absolute image
maskedSection_photo = zeros(8,360,360);
for i = 1:8
    maskedSection_photo(i, :, :) = squeeze(maskedImg_max(i, :, :)) .* maxMask;
end

maskedSection_photo_reshaped = reshape(maskedSection_photo, [8, h*w]);  % -> 8 x 129600

photo_sum = sum(maskedSection_photo_reshaped,2, "omitnan");
curM = circ_mean(theta', photo_sum);
    
curM = mod(curM,2*pi);%from [-pi, pi] to [0 2pi]
curM = curM./2;%range 0 to pi.

curM_deg = rad2deg(curM);

figure; imagesc(squeeze(maskedSection_photo(1,:,:)));

%%
for channel = 1:3
    currentChannel_max = photoSP_img(:,:,channel);  % Extract the channel
    currentChannel_max(~maxMask) = 255;  % Set background to white
    maskedImg_max(:,:,channel) = currentChannel_max; % Store it back in rgbImg
    currentChannel_min = photoSP_img(:,:,channel);  % Extract the channel
    currentChannel_min(~minMask) = 255;  % Set background to white
    maskedImg_min(:,:,channel) = currentChannel_min; % Store it back in rgbImg
end

maskedImg_max = uint8(maskedImg_max);
maskedImg_min = uint8(maskedImg_min);
croppedMaskedImg_max = maskedImg_max(minRow_max:maxRow_max, minCol_max:maxCol_max, :);
croppedMaskedImg_min = maskedImg_min(minRow_min:maxRow_min, minCol_min:maxCol_min, :);

figure; imagesc(croppedMaskedImg_max)
figure; imagesc(croppedMaskedImg_min)

outputFile_max = fullfile(['img',num2str(imgIdx),'_photoSP_max.png']);
imwrite(croppedMaskedImg_max, outputFile_max, 'Alpha', croppedAlphaChannel_max);
outputFile_min = fullfile(['img',num2str(imgIdx),'_photoSP_min.png']);
imwrite(croppedMaskedImg_min, outputFile_min, 'Alpha', croppedAlphaChannel_min);


%% LD steerable filter - orientation map

filename = ['pyrImg' num2str(imgIdx) '.mat'];
load(fullfile(modelOriLD_folder, filename), 'modelOri');
curModelOri = cat(4, modelOri{:});
curM = squeeze(mean(curModelOri,4));
curM_cropped = zeros(8,imgWidth, imgHeight);
for j = 1:8
    img = squeeze(curM(j, :, :));
    img_cropped = imcrop(img, [topLeft(1), topLeft(2), imgWidth-1, imgHeight-1]);
    curM_cropped(j, :, :) = img_cropped;
end
%curM_cropped=permute(curM_cropped, [2 3 1]);
curM_cropped = curM_cropped - min(curM_cropped(:));

% reshape data so each pixel is a column
[nTheta, h, w] = size(curM_cropped);
curM_reshaped = reshape(curM_cropped, [nTheta, h*w]);  % -> 8 x 129600

% compute weighted circular mean for each pixel
thisImgM = [];
for curP = 1: size(curM_reshaped,2)
    thisImgM(curP) = circ_mean(theta', curM_reshaped(:,curP));
end
    
thisImgM = reshape(thisImgM, h, w);
thisImgM = mod(thisImgM,2*pi);%from [-pi, pi] to [0 2pi]
thisImgM = thisImgM./2;%range 0 to pi.

thisImgM_deg = rad2deg(thisImgM);
thisImgM_deg = mod(thisImgM_deg, 180);

colorIdxImg = round(thisImgM_deg / 180 * (numCols/2) + 1);

% colour scheme
fig = figure; imagesc(colorIdxImg);
set(gca, 'XTick', [], 'YTick', []);
box on;colormap(cmap);
caxis([1 numCols]); % ensures proper color scaling
frame = getframe(gca);
LDSP_img = frame.cdata;
LDSP_img = imresize(LDSP_img, [imgWidth,imgHeight]);
% figure; imshow(LD_img)
imwrite(LDSP_img, [mainFileName, '_LDSP_map.png']);


%% Save LD steerable filter patch

for channel = 1:3
    currentChannel_max = LDSP_img(:,:,channel);  % Extract the channel
    currentChannel_max(~maxMask) = 255;  % Set background to white
    maskedImg_max(:,:,channel) = currentChannel_max; % Store it back in rgbImg
    currentChannel_min = LDSP_img(:,:,channel);  % Extract the channel
    currentChannel_min(~minMask) = 255;  % Set background to white
    maskedImg_min(:,:,channel) = currentChannel_min; % Store it back in rgbImg
end

maskedImg_max = uint8(maskedImg_max);
maskedImg_min = uint8(maskedImg_min);
croppedMaskedImg_max = maskedImg_max(minRow_max:maxRow_max, minCol_max:maxCol_max, :);
croppedMaskedImg_min = maskedImg_min(minRow_min:maxRow_min, minCol_min:maxCol_min, :);

figure; imagesc(croppedMaskedImg_max)
figure; imagesc(croppedMaskedImg_min)

outputFile_max = fullfile(['img',num2str(imgIdx),'_LDSP_max.png']);
imwrite(croppedMaskedImg_max, outputFile_max, 'Alpha', croppedAlphaChannel_max);
outputFile_min = fullfile(['img',num2str(imgIdx),'_LDSP_min.png']);
imwrite(croppedMaskedImg_min, outputFile_min, 'Alpha', croppedAlphaChannel_min);


%% Contour - orientation map by orientation

filename = ['oriImg' num2str(imgIdx) '.mat'];
load(fullfile(modelOriVecLD_folder, filename), 'modelOri'); %→ 8 x W x H

curM_cropped = zeros(8,imgWidth, imgHeight);
for j = 1:size(modelOri,1)
    img = squeeze(modelOri(j, :, :));
    img_cropped = imcrop(img, [topLeft(1), topLeft(2), imgWidth-1, imgHeight-1]);
    curM_cropped(j, :, :) = img_cropped;
end

%thisImg = squeeze(curM_cropped(1,1,:,:));

allVals = [];
for j = 1:8
        curImg = squeeze(curM_cropped(j,:,:));
        allVals = [allVals; curImg(:)];
end
clim = [min(allVals), max(allVals)];

% Plot all images with same color limits
fig = figure;
t = tiledlayout(8, 1, ...
    'Padding', 'compact', ...   % minimal padding around edges
    'TileSpacing', 'tight');    % even smaller gaps between tiles
for j = 1:8
        nexttile
        imagesc(squeeze(curM_cropped(j,:,:)), clim); % apply global color scale
        axis image off;  % remove ticks and preserve aspect ratio
        box on;
end

set(gcf, 'Position', [100 100 900 1000]);
exportgraphics(fig, [mainFileName, '_contour_map.png'])
exportgraphics(fig, [mainFileName, '_contour_map.pdf'])


%% photo steerable filter - orientation map by spatical frequency & orientation

filename = ['pyrImg' num2str(imgIdx) '.mat'];
load(fullfile(modelOriPhoto_folder, filename), 'modelOri');
curModelOri = cat(4, modelOri{:}); % 8 orientation x W X H x 7 spatial frequency
curModelOri = permute(curModelOri, [1 4 2 3]); % reorder → 8 x 7 x W x H

curM_cropped = zeros(8,7, imgWidth, imgHeight);
for j = 1:size(curModelOri,1)
    for k = 1:size(curModelOri,2)
        img = squeeze(curModelOri(j, k, :, :));
        img_cropped = imcrop(img, [topLeft(1), topLeft(2), imgWidth-1, imgHeight-1]);
        curM_cropped(j,k, :, :) = img_cropped;
    end
end

%thisImg = squeeze(curM_cropped(1,1,:,:));


allVals = [];
for j = 1:8
    for k = 1:7
        curImg = squeeze(curM_cropped(j,k,:,:));
        allVals = [allVals; curImg(:)];
    end
end
clim = [min(allVals), max(allVals)];

% Plot all images with same color limits
fig = figure;
t = tiledlayout(8, 7, ...
    'Padding', 'compact', ...   % minimal padding around edges
    'TileSpacing', 'tight');    % even smaller gaps between tiles
for j = 1:8
    for k = 1:7
        nexttile
        imagesc(squeeze(curM_cropped(j,k,:,:)), clim); % apply global color scale
        axis image off;  % remove ticks and preserve aspect ratio
        box on;
    end
end
cb = colorbar;
cb.Layout.Tile = 'east';
set(gcf, 'Position', [100 100 900 1000]);
exportgraphics(fig, [mainFileName, '_photoSP_fullSF_map.png'])
exportgraphics(fig, [mainFileName, '_photoSP_fullSF_map.pdf'])



new = zeros(imgWidth, imgHeight);
for k = 1:7
        curImg = squeeze(curM_cropped(5,k,:,:));
        new = new + curImg;
end

figure; imagesc(new)


new = zeros(imgWidth, imgHeight);
for j = 1:8
        curImg = squeeze(curM_cropped(j,5,:,:));
        new = new + curImg;
end

figure; imagesc(new)
%% Lined drawing steerable filter - orientation map by spatical frequency & orientation

filename = ['pyrImg' num2str(imgIdx) '.mat'];
load(fullfile(modelOriLD_folder, filename), 'modelOri');
curModelOri = cat(4, modelOri{:}); % 8 orientation x W X H x 7 spatial frequency
curModelOri = permute(curModelOri, [1 4 2 3]); % reorder → 8 x 7 x W x H

curM_cropped = zeros(8,7, imgWidth, imgHeight);
for j = 1:size(curModelOri,1)
    for k = 1:size(curModelOri,2)
        img = squeeze(curModelOri(j, k, :, :));
        img_cropped = imcrop(img, [topLeft(1), topLeft(2), imgWidth-1, imgHeight-1]);
        curM_cropped(j,k, :, :) = img_cropped;
    end
end

%thisImg = squeeze(curM_cropped(1,1,:,:));


allVals = [];
for j = 1:8
    for k = 1:7
        curImg = squeeze(curM_cropped(j,k,:,:));
        allVals = [allVals; curImg(:)];
    end
end
clim = [min(allVals), max(allVals)];

% Plot all images with same color limits
fig = figure;
t = tiledlayout(8, 7, ...
    'Padding', 'compact', ...   % minimal padding around edges
    'TileSpacing', 'tight');    % even smaller gaps between tiles
for j = 1:8
    for k = 1:7
        nexttile
        imagesc(squeeze(curM_cropped(j,k,:,:)), clim); % apply global color scale
        axis image off;  % remove ticks and preserve aspect ratio
        box on;
    end
end
cb = colorbar;
cb.Layout.Tile = 'east';
set(gcf, 'Position', [100 100 900 1000]);
exportgraphics(fig, [mainFileName, '_LDSP_fullSF_map.png'])
exportgraphics(fig, [mainFileName, '_LDSP_fullSF_map.pdf'])

