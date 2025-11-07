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
n = 256;
n1 = floor(n/2);     % white -> orange
n2 = n - n1;         % orange -> red

% segment 1: white (1,1,1) -> orange (1,0.6,0)
seg1 = [ones(n1,1), ...                 % R
    linspace(1,0.6,n1)', ...             % G
    linspace(1,0,n1)' ];                 % B

% segment 2: orange (1,0.6,0) -> red (1,0,0)
seg2 = [ones(n2,1), ...                  % R
    linspace(0.6,0,n2)', ...             % G
    zeros(n2,1) ];                       % B

% combine
cmap_wr = [seg1; seg2];


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

imgList = {'images01/img4017','images01/img5178','images01/img8687','images01/img8727'};
imgNum = 4;
[~, mainFileName, ~] = fileparts(imgList{imgNum});
numStr = regexp(mainFileName, '\d+', 'match');  % returns a cell array
imgIdx = str2double(numStr{1});              % convert to number

%% Photograph
photoImg = imread(['/bwdata/NSDData/stimuli/',imgList{imgNum},'.png']);
photoImg = rgb2gray(photoImg);
photoImg = imresize(photoImg, [imgWidth,imgHeight]);
figure;imshow(photoImg)
imwrite(photoImg, [mainFileName, '_grey.png']); 

%% Line drawing 

load(['/bwdata/NSDData/stimuli/vecLD/',imgList{imgNum},'.mat']);
imgLD = renderLinedrawing(vecLD);
imgLD = squeeze(imgLD(:,:,1)); % use grayscale encoding
fig=figure; drawLinedrawing(vecLD); box on; set(gca, 'XTick', [],'YTick', []);      % (optional)
frame = getframe(gca);
LD_img = frame.cdata;
LD_img = imresize(LD_img, [imgWidth,imgHeight]);
imwrite(LD_img,strcat(mainFileName, '_LD.png'));


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
    'TileSpacing', 'compact');    % even smaller gaps between tiles
for j = 1:8
        nexttile
        imagesc(squeeze(curM_cropped(j,:,:)), clim); % apply global color scale
        axis image;                   % keep correct aspect ratio, don't hide axes
        set(gca, 'XTick', [], ...     % remove X tick marks
            'YTick', []);      % (optional)
        box on;
end
redmap = [ones(256,1), linspace(1,0,256)', linspace(1,0,256)']; % white→red
colormap(fig, cmap_wr);
set(gcf, 'Position', [100 100 500 1000]);
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

% --- Apply root-law compression ---
allVals_log = allVals .^ 0.5;
clim = [min(allVals_log), max(allVals_log)];

% Plot all images with same color limits
fig = figure;
t = tiledlayout(8, 7, ...
    'Padding', 'compact', ...   % minimal padding around edges
    'TileSpacing', 'compact');    % even smaller gaps between tiles
for j = 1:8
    for k = 1:7
        nexttile
        curImg = squeeze(curM_cropped(j,k,:,:));
        rootImg = curImg.^ 0.5;     % apply to each image
        imagesc(rootImg, clim);
        axis image;
        set(gca, 'XTick', [], 'YTick', []);
        box on;
    end
end
colormap(fig, cmap_wr);
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Ticks = [];
cb.TickLabels = {};
set(gcf, 'Position', [100 100 900 1000]);
exportgraphics(fig, [mainFileName, '_photoSP_fullSF_map.png'])
exportgraphics(fig, [mainFileName, '_photoSP_fullSF_map.pdf'])


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

% thisImg = squeeze(curM_cropped(1,6,:,:));
% figure; imagesc(thisImg)

allVals = [];
for j = 1:8
    for k = 1:7
        curImg = squeeze(curM_cropped(j,k,:,:));
        allVals = [allVals; curImg(:)];
    end
end
% --- Apply root-law compression ---
allVals_log = allVals .^ 0.5;
clim = [min(allVals_log), max(allVals_log)];


% Plot all images with same color limits
fig = figure;
t = tiledlayout(8, 7, ...
    'Padding', 'compact', ...   % minimal padding around edges
    'TileSpacing', 'compact');    % even smaller gaps between tiles
for j = 1:8
    for k = 1:7
        nexttile
        curImg = squeeze(curM_cropped(j,k,:,:));
        rootImg = curImg.^ 0.5;     % apply to each image
        imagesc(rootImg, clim);
         axis image;                   % keep correct aspect ratio, don't hide axes
        set(gca, 'XTick', [], 'YTick', []);  
        box on;
    end
end
colormap(fig, cmap_wr);
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Ticks = [];
cb.TickLabels = {};
set(gcf, 'Position', [100 100 900 1000]);

exportgraphics(fig, [mainFileName, '_LDSP_fullSF_map.png'])
exportgraphics(fig, [mainFileName, '_LDSP_fullSF_map.pdf'])








