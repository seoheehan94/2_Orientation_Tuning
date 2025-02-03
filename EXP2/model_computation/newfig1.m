% nsdStim.m
%
% associated with the following publication: Roth, ZN, Kay, K, and Merriam, EP (2022).
% Massive natural scene sampling reveals reliable coarse-scale orientation tuning in human V1
% DOI:
%
%   usage: nsdStim()
%   by: zvi roth
%   date: 7/29/2022
%   purpose: run natural scene stimuli through the models, get filter energy responses,
%   and save output
%   creates files used by: prfSampleModel.m

% uses the steerable pyramid: https://github.com/elimerriam/stimulusVignetting

close all
clear all

pyramidfolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/figures/';%to save model outputs
pyramidfilename = 'ForestExp.mat';
imgFile = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/figures/ForestExp.png';
origImg = imread(imgFile);
origImg = double(origImg);

%%
interpImgSize = 714;
backgroundSize = 1024;
imgScaling = 0.5;

numOrientations = 8;
bandwidth = 1;
dims = [backgroundSize backgroundSize];
dims = dims*imgScaling;
numLevels = maxLevel(dims,bandwidth);
[freqRespsImag, freqRespsReal, pind] = makeQuadFRs(dims, numLevels, numOrientations, bandwidth);

imgSizeX = size(origImg,1);
imgSizeY = size(origImg,2);

backgroundColor(1,1,:) = uint8([127,127,127]);
fixPoint(1,1,:) = [255 0 0];


%% photo - filter
[Xq, Yq] = meshgrid(linspace(1,imgSizeX, interpImgSize), linspace(1,imgSizeY, interpImgSize));
for irgb=1:3
    interpImg(:,:,irgb) = interp2(squeeze(origImg(:,:,irgb)), Xq, Yq);
end
%add red semi-transparent fixation point

interpImg(interpImgSize/2-8:interpImgSize/2+8,interpImgSize/2-8:interpImgSize/2+8,:) = ...
    (interpImg(interpImgSize/2-8:interpImgSize/2+8,interpImgSize/2-8:interpImgSize/2+8,:) + ...
    repmat(fixPoint,17,17,1))/2;

%%add background
bigImg = repmat(backgroundColor,backgroundSize,backgroundSize,1);
bigImg(1+backgroundSize/2-interpImgSize/2 : backgroundSize/2+interpImgSize/2, 1+backgroundSize/2-interpImgSize/2 : backgroundSize/2+interpImgSize/2,:) = interpImg(:,:,:);

%axis image
% colormap gray

%change to grayscale
% for now, simply by averaging across RGB channels
bigImg = mean(bigImg,3);

%DOWNSAMPLE
bigImg = imresize(bigImg,imgScaling);

[pyr, pind] = buildQuadBands(bigImg, freqRespsImag, freqRespsReal);
sumOri = cell(numLevels,1);
modelOri = cell(numLevels,1);
for ilev = 1:numLevels
    % loop over levels and orientations of the pyramid
    % initialize output
    sumOri{ilev}(:,:) = zeros(dims(1), dims(2));
    modelOri{ilev} = zeros(numOrientations, dims(1), dims(2));
    for orientation = 1:numOrientations
        thisBand = abs(accessSteerBand(pyr, pind, numOrientations,ilev, orientation)).^2;
        sumOri{ilev}(:,:) = sumOri{ilev}(:,:) + thisBand;
        modelOri{ilev}(orientation,:,:) = thisBand;
    end
end

figure;
oriList = [1, 3, 5, 7];
t = tiledlayout(4, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:12
    ori = ceil(i / 3); % Cycle through ori values 1 to 4
    sq = mod(i-1, 3) + 1;       % Cycle through sq values 1 to 3

    % Select the subplot position
    nexttile;
    imagesc(squeeze(modelOri{sq}(oriList(ori), :, :))); % Display the image
    colormap gray; % Use grayscale colormap
    axis off; % Turn off axes for a cleaner look
end


% Set figure properties
 set(gcf, 'Position', [100, 100, 750, 900]); % Set the figure size if desired

% Save the figure as a PNG file
saveas(gcf, '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/figures/energy_photofilter.png');

%% line drawing - filter
vecLD = traceLineDrawingFromRGB(imgFile);
imgLD = renderLinedrawing(vecLD,[],[357,357]);
imgLD = rgb2gray((1-imgLD))*255;
bigImg = zeros(512);
offset = round((size(bigImg,1)-size(imgLD,1))/2);
bigImg(offset+[1:size(imgLD,1)],offset+[1:size(imgLD,2)]) = imgLD;

% pass image through steerable pyramid
[pyr, pind] = buildQuadBands(bigImg, freqRespsImag, freqRespsReal);
sumOri = cell(numLevels,1);
modelOri = cell(numLevels,1);
for ilev = 1:numLevels
    % loop over levels and orientations of the pyramid
    % initialize output
    sumOri{ilev}(:,:) = zeros(dims(1), dims(2));
    modelOri{ilev} = zeros(numOrientations, dims(1), dims(2));
    for orientation = 1:numOrientations
        thisBand = abs(accessSteerBand(pyr, pind, numOrientations,ilev, orientation)).^2;
        sumOri{ilev}(:,:) = sumOri{ilev}(:,:) + thisBand;
        modelOri{ilev}(orientation,:,:) = thisBand;
    end
end
            
figure;
oriList = [1, 3, 5, 7];
t = tiledlayout(4, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:12
    ori = ceil(i / 3); % Cycle through ori values 1 to 4
    sq = mod(i-1, 3) + 1;       % Cycle through sq values 1 to 3

    % Select the subplot position
    nexttile;
    imagesc(squeeze(modelOri{sq}(oriList(ori), :, :))); % Display the image
    colormap gray; % Use grayscale colormap
    axis off; % Turn off axes for a cleaner look
end


% Set figure properties
 set(gcf, 'Position', [100, 100, 750, 900]); % Set the figure size if desired

% Save the figure as a PNG file
saveas(gcf, '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/figures/energy_LDfilter.png');


%% contour
backgroundSize = [512 512];
renderSize = [357,357];
numOrientations = 8;
bandwidth = 1;
dims = backgroundSize;
numLevels = 1;
[freqRespsImag, freqRespsReal, pind] = makeQuadFRs(dims, numLevels, numOrientations, bandwidth);

vecLD = computeContourProperties(vecLD);
vecLD = getContourPropertiesStats(vecLD);
drawLinedrawing(vecLD)

for c = 1:vecLD.numContours
    vecLD.orientations{c} = vecLD.orientations{c}-90;
end

oriMap = generateOrientationMap(vecLD, NaN, backgroundSize, renderSize);

binWidth2 = 90 / length(vecLD.orientationBins);
horIdx = (oriMap > (180-binWidth2));
oriMap(horIdx) = oriMap(horIdx) - 180;

%vecLD.orientationBins
for binIdx = 1: length(vecLD.orientationBins)
    oribinmap{1,binIdx} = (abs(oriMap-vecLD.orientationBins(binIdx)) <= binWidth2);

end
modelOri=cat(3,oribinmap{:});
modelOri = permute(modelOri,[3 1 2]);

figure;
oriList = [1, 7, 5, 3];
t = tiledlayout(4, 1,'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:4
    % Select the subplot position
    nexttile;
    imagesc(squeeze(modelOri(oriList(i), :, :))); % Display the image
    colormap gray; % Use grayscale colormap
    axis off; % Turn off axes for a cleaner look
end


% Set figure properties
 set(gcf, 'Position', [100, 100, 250, 900]); % Set the figure size if desired

% Save the figure as a PNG file
saveas(gcf, '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/figures/energy_vecLD.png');
