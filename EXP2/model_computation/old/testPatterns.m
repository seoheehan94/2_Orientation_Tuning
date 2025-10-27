clear all;

addpath(imgFolder);
addpath('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli');

numOrientations = 8;
bandwidth = 1;
dims = backgroundSize;
numLevels = 1;

imgFolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/testPattern/';
imFilePath = dir(strcat(imgFolder, '*.svg'));
curFile = imFilePath(1).name;
[dirPath,mainFileName,~] = fileparts(curFile);


%% get vecLD
vecLD = importSVG(curFile);
vecLD = computeContourProperties(vecLD);
for c = 1:length(vecLD.numContours)
    vecLD.orientations{c} = 360-(vecLD.orientations{c});
end
vecLD = getContourPropertiesStats(vecLD);

[vecLD,MAT] = computeAllMATfromVecLD(vecLD);

%% nsdstim -vecLD
backgroundSize = [512 512];
renderSize = [357,357];
orifilename = ['vecLD' mainFileName '.mat'];


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

save(fullfile(imgFolder, orifilename),...
    'numOrientations','bandwidth','dims','modelOri','numLevels', 'oriMap');

%% nsdstim - control
interpImgSize = 714;
backgroundSize = 1024;
imgScaling = 0.5;

numOrientations = 8;
bandwidth = 1;
dims = [backgroundSize backgroundSize];
dims = dims*imgScaling;
numLevels = maxLevel(dims,bandwidth);
[freqRespsImag, freqRespsReal, pind] = makeQuadFRs(dims, numLevels, numOrientations, bandwidth);

pyramidfilename = ['LD' mainFileName '.mat'];

imgLD = renderLinedrawing(vecLD,[],[357,357]);
imgLD = rgb2gray((1-imgLD))*255;
bigImg = zeros(512);
offset = round((size(bigImg,1)-size(imgLD,1))/2);
bigImg(offset+[1:size(imgLD,1)],offset+[1:size(imgLD,2)]) = imgLD;

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

save(fullfile(imgFolder, pyramidfilename), 'interpImgSize','backgroundSize','imgScaling',...
    'numOrientations','bandwidth','dims','bigImg','sumOri','modelOri','numLevels');