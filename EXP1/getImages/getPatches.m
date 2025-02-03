clear all;
rng(4228);

%% Load data
original_folder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/experiment/getImages';
modelOriPhoto_folder = '/bwdata/NSDData/Seohee/stimuli/pyramid'; % Original image file folder
modelOriVecLD_folder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli/orientationfilter'; % Original image file folder
img_folder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/experiment/rawStimuli';
save_folder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/experiment/stimuli/candidates';
pracImg_folder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/experiment/rawPracticeStimuli';
pracSave_folder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/experiment/practiceStimuli';

cd(original_folder);
load("top550.mat")
imgNum=top550.imgNum;
imgFiles = dir(fullfile(img_folder, '*.png')); 

%img size
topLeft = [77, 77];
imgWidth =360;
imgHeight = 360;

% Create a circular mask
newRows = imgWidth/3;
newCols = imgHeight/3;
centerX = newCols / 2;
centerY = newRows / 2;
radius = newRows/2; % Set the radius of the circle
    
% Create a binary mask for the circle
[xGrid, yGrid] = meshgrid(1:imgHeight, 1:imgWidth);
circularMask = ((xGrid - centerX).^2 + (yGrid - centerY).^2) <= radius^2;
% circularMask = double(circularMask);
% circularMask(circularMask ==0) = NaN;

numAngles=8;
theta = linspace(0,2*pi,numAngles+1);%for circular calculation
theta = theta(1:end-1);

results = cell(length(imgNum), 8);
%%
for imgIdx = 1:550
    fprintf('%d. ...\n', imgIdx);
% load modelOri for photo
    filename = ['pyrImg' num2str(imgNum(imgIdx)) '.mat'];
    % filename = ['pyrImg' num2str(imgIdx) '.mat'];
    load(fullfile(modelOriPhoto_folder, filename), 'modelOri');
    curModelOri = cat(4, modelOri{:});
    curM = squeeze(mean(curModelOri,4));
    curM_cropped = zeros(8,imgWidth, imgHeight);
    for j = 1:8
        img = squeeze(curM(j, :, :));
        img_cropped = imcrop(img, [topLeft(1), topLeft(2), imgWidth-1, imgHeight-1]);
        curM_cropped(j, :, :) = img_cropped;
    end
    model.photo = curM_cropped;

% load modelOri for vecLD
    filename = ['oriImg' num2str(imgNum(imgIdx)) '.mat'];
    % filename = ['oriImg' num2str(imgIdx) '.mat'];
    load(fullfile(modelOriVecLD_folder, filename), 'oriMap');
    oriMap_cropped = imcrop(oriMap, [topLeft(1), topLeft(2), imgWidth-1, imgHeight-1]);
    reverseIdx = (oriMap_cropped < 0); % Identify values modified in original operation
    oriMap_cropped(reverseIdx) = oriMap_cropped(reverseIdx) + 180; % Add 180 back to reverse the previous operation
    model.vecLD = oriMap_cropped;
    % model.vecLD(model.vecLD ==0) = NaN;

    % get difference
    % diff = model.photo - model.vecLD;

    % Find the maximum diff location of the mask (within all possible positions)
    maxDiff = 0;
    maxPos = [];
    minDiff = 1;
    minPos = [];

    % Define the range to slide the circular mask
    for xShift = 1:(imgWidth - 2*radius)
        for yShift = 1:(imgHeight - 2*radius)
            % Create the shifted mask
            shiftedMask = circshift(circularMask, [xShift, yShift]);
            shiftedMask = double(shiftedMask); % Convert to double
            shiftedMask(shiftedMask == 0) = NaN;
            
            % Apply the shifted mask to the absolute image
            maskedSection_photo = zeros(size(model.photo)); 
            for i = 1:size(model.photo, 1)  
                maskedSection_photo(i, :, :) = squeeze(model.photo(i, :, :)) .* shiftedMask;
            end
            maskedSection_vecLD = model.vecLD .* shiftedMask;

            % Calculate the mean of the masked section
            allValues_photo = squeeze(mean(maskedSection_photo,[2 3], "omitnan"));
            curMean_photo = circ_mean(theta',allValues_photo);
            curMean_photo = mod(curMean_photo,2*pi);%from [-pi, pi] to [0 2pi]
            curMean_photo = curMean_photo./2;%range 0 to pi.
            
            allValues_vecLD = maskedSection_vecLD(~isnan(maskedSection_vecLD)); 
            allValues_vecLD_rad = 2*deg2rad(allValues_vecLD);
            curMean_vecLD = circ_mean(allValues_vecLD_rad);
            curMean_vecLD = mod(curMean_vecLD,2*pi);%from [-pi, pi] to [0 2pi]
            curMean_vecLD = curMean_vecLD./2;%range 0 to pi.
            
            diff = curMean_photo-curMean_vecLD;
            diff_converted = mod(diff + pi/2, pi) - pi/2;
   
            % Check if this is the maximum sum
            if nnz(~isnan(maskedSection_vecLD)) >= 50
                if abs(diff_converted) > abs(maxDiff)
                    maxDiff = diff_converted;
                    maxPos = [xShift, yShift];  % Store the position of the maximum
                    maxMean_photo = curMean_photo;
                    maxMean_vecLD = curMean_vecLD;
                    % fprintf('%d.%d ...\n', xShift,yShift);
                end

                if abs(diff_converted) < abs(minDiff)
                    minDiff = diff_converted;
                    minPos = [xShift, yShift];  % Store the position of the maximum
                    minMean_photo = curMean_photo;
                    minMean_vecLD = curMean_vecLD;
                    % fprintf('%d.%d ...\n', xShift,yShift);
                end
            end
        end
    end

    maxMean_photo_deg=rad2deg(maxMean_photo);
    maxMean_vecLD_deg=rad2deg(maxMean_vecLD);
    maxDiff_deg = rad2deg(maxDiff);
    minMean_photo_deg=rad2deg(minMean_photo);
    minMean_vecLD_deg=rad2deg(minMean_vecLD);
    minDiff_deg = rad2deg(minDiff);


    maxMask = circshift(circularMask, maxPos);
    minMask = circshift(circularMask, minPos);
    % maxMask = double(shiftedMask); % Convert to double
    % maxMask(shiftedMask == 0) = NaN;
    % maskedSection_photo = zeros(size(model.photo));
    % for i = 1:size(model.photo, 1)
    %     maskedSection_photo(i, :, :) = squeeze(model.photo(i, :, :)) .* maxMask;
    % end
    % maskedSection_vecLD = model.vecLD .* maxMask;

            
    % maxMask(isnan(maxMask)) = 0;        
    maxMask_resize = imresize(maxMask, [425, 425]);
    [row_max, col_max] = find(maxMask_resize);
    minRow_max = min(row_max);
    maxRow_max = max(row_max);
    minCol_max = min(col_max);
    maxCol_max = max(col_max);
    minMask_resize = imresize(minMask, [425, 425]);
    [row_min, col_min] = find(minMask_resize);
    minRow_min = min(row_min);
    maxRow_min = max(row_min);
    minCol_min = min(col_min);
    maxCol_min = max(col_min);

    curi = imread(fullfile(img_folder, ['img',num2str(imgNum(imgIdx)),'.png']));

    % Convert the image to greyscale
    % colormap gray
    gsImg = mean(curi, 3);

    % Create an RGB version of the grayscale image (3 identical channels)
    rgbImg = repmat(uint8(gsImg), [1 1 3]);

    % Create the alpha channel for transparency
    alphaChannel_max = uint8(255 * maxMask_resize); % Alpha channel (255 inside circle, 0 outside)
    croppedAlphaChannel_max = alphaChannel_max(minRow_max:maxRow_max, minCol_max:maxCol_max, :);
    alphaChannel_min = uint8(255 * minMask_resize);
    croppedAlphaChannel_min = alphaChannel_min(minRow_min:maxRow_min, minCol_min:maxCol_min, :);

    for channel = 1:3
        currentChannel_max = rgbImg(:,:,channel);  % Extract the channel
        currentChannel_max(~maxMask_resize) = 255;  % Set background to white
        maskedImg_max(:,:,channel) = currentChannel_max; % Store it back in rgbImg
        currentChannel_min = rgbImg(:,:,channel);  % Extract the channel
        currentChannel_min(~minMask_resize) = 255;  % Set background to white
        maskedImg_min(:,:,channel) = currentChannel_min; % Store it back in rgbImg
    end

    croppedMaskedImg_max = maskedImg_max(minRow_max:maxRow_max, minCol_max:maxCol_max, :);
    croppedMaskedImg_min = maskedImg_min(minRow_min:maxRow_min, minCol_min:maxCol_min, :);

    outputFile_max = fullfile(save_folder, ['img',num2str(imgNum(imgIdx)),'_max.png']);
    imwrite(croppedMaskedImg_max, outputFile_max, 'Alpha', croppedAlphaChannel_max);
    outputFile_min = fullfile(save_folder, ['img',num2str(imgNum(imgIdx)),'_min.png']);
    imwrite(croppedMaskedImg_min, outputFile_min, 'Alpha', croppedAlphaChannel_min);

    % Store results in the cell array
    results{imgIdx, 1} = sprintf('img%d.png', imgNum(imgIdx));  % imgIdx.name
    results{imgIdx, 2} = maxMean_photo;                       
    results{imgIdx, 3} = maxMean_photo_deg;                   
    results{imgIdx, 4} = maxMean_vecLD;                       
    results{imgIdx, 5} = maxMean_vecLD_deg;                  
    results{imgIdx, 6} = maxPos;                                
    results{imgIdx, 7} = maxDiff;                               
    results{imgIdx, 8} = rad2deg(maxDiff);                      
    results{imgIdx, 9} = minMean_photo;                      
    results{imgIdx, 10} = minMean_photo_deg;                   
    results{imgIdx, 11} = minMean_vecLD;                       
    results{imgIdx, 12} = minMean_vecLD_deg;                   
    results{imgIdx, 13} = minPos;                                
    results{imgIdx, 14} = minDiff;                              
    results{imgIdx, 15} = rad2deg(minDiff);                     

    % save(fullfile(original_folder, 'cropImagesInfo.mat'), 'results');

end


% Save the table to a file
resultsTable = cell2table(results, 'VariableNames', {'imgIdx_name', ...
    'maxMean_photo', 'maxMean_photo_deg', 'maxMean_vecLD', 'maxMean_vecLD_deg', 'maxPos', 'maxDiff', 'maxDiff_deg',...
    'minMean_photo', 'minMean_photo_deg','minMean_vecLD', 'minMean_vecLD_deg', 'minPos', 'minDiff', 'minDiff_deg'});
% writetable(resultsTable, fullfile(original_folder, 'results_table2.csv'));

imgnumColumn = top550.imgNum;
resultsTable = [table(imgnumColumn) resultsTable];
resultsTable.Properties.VariableNames{1} = 'imgnum';
save(fullfile(original_folder, 'results_tableTotal.mat'),'resultsTable');

% figure;imshow(finalMask)
% figure;imshow(maskedImg)
 % figure;imshow(croppedMaskedImg_min)
% 
% figure;imshow(curi)
% 
% figure;imagesc(maskedSection_vecLD)
% figure;imagesc(maskedSection_photo)
% 
% figure;imagesc(model.vecLD)
% % figure;imagesc(model.photo)
% photomodel=model.photo;
% photomodel(isnan(photomodel)) = 0;
% figure;imagesc(photomodel)
% figure;imagesc(squeeze(curM(8,:,:)))





% cand_folder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/experiment/stimuli/';
% save_folder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/experiment/practiceStimuli';
% 
% for k = 460:515
%     imgName = resultsTable.imgnum(k);
%     curImg = dir(fullfile([cand_folder, 'img',num2str(imgName), '*']));
%      % Loop through all matching files
%     for i = 1:length(curImg)
%         srcFile = fullfile(cand_folder, curImg(i).name);  % Full path of the source file
% 
%         % Check if the file exists
%         if exist(srcFile, 'file') == 2
%             % Move the file to the candidate folder
%             % movefile(srcFile, save_folder);
%             delete(srcFile);
%         else
%             disp(['File does not exist: ', curImg(i).name]);
%         end
%     end
% 
% end
% 
