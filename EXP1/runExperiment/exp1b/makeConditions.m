clear all;
rng(4228);

imgFolder = '/Users/seoheehan/Library/CloudStorage/GoogleDrive-osooni94@gmail.com/My Drive/School/phD/Projects/Perceptual grouping/Orientation/experiment/Stimuli Set/stimuli';
imgFiles_max = dir(fullfile(imgFolder, '*max.png')); 
imgList_max = {imgFiles_max(:).name};
imgFiles_min = dir(fullfile(imgFolder, '*min.png'));
imgList_min = {imgFiles_min(:).name};
gratingFolder = '/Users/seoheehan/Library/CloudStorage/GoogleDrive-osooni94@gmail.com/My Drive/School/phD/Projects/Perceptual grouping/Orientation/experiment/Stimuli Set/gratingStimuli';
gratingFiles = dir(fullfile(gratingFolder, '*.png')); 
gratingList = {gratingFiles(:).name};

% Ensure the output folder exists
outputFolder = 'participant_conditions';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Experiment setup
numParticipants = 10;  % Total number of participants
numAngle = 8;
angleRepeat = 5;
imagePosition = 2; % left or right
responseCondition = 2; % bar or grating
trialPerBlk = numAngle * angleRepeat * imagePosition;
imgPerBlk = trialPerBlk / 2;  % condition: max/min
attentionCheckPerBlk = 2; % 2 attention check trials per block
nblock = 8; % Total blocks (5 bar + 5 grating)
numTrials = trialPerBlk*nblock+16;       % Total trials per participant 


% Generate the rotation angles
angles = linspace(0, 315, numAngle); % 8 angles from 0 to 360 (0, 45, 90, ..., 315)
all_angles = repmat(angles, 1, angleRepeat); 

% Block designs
evenDesign = {'bar', 'grating', 'grating', 'bar', 'bar', 'grating', 'grating', 'bar'}; % ABBAABBA
oddDesign = {'grating', 'bar', 'bar', 'grating', 'grating', 'bar', 'bar', 'grating'}; % BAABBAAB

% Assign images and conditions to each participant
for p = 41:numParticipants+40
    fprintf('Generating conditions for participant %d...\n', p);

    % Determine block design based on participant number
    if mod(p, 2) == 0
        blockDesign = evenDesign; % Even participant number
    else
        blockDesign = oddDesign; % Odd participant number
    end

    % Shuffle arrays
    shuffledAttentionCheck =gratingList(randperm(length(gratingList)));
    shuffledImg_max =imgList_max(randperm(length(imgList_max)));
    shuffledImg_min =imgList_min(randperm(length(imgList_min)));

    % Initialize trial data
    Trials.trials = [];

    % Generate trials for each block
    for block = 1:nblock

        participantTrials = struct('block',{}, 'imageName', {}, 'rotationAngle', {}, 'imagePosition', {}, 'imageCondition', {}, 'responseType', {});
        % Determine the condition for this block
        responseType = blockDesign{block};

        % Shuffle angles and image positions for this block
        blockAngles = [all_angles(randperm(length(all_angles))),all_angles(randperm(length(all_angles)))];
        imagePositions = repmat({'left', 'right'}, 1, numAngle * angleRepeat / 2);
        imagePositions = [imagePositions(randperm(length(imagePositions))),imagePositions(randperm(length(imagePositions)))];
        attentionCheckPositions = {'left', 'right'};

         % Assign images for this block
        blockImages = [shuffledImg_max((block-1)*imgPerBlk + 1:block*imgPerBlk), ...
                       shuffledImg_min((block-1)*imgPerBlk + 1:block*imgPerBlk)];

        % Add 2 attention check trials to each block
        attentionCheckTrials = shuffledAttentionCheck((block-1)*attentionCheckPerBlk + 1:block*attentionCheckPerBlk);

         % Create trials for this block
        for trial = 1:trialPerBlk
            participantTrials(end+1).block = block;
            participantTrials(end).imageName = blockImages{trial};
            participantTrials(end).rotationAngle = blockAngles(trial);
            if floor((trial-1)/(trialPerBlk/2)) == 0
                participantTrials(end).imageCondition = 'max';
            elseif floor((trial-1)/(trialPerBlk/2)) == 1
                participantTrials(end).imageCondition = 'min';
            end
            participantTrials(end).imagePosition = imagePositions{trial};
            participantTrials(end).responseType = responseType; % bar or grating
        end

        % Add attention check trials
        for trial = 1:attentionCheckPerBlk
            participantTrials(end+1).block = block;
            participantTrials(end).imageName = attentionCheckTrials{trial};
            participantTrials(end).rotationAngle = blockAngles(end - attentionCheckPerBlk + trial); % Use remaining angles
            participantTrials(end).imageCondition = 'attentionCheck'; % Attention check trial
            participantTrials(end).imagePosition = attentionCheckPositions{trial};
            participantTrials(end).responseType = responseType; % bar or grating
        end

        % Shuffle participantTrials
        shuffledIndices = randperm(length(participantTrials)); % Generate random indices
        participantTrials = participantTrials(shuffledIndices); % Reorder participantTrials

        Trials.trials = [Trials.trials, participantTrials];
    end

    
% Save the condition to a .JSON file
    saveJSONfile(Trials, fullfile(outputFolder, sprintf('sub%02d_conditions.json', p)));
end







%% save practice stimuli list
% imgFolder = 'practiceStimuli';
% imgFiles = dir(fullfile(imgFolder, '*.png')); % Change to your image format
% 
% practiceImages.practiceImages = cell2struct({imgFiles(:).name}, 'imageName');
% saveJSONfile(practiceImages, 'practiceImages.json');
% 
% 
% % %
% cand_folder = '/Users/seoheehan/orijudgeexp/stimuli/';
% save_folder = '/Users/seoheehan/Library/CloudStorage/GoogleDrive-osooni94@gmail.com/My Drive/School/phD/Projects/Perceptual grouping/Orientation/experiment/rawStimuli/candidates';
% 
% for k = 461:515
%     imgName = resultsTable.imgnum(k);
%     curImg = dir(fullfile([cand_folder, 'img',num2str(imgName), '*']));
%      % Loop through all matching files
%     for i = 1:length(curImg)
%         srcFile = fullfile(cand_folder, curImg(i).name);  % Full path of the source file
% 
%         % Check if the file exists
%         if exist(srcFile, 'file') == 2
%             % Move the file to the candidate folder
%              %movefile(srcFile, save_folder);
%             delete(srcFile);
%         else
%             disp(['File does not exist: ', curImg(i).name]);
%         end
%     end
% 
% end
% 
% 
