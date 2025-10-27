clear all;
rng(4228);

imgFolder = 'Stimuli Set/stimuli/';
imgFiles_max = dir(fullfile(imgFolder, '*max.png')); 
imgList_max = {imgFiles_max(:).name};
imgFiles_min = dir(fullfile(imgFolder, '*min.png'));
imgList_min = {imgFiles_min(:).name};
gratingFolder = 'gratingStimuli';
gratingFiles = dir(fullfile(gratingFolder, '*.png')); 
gratingList = {gratingFiles(:).name};

% Ensure the output folder exists
outputFolder = 'participant_conditions';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Experiment setup
numParticipants = 20;  % Total number of participants
numAngle = 8;
angleRepeat = 7;
trialPerBlk=numAngle*angleRepeat;
imgPerBlk = trialPerBlk / 2;  % condition: max/min
grating = 16;
nblock = 8;
numTrials = trialPerBlk*nblock+16;       % Total trials per participant (each participant sees 200 images)


% Generate the rotation angles
angles = linspace(0, 315, numAngle); % 8 angles from 0 to 360 (0, 45, 90, ..., 315)
all_angles = repmat(angles, 1, angleRepeat); 


% Assign images and conditions to each participant
for p = 71:70+numParticipants
    p
    
    % Shuffle arrays
    shuffledGrating =gratingList(randperm(length(gratingList)));
    shuffledImg_max =imgList_max(randperm(length(imgList_max)));
    shuffledImg_min =imgList_min(randperm(length(imgList_min)));
    shuffledImg = [];
    for i = 1:nblock
        max_block = shuffledImg_max((i-1)*imgPerBlk + 1:i*imgPerBlk);
        min_block = shuffledImg_min((i-1)*imgPerBlk + 1:i*imgPerBlk);
        grating_block = shuffledGrating(2*i-1:2*i);
        img_block = [max_block min_block grating_block];
        img_block =img_block(randperm(length(img_block)));
        shuffledImg = [shuffledImg, img_block];
    end

    shuffled_angles = [];
    for i = 1:nblock
        angle_block = all_angles(randperm(length(all_angles)));
        angle_block = [angle_block angle_block(1:2)]; % Add two random angles for gratings
        shuffled_angles = [shuffled_angles, angle_block];
    end

   % Create a cell array to store both image names and conditions
   participantTrials = struct('imageName', [], 'rotationAngle', []);
   for trial = 1:numTrials
       participantTrials(trial).imageName = shuffledImg{trial};  % Image name
       participantTrials(trial).rotationAngle = shuffled_angles(trial);  % Angle rotation
   end

   participantTrials = participantTrials';

    % Save the condition to a .mat file
    % save(fullfile(outputFolder, sprintf('sub%02d_conditions.mat', p)), 'participantTrials');

    % Save the condition to a .JSON file
    Trials.trials = participantTrials;
    saveJSONfile(Trials, fullfile(outputFolder, sprintf('sona_sub%02d_conditions.json', p)));


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
