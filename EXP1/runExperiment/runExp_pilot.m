% Name          : Orientation experiment
% Author        : Seohee Han
% Environment   :
% Windows 11 64bit, Matlab R2024b 64bit, PTB 3.0.19.14
% Version       : 2024.09
% Note          : Orientation judgment with 30/150ms

%% ---------- Initialize ----------%
clear all; % clear all pre-existing variables (interference prevention)
try % error handling start (error -> automatically close all screens)
    rng('shuffle');
    echo off; % do not re-show my command
    clc; % clear command window
    fclose('all'); % close all pre-opened files (interference prevention)
    GetSecs; WaitSecs(0); % pre-load timing function for precision timing

    % Set up Psychtoolbox
    Screen('Preference', 'SkipSyncTests', 1); % Skip sync tests for demo purposes, 0 for experiment
    PsychDefaultSetup(2);

    expName     = 'oriJudge'; % experiment name for data files
    participantID = input('Enter participant ID: ', 's');
    participantID = sprintf('%02d', str2double(participantID));  % Format to 2-digit string
        
    % Check for saved state
    saveFileName = ['results/', expName, '_results_pilot' participantID '.mat'];
    if exist(saveFileName, 'file')
        load(saveFileName, 'data', 'lastTrialNumber', 'lastPhase', 'maskPath_prac', 'practiceAngles','practiceTimes');
        startTrial = lastTrialNumber + 1;
        startPhase = lastPhase;
        outputFile = fopen(['results/', expName, '_results_pilot' participantID '.csv'], 'a');
    else
        startTrial = 1;
        startPhase = 'practice'; % Start with practice phase if no previous state
        data = struct('participantID', {},'trialPhase', {}, 'trialNumber', {}, 'imageName', {}, 'maskName', {}, 'presentationTime', {},'imageDuration', {}, 'rotationAngle', {}, 'startOrientation', {},'responseOrientation', {}, 'barResponseTime', {}, 'confidenceLevel', {});
        outputFile = fopen(['results/', expName, '_results_pilot' participantID '.csv'], 'w');
        fprintf(outputFile, 'participantID,trial_phase,trial,image_name,maskName,presentationTime,imageDuration,rotationAngle,startOrientation,response_orientation,bar_response_time,confidenceLevel\n'); % Header

    end

    % Load participant condition file
    conditionFile = fullfile('participant_conditions', sprintf('pilot_%02d_conditions.mat', str2double(participantID)));
    if exist(conditionFile, 'file')
        load(conditionFile, 'participantTrials','randomOrderIdx'); % Load the conditions file
    else
        error('Condition file for participant %s not found!', participantID);
    end

    %% ---------- Configuration ----------%
    % Display setup
    ListenChar(2);
    screens = Screen('Screens');
    screenNumber = max(screens);
    black = BlackIndex(screenNumber);
    white = WhiteIndex(screenNumber);
    grey = white / 2; % Define grey color
    % [w, rect] = PsychImaging('OpenWindow', screenNumber, grey, [0,0,1200,800]);
    [w, rect] = PsychImaging('OpenWindow', screenNumber, grey);
    [screenXpixels, screenYpixels] = Screen('WindowSize', w);
    ifi = Screen('GetFlipInterval', w);
    topPriorityLevel = MaxPriority(w);
    Priority(topPriorityLevel);
    % Set blend function for transparency
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    % Notice
    notice1 = 'Welcome to the experiment'; % notice message for phase 1
    notice2 = ['In this experiment, you will view a series of images.\n' ...
        'Your task is to report the predominant orientation of each image.\n\n' ...
        'You will have a limited time to view each image,\n' ...
        'and then you will be asked to adjust a bar to match\n' ...
        'the mean orientation you perceived.\n' ...
        'There will be 3 blocks in total.\n\n' ...
        'Press ''spacebar'' to begin practice trials'];
    practiceEnd = 'Practice phase is done.';
    mainStart = ['Remember, your task is to report the predominant orientation of each image.\n' ...
        'Press ''spacebar'' to begin the experiment'];
    noticeEnd = 'Thank you for your participation! Please call the research assistant.';
    confidence = 'Rate your confidence (1=low, 5=high)\n 1 2 3 4 5';


    % Key Setup
    goKey       = KbName('space'); % response key to proceed
    RAKey       = KbName('z'); % RA key to proceed
    quitKey      = KbName('q'); % Quit key
    confidenceKeyList   = {'1!', '2@', '3#', '4$', '5%'}; % keys for surprise WM task
    for k = 1:numel(confidenceKeyList)
        confidenceKeys(k) = KbName(confidenceKeyList{k});
    end

    % Number of trials and blocks
    nTrials = length(participantTrials); % Should be 200, as per your setup
    nTrials_prac = 18;
    nBlocks = 3; % Number of blocks
    trialsPerBlock = nTrials / nBlocks; % Number of trials per block
    breakInterval = trialsPerBlock; % Number of trials before each break
    numAngle = 8;

    % image settings
    fixPx     = 12;
    stimSize = [425 425];

    % Image presentation times
    shortPresentation = 0.03; % 30 ms
    longPresentation = 0.15; % 150 ms
    unlimitPresentation = 0.2; % 200 ms
    inputDelay    = 1; % form input delay to avoid skipping (seconds)
    fixTime       = 0.5; % fixation duration
    maskTime = 0.4;
    timeout = 10; % Timeout after 10 seconds

    % Define the size and coordinates of the fixation cross
    fixCrossDimPix = 20; % Dimension of the fixation cross
    xCenter = rect(3)/2; % X-coordinate for screen center
    yCenter = rect(4)/2; % Y-coordinate for screen 10center
    xCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
    yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
    allCoords = [xCoords; yCoords];
    lineWidthPix = 3; % Thickness of the fixation cross lines

    % Preload images and convert them to textures
    imgFolder = 'stimuli';
    preloadedTextures = containers.Map; % Store image textures to avoid reloading
    for i = 1:length(participantTrials)
        imgName = participantTrials{i, 1}; % Get image name
        if ~isKey(preloadedTextures, imgName)
            imgPath = fullfile(imgFolder, imgName);
            [img, ~, alpha] = imread(imgPath); % Read image and alpha channel
            img(:, :, 4) = alpha;
            imgTexture = Screen('MakeTexture', w, img); % Create texture with RGBA image
            preloadedTextures(imgName) = imgTexture; % Store texture in map
        end
    end

    imgFolder_prac = 'practiceStimuli';
    imgPath_prac = dir(fullfile(imgFolder_prac, '*.png'));
    imgPath_prac = imgPath_prac(randperm(length(imgPath_prac)));
    preloadedTextures_prac = containers.Map; % Store image textures to avoid reloading
    for i = 1:nTrials_prac
        imgName = imgPath_prac(i).name; % Get image name
        if ~isKey(preloadedTextures_prac, imgName)
            imgPath = fullfile(imgFolder_prac, imgName);
            [img, ~, alpha] = imread(imgPath); % Read image and alpha channel
            img(:, :, 4) = alpha;
            imgTexture = Screen('MakeTexture', w, img); % Create texture with RGBA image
            preloadedTextures_prac(imgName) = imgTexture; % Store texture in map
        end
    end

    maskFolder = 'masks';
    maskPath = dir(fullfile(maskFolder, '*.png'));
    preloadedTextures_mask = containers.Map; % Store image textures to avoid reloading
    for i = 1:length(participantTrials)
        imgName = participantTrials{i, 4}; % Get image name
        if ~isKey(preloadedTextures_mask, imgName)
            imgPath = fullfile(maskFolder, imgName);
            [img, ~, alpha] = imread(imgPath); % Read image and alpha channel
            img(:, :, 4) = alpha;
            imgTexture = Screen('MakeTexture', w, img); % Create texture with RGBA image
            preloadedTextures_mask(imgName) = imgTexture; % Store texture in map
        end
    end

    % choose 18 random masks
     if ~exist('maskPath_prac', 'var')
        maskPath_prac = maskPath(randperm(length(maskPath)));
        maskPath_prac = maskPath_prac(1:nTrials_prac);
     end

    % choose 18 random rotation angle
     if ~exist('practiceAngles', 'var')
        angles = linspace(0, 315, numAngle); % 8 angles from 0 to 360 (0, 45, 90, ..., 315)
        practiceAngles = datasample(angles,nTrials_prac,'Replace',true);
     end

    % Randomize practice presentation times
     if ~exist('practiceTimes', 'var')
         condition30 = [repmat(0.03, 1, nTrials_prac/3)];
         condition150 = [repmat(0.15, 1, nTrials_prac/3)];
         conditionUnlimit = [repmat(10, 1, nTrials_prac/3)];
         conditions = {condition30, condition150, conditionUnlimit};
         shuffledConditions = conditions(randomOrderIdx);
         practiceTimes = [shuffledConditions{:}];
     end

      expTimer    = GetSecs;
    %% ---------- Notice phase ----------%
    % Display notice screens
    displayNotice(w, notice1, white);
    RestrictKeysForKbCheck(RAKey);
    FlushEvents('keyDown'); % flush pre-pressed key events
    while 1
        [keyIsDown, ~, ~] = KbCheck();
        if keyIsDown
            break
        end
    end

    if strcmp(startPhase, 'practice')
        displayNotice(w, notice2, white);
        RestrictKeysForKbCheck(goKey);
        FlushEvents('keyDown'); % Flush all previous key presses
        while 1
            [keyIsDown, ~, ~] = KbCheck();
            if keyIsDown
                break
            end
        end
    end

    %% ---------- Practice phase ----------%
    if strcmp(startPhase, 'practice')        
        % Continue from where it was interrupted
        startTrial = max(startTrial, 1); % Ensure trial index starts from 1

        for trial = startTrial:nTrials_prac
            % Get the practice image texture and presentation time
            imgName = imgPath_prac(trial).name;  % Practice image name
            presentationTime = practiceTimes(trial); % Get randomized presentation time
            rotationAngle = practiceAngles(trial);
            maskName = maskPath_prac(trial).name;

            % Load image texture from map
            imgTexture = preloadedTextures_prac(imgName);
            maskTexture = preloadedTextures_mask(maskName);

            % Calculate position to draw the stimulus centered on the screen
            dstRect = CenterRectOnPointd([0 0 stimSize], screenXpixels / 2, screenYpixels / 2);
            
            % Display the fixation
            Screen('DrawLines', w, allCoords, lineWidthPix, [0 0 0], [xCenter yCenter], 2);
            fixationFlipTime = Screen('Flip', w);
            WaitSecs(fixTime);

            if presentationTime ~= 10
                % Display the image
                Screen('DrawTexture', w, imgTexture, [], dstRect, rotationAngle);
                imageFlipTime = Screen('Flip', w);
                WaitSecs(presentationTime);

                % Display a mask
                Screen('DrawTexture', w, maskTexture, [], dstRect);
                maskFlipTime = Screen('Flip', w);
                WaitSecs(maskTime); % Mask duration

                imageDuration = maskFlipTime - imageFlipTime;

                % Collect the bar orientation using mouse
                [startOrientation, responseOrientation, barResponseTime] = collectBarResponse(w, screenXpixels, screenYpixels, grey, timeout, quitKey);

            elseif presentationTime == 10
                % Display the image
                Screen('DrawTexture', w, imgTexture, [], dstRect, rotationAngle);
                imageFlipTime = Screen('Flip', w);

                % Wait for 0.2 seconds with the image on the screen before collecting the response
                WaitSecs(unlimitPresentation);

                % Collect bar orientation response, with red bar drawn on top of the image
                [startOrientation, responseOrientation, barResponseTime] = collectBarResponseWithRedBar(w, imgTexture, dstRect, rotationAngle, screenXpixels, screenYpixels, grey, timeout - unlimitPresentation, quitKey);
                imageDuration = barResponseTime + unlimitPresentation;
            end

            % get confidence level
            displayNotice(w, confidence, white);
            RestrictKeysForKbCheck(confidenceKeys);
            FlushEvents('keyDown'); % Flush all previous key presses
            while 1
                [keyIsDown, ~, keyCode] = KbCheck();
                if keyIsDown
                    keyPressed = KbName(keyCode);
                    % Check if the key pressed is a number between 1 and 5
                    if ismember(keyPressed, {'1!', '2@', '3#', '4$', '5%'}) % Handles different keyboard layouts
                        confidenceLevel = str2double(keyPressed(1)); % Convert first character to a number
                    end
                    break
                end
            end
            
            % Save mat
            lastTrialNumber = trial;
            lastPhase = 'practice';
            trialData = struct('participantID', participantID, 'trialPhase', 'practice', 'trialNumber', trial, ...
                'imageName', imgName, 'maskName', maskName,'presentationTime', presentationTime, ...
                'imageDuration', imageDuration, 'rotationAngle', rotationAngle, ...
                'startOrientation', startOrientation, 'responseOrientation', responseOrientation, 'barResponseTime', barResponseTime, 'confidenceLevel', confidenceLevel);
            data = [data; trialData];
            save(saveFileName, 'data','lastTrialNumber','lastPhase', 'maskPath_prac', 'practiceAngles','practiceTimes'); % Save data to .mat file

            % Save csv
            fprintf(outputFile, '%s,practice,%d,%s,%s,%.2f,%.4f,%d,%.4f,%.4f,%.4f,%d\n', ...
                participantID, trial, imgName, maskName, presentationTime, imageDuration, rotationAngle, startOrientation, responseOrientation, barResponseTime, confidenceLevel);
        end

        lastPhase = 'main';
        lastTrialNumber = 0;
        save(saveFileName, 'data','lastTrialNumber','lastPhase', 'maskPath_prac', 'practiceAngles','practiceTimes'); % Save data to .mat file
        
        % End of practice phase
        displayNotice(w, practiceEnd, white);
        RestrictKeysForKbCheck(RAKey);
        FlushEvents('keyDown'); % Flush all previous key presses
        while 1
            [keyIsDown, ~, ~] = KbCheck();
            if keyIsDown
                break
            end
        end
    end

    displayNotice(w, mainStart, white);
    RestrictKeysForKbCheck(goKey);
    FlushEvents('keyDown'); % Flush all previous key presses
    while 1
        [keyIsDown, ~, ~] = KbCheck();
        if keyIsDown
            break
        end
    end

    %% ---------- Main Experiment ----------%
    startIndex = nTrials_prac; % Adjust index after practice trials

    for trial = startTrial:nTrials
        % Get image and presentation time for this trial
        imgName = participantTrials{trial, 1};  % Image name
        presentationTime = participantTrials{trial, 2};
        rotationAngle = participantTrials{trial, 3};
        maskName = participantTrials{trial, 4};

        % Load image texture from map
        imgTexture = preloadedTextures(imgName);
        maskTexture = preloadedTextures_mask(maskName);

        % Calculate position to draw the stimulus centered on the screen
        dstRect = CenterRectOnPointd([0 0 stimSize], screenXpixels / 2, screenYpixels / 2);

        % Display the fixation
        Screen('DrawLines', w, allCoords, lineWidthPix, [0 0 0], [xCenter yCenter], 2);
        fixationFlipTime = Screen('Flip', w);
        WaitSecs(fixTime);

        if presentationTime ~= 10
            % Display the image and record flip time
            Screen('DrawTexture', w, imgTexture, [], dstRect, rotationAngle);
            imageFlipTime = Screen('Flip', w);
            WaitSecs(presentationTime);

            % Display a mask and record mask flip time
            Screen('DrawTexture', w, maskTexture, [], dstRect);
            maskFlipTime = Screen('Flip', w);
            WaitSecs(maskTime); % Mask duration

            imageDuration = maskFlipTime - imageFlipTime;

            % Collect the bar orientation and measure response times
            [startOrientation, responseOrientation, barResponseTime] = collectBarResponse(w, screenXpixels, screenYpixels, grey, timeout, quitKey);

        elseif presentationTime == 10
            % Display the image
            Screen('DrawTexture', w, imgTexture, [], dstRect, rotationAngle);
            imageFlipTime = Screen('Flip', w);

            % Wait for 0.2 seconds with the image on the screen before collecting the response
            WaitSecs(unlimitPresentation);

            % Collect bar orientation response, with red bar drawn on top of the image
            [startOrientation, responseOrientation, barResponseTime] = collectBarResponseWithRedBar(w, imgTexture, dstRect, rotationAngle, screenXpixels, screenYpixels, grey, timeout - unlimitPresentation, quitKey);
            imageDuration = barResponseTime + unlimitPresentation;
        end

        % get confidence level
        displayNotice(w, confidence, white);
        RestrictKeysForKbCheck(confidenceKeys);
        FlushEvents('keyDown'); % Flush all previous key presses
        while 1
            [keyIsDown, ~, keyCode] = KbCheck();
            if keyIsDown
                keyPressed = KbName(keyCode);
                % Check if the key pressed is a number between 1 and 5
                if ismember(keyPressed, {'1!', '2@', '3#', '4$', '5%'}) % Handles different keyboard layouts
                    confidenceLevel = str2double(keyPressed(1)); % Convert first character to a number
                end
                break
            end
        end

        % Save mat
        lastTrialNumber = trial;
        lastPhase = 'main';
        trialData = struct('participantID', participantID, 'trialPhase', 'main', 'trialNumber', trial, ...
                'imageName', imgName, 'maskName', maskName,'presentationTime', presentationTime, ...
                'imageDuration', imageDuration, 'rotationAngle', rotationAngle, ...
                'startOrientation', startOrientation, 'responseOrientation', responseOrientation, 'barResponseTime', barResponseTime,'confidenceLevel', confidenceLevel);
        data = [data; trialData];
        save(saveFileName, 'data','lastTrialNumber','lastPhase', 'maskPath_prac', 'practiceAngles','practiceTimes'); % Save data to .mat file

        % Save csv
        fprintf(outputFile, '%s,main,%d,%s,%s,%.2f,%.4f,%d,%.4f,%.4f,%.4f,%d\n', ...
            participantID, trial, imgName, maskName, presentationTime, imageDuration, rotationAngle,  startOrientation,responseOrientation, barResponseTime,confidenceLevel);

        % Insert a break every 'breakInterval' trials
        if mod(trial, breakInterval) == 0 && trial < nTrials
            blocksLeft = (nTrials - trial) / trialsPerBlock;
            breakMessageWithBlocks = sprintf('You are now on a break. Please take a short rest.\nBlocks left: %d\nPress ''spacebar'' to continue', blocksLeft);
            displayNotice(w, breakMessageWithBlocks, white);
            RestrictKeysForKbCheck([goKey, quitKey]); % Restrict key press to continue key or quit key
            FlushEvents('keyDown'); % Flush all previous key presses
            while 1
                [keyIsDown, ~, keyCode] = KbCheck();
                if keyIsDown
                    % Check if quitKey was pressed
                    if keyCode(quitKey)
                        % Quit the experiment
                        sca; % Close all screens
                        fclose('all'); % Close all open files
                        error('Experiment terminated by user during the break.');
                    end
                    % Check if goKey (spacebar) was pressed to continue
                    if keyCode(goKey)
                        break; % Exit the loop and continue the experiment
                    end
                end
            end
        end
    end
    
    % Close the file
    EXPTIME = GetSecs - expTimer;
    fclose(outputFile);

    % Save the results
    save(saveFileName, 'data', 'EXPTIME');

    % End of experiment notice
    displayNotice(w, noticeEnd, white);
    RestrictKeysForKbCheck(RAKey);
    FlushEvents('keyDown'); % Flush all previous key presses
    while 1
        [keyIsDown, ~, ~] = KbCheck();
        if keyIsDown
            break
        end
    end
        
    % Close the screen
    ListenChar(0);
    sca;

    
catch ME % if an error occured, run following script and end the experiment
Screen('CloseAll');
fclose('all');
ShowCursor;
ListenChar(0);
rethrow(ME); % show me the error
end

%% ---------- Functions ----------%
% Function to display a notice screen
function displayNotice(window, message, textColor)
Screen('TextSize', window, 40); % Set text size
DrawFormattedText(window, message, 'center', 'center', textColor);
Screen('Flip', window);
end

% Function to collect bar orientation response from participant
function [startOrientation, orientation, barResponseTime] = collectBarResponse(window, screenXpixels, screenYpixels, grey, timeout, quitKey)
    % Randomize the initial bar orientation between 0 and 180 degrees
    startOrientation = rand() * 180;
    barOrientation = startOrientation;
    barLength = 200;
    barWidth = 10;
    sensitivity = 0.2; % Lower sensitivity factor to make movement smoother
    RestrictKeysForKbCheck([]);

    % Set up the mouse
    HideCursor;
    Screen('FillRect', window, grey);
    Screen('Flip', window);

    % Set the initial position of the mouse
    [startX, startY] = GetMouse(window);

    % Record the start time for bar movement
    responseStartTime = GetSecs;

    % Loop until a mouse click to confirm
    while 1 
        % Check if the timeout duration(10 sec) has passed
        currentTime = GetSecs;
        if currentTime - responseStartTime > timeout
            barResponseTime = NaN;
            orientation = NaN;
            break; % Exit the loop and move to the next trial
        end

        % Get the current mouse position
        [mouseX, mouseY, buttons] = GetMouse(window);

        % Calculate the change in mouse position
        deltaX = mouseX - startX;

        % Adjust the bar orientation based on mouse movement, scaled by sensitivity
        barOrientation = barOrientation - sensitivity * deltaX;  % Reverse to ensure clockwise is increasing angle
        barOrientation = mod(barOrientation, 180); % Keep orientation in range [0, 180)

        % Clear the screen and draw the bar
        Screen('FillRect', window, grey);
        drawRotatedBar(window, screenXpixels, screenYpixels, barOrientation, barLength, barWidth, [0,0,0]);
        Screen('Flip', window);

        % Check for a mouse button press to confirm the response
        if any(buttons)
            barResponseTime = GetSecs - responseStartTime; % Record the mouse click time
            break; % Exit the loop when the mouse button is clicked
        end

        % Check for the quit key press
        [keyIsDown, ~, keyCode] = KbCheck();
        if keyIsDown && keyCode(quitKey)
            % Force quit the experiment
            sca; % Close all screens
            fclose('all'); % Close all files
            error('Experiment terminated by user using quit key.');
        end

        % Update the starting position for the next frame
        startX = mouseX;

        WaitSecs(0.05); % Small pause to avoid multiple fast updates
    end

    % If no timeout, store the final bar orientation
    if ~isnan(barResponseTime)
        orientation = barOrientation;
    end
end

% Function to collect bar orientation response with red bar drawn on top of the image
function [startOrientation, orientation, barResponseTime] = collectBarResponseWithRedBar(window, imgTexture, dstRect, rotationAngle, screenXpixels, screenYpixels, grey, timeout, quitKey)
    % Randomize the initial bar orientation between 0 and 180 degrees
    startOrientation = rand() * 180;
    barOrientation = startOrientation;
    barLength = 200;
    barWidth = 10;
    sensitivity = 0.2; % Lower sensitivity factor to make movement smoother
    RestrictKeysForKbCheck([]);
    
    % Set up the mouse
    HideCursor;
    Screen('Flip', window);
    
    % Set the initial position of the mouse
    [startX, startY] = GetMouse(window);
    
    % Record the start time for bar movement
    responseStartTime = GetSecs;
    
    % Loop until a mouse click to confirm or timeout
    while 1
        % Check if the timeout duration has passed
        currentTime = GetSecs;
        if currentTime - responseStartTime > timeout
            barResponseTime = NaN;
            orientation = NaN;
            break; % Exit the loop and move to the next trial
        end
        
        % Get the current mouse position
        [mouseX, mouseY, buttons] = GetMouse(window);
        
        % Calculate the change in mouse position
        deltaX = mouseX - startX;
        
        % Adjust the bar orientation based on mouse movement, scaled by sensitivity
        barOrientation = barOrientation - sensitivity * deltaX;  % Reverse to ensure clockwise is increasing angle
        barOrientation = mod(barOrientation, 180); % Keep orientation in range [0, 180)
        
        % Clear the screen and draw the image with the red bar on top
        Screen('DrawTexture', window, imgTexture, [], dstRect, rotationAngle);
        drawRotatedBar(window, screenXpixels, screenYpixels, barOrientation, barLength, barWidth, [255 0 0]); % Red bar
        Screen('Flip', window);
        
        % Check for a mouse button press to confirm the response
        if any(buttons)
            barResponseTime = GetSecs - responseStartTime; % Record the mouse click time
            break; % Exit the loop when the mouse button is clicked
        end
        
        % Check for the quit key press
        [keyIsDown, ~, keyCode] = KbCheck();
        if keyIsDown && keyCode(quitKey)
            % Force quit the experiment
            sca; % Close all screens
            fclose('all'); % Close all files
            error('Experiment terminated by user using quit key.');
        end
        
        % Update the starting position for the next frame
        startX = mouseX;
        
        WaitSecs(0.05); % Small pause to avoid multiple fast updates
    end
    
    % If no timeout, store the final bar orientation
    if ~isnan(barResponseTime)
        orientation = barOrientation;
    end
end

% Function to draw a rotated bar for orientation adjustment
function drawRotatedBar(window, screenXpixels, screenYpixels, angle, length, width, color)
    % Calculate the center points of the screen
    xCenter = screenXpixels / 2;
    yCenter = screenYpixels / 2;

    % Define the two endpoints of the line (centered at the middle)
    x1 = 0;
    y1 = -length / 2;
    x2 = 0;
    y2 = length / 2;

    % Create the rotation matrix for the given angle
    rotateMat = [cosd(angle), -sind(angle); sind(angle), cosd(angle)];

    % Apply the rotation to the endpoints
    point1 = rotateMat * [x1; y1];
    point2 = rotateMat * [x2; y2];

    % Translate the points to be centered on the screen
    x1Rot = point1(1) + xCenter;
    y1Rot = point1(2) + yCenter;
    x2Rot = point2(1) + xCenter;
    y2Rot = point2(2) + yCenter;

    % Draw the rotated bar as a black line
    Screen('DrawLine', window, color, x1Rot, y1Rot, x2Rot, y2Rot, width);
end
