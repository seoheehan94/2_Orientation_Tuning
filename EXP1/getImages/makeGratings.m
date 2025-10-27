makeGratings.m
%
% purpose: generate grating images used for experiment stimuli
% 




interpDeg = 8.4;
width = 142; % Image width and height (square grid)
pixPerDeg = width / interpDeg;
totalDeg = width / pixPerDeg; % Including gray background

numFreqs = 4;
minCycles = 1;
maxCycles = width / 4; % Reduce maximum cycles to avoid extreme frequencies
minCpd = minCycles / totalDeg;
maxCpd = maxCycles / totalDeg;
cpds = logspace(log10(minCpd), log10(maxCpd), numFreqs); % Spatial frequencies
cyclesPerImage = cpds * totalDeg;

angle = 0; % Orientation angle (horizontal)

% Phase
phase = 0;
phase = pi * phase / 180;

% Get grid of x and y coordinates
x = -(width - 1) / 2 : (width - 1) / 2;
y = -(width - 1) / 2 : (width - 1) / 2;
[xMesh, yMesh] = meshgrid(x, y);

% Convert spatial frequencies to cycles per pixel
freqs = cpds / pixPerDeg;

% Brightness levels (4 levels between 0.25 and 1)
brightnessLevels = linspace(0.25, 1, 4);

% Create a circular mask
radius = (width - 1) / 2; % Circle radius
circleMask = (xMesh.^2 + yMesh.^2) <= radius^2; % Logical circular mask

% Create a directory to save images if it doesn't exist
outputDir = 'gratingStimuli';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

for isf = 1 : length(freqs)
    sf = freqs(isf); % Spatial frequency
    
    % Calculate orientation (in radians)
    angleRad = pi * angle / 180;
    a = cos(angleRad) * sf * 2 * pi;
    b = sin(angleRad) * sf * 2 * pi;
    
    % Compute base grating pattern (cosine grating)
    grating = cos(a * xMesh + b * yMesh + phase);
    
    % Apply circular mask and save as PNG for each brightness level
    for ib = 1 : length(brightnessLevels)
        brightness = brightnessLevels(ib);
        
        % Create an RGBA image
        rgbaImage = zeros(width, width, 3); % Initialize an RGBA image
        adjustedGrating = (grating + 1) / 2 * brightness; % Scale between 0 and brightness level
        
        % Set the RGB channels
        rgbaImage(:, :, 1) = adjustedGrating; % Red channel
        rgbaImage(:, :, 2) = adjustedGrating; % Green channel
        rgbaImage(:, :, 3) = adjustedGrating; % Blue channel
        
        % Set the alpha channel (transparency)
        alphaChannel = zeros(width, width); % Initialize alpha channel
        alphaChannel(circleMask) = 1; % Set inside the circle to opaque
        % rgbaImage(:, :, 4) = alphaChannel; % Assign alpha channel
        
        % Save the image as PNG
        filename = sprintf('%s/grating_freq_%.2f_brightness_%.2f.png', outputDir, cpds(isf), brightness);
        imwrite(rgbaImage, filename, 'Alpha', alphaChannel);
    end
end
