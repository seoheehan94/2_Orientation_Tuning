function labelData = read_label(filePath)
    % Open the file
    fid = fopen(filePath, 'r');
    if fid == -1
        error('Cannot open label file: %s', filePath);
    end
    
    % Skip the first line (comment)
    fgetl(fid);
    
    % Read the number of vertices
    numVertices = fscanf(fid, '%d', 1);
    
    % Read the rest of the data
    data = fscanf(fid, '%d %f %f %f %f', [5, numVertices])';
    
    % Close the file
    fclose(fid);
    
    % Parse the data
    labelData.vertexIndex = data(:, 1) + 1; % Convert to MATLAB's 1-based indexing
    labelData.coordinates = data(:, 2:4);   % x, y, z coordinates
    labelData.values = data(:, 5);         % Measurement values
end