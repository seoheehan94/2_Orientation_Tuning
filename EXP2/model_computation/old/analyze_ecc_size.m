
clear all;

filedir = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/surfaceData/';
condition = {'old', 'control','ori'};

% addpath(genpath('/usr/local/freesurfer/7.4.1/matlab'));
[~,M,mr] = load_mgh('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/surfaceData/sub1/oldBrain_sub1_lh_fsaverage.mgh');


%% bin eccentricity 

for isub=1:8
    %ecc
    ecc_lh = load_mgh([filedir, 'sub', num2str(isub), '/', 'prf_eccentricity_lh_fsaverage.mgh']);
    ecc_rh = load_mgh([filedir, 'sub', num2str(isub), '/', 'prf_eccentricity_rh_fsaverage.mgh']);

    Allecc_lh(:,isub) = ecc_lh;
    Allecc_rh(:,isub) = ecc_rh;
    

    edges = [4.2, 3, 2, 1, 0]; % Define edges in descending order
    values = [4, 3, 2, 1];     % Define corresponding values for the ranges
    % image size 8.4° of visual angle. Cut off should be half of the size
    ecc_lh(ecc_lh > 4.2) = 5;
    ecc_lh(ecc_lh == edges(1)) = 4;
    ecc_rh(ecc_rh > 4.2) = 5;
    ecc_rh(ecc_rh == edges(1)) = 4;

    for i = 1:length(values)
        ecc_lh(ecc_lh >= edges(i+1) & ecc_lh < edges(i)) = values(i);
        ecc_rh(ecc_rh >= edges(i+1) & ecc_rh < edges(i)) = values(i);
    end

    binnedEcc_lh(:,isub) = ecc_lh;
    binnedEcc_rh(:,isub) = ecc_rh;
end

    % %size
    % size_lh = load_mgh([filedir, 'sub', num2str(isub), '/', 'prf_size_lh_fsaverage.mgh']);
    % size_rh = load_mgh([filedir, 'sub', num2str(isub), '/', 'prf_size_rh_fsaverage.mgh']);
    % 
    % Allsize_lh(:,isub) = size_lh;
    % Allsize_rh(:,isub) = size_rh;
    % 
    % edges = [8.4, 7, 6, 5, 4, 3, 2, 1, 0];
    % values = [8, 7, 6, 5, 4, 3, 2, 1];
    % size_lh(size_lh > 8.4) = 9;
    % size_lh(size_lh == edges(1)) = 8;
    % size_rh(size_rh > 8.4) = 9;
    % size_rh(size_rh == edges(1)) = 8;
    % 
    % for i = 1:length(values)
    %     size_lh(size_lh >= edges(i+1) & size_lh < edges(i)) = values(i);
    %     size_rh(size_rh >= edges(i+1) & size_rh < edges(i)) = values(i);
    % end  
 

%% R2 by ecc
R2_ecc_mean = zeros(3, 5); % Rows: conditions, Columns: eccentricity bins

% Loop over conditions and eccentricity bins
for curcond = 1:3
    % Initialize temporary storage for R² values across subjects
    R2_ecc{curcond} = cell(1, 5); % 5 eccentricity bins

    for isub = 1:8
        % Load current subject's data
        ecc_lh = binnedEcc_lh(:, isub);
        ecc_rh = binnedEcc_rh(:, isub);
        R2_lh = load_mgh([filedir, 'sub', num2str(isub), '/', condition{curcond}, 'BrainR2_sub', num2str(isub), '_lh_fsaverage.mgh']);
        R2_rh = load_mgh([filedir, 'sub', num2str(isub), '/', condition{curcond}, 'BrainR2_sub', num2str(isub), '_rh_fsaverage.mgh']);

        % Loop over eccentricity bins
        for i = 1:5
            % Collect R² values for left and right hemisphere
            R2_ecc{curcond}{i} = [R2_ecc{curcond}{i}; R2_lh(ecc_lh == i); R2_rh(ecc_rh == i)];
        end
    end

    % Compute mean R² values by eccentricity
    for i = 1:5
        R2_ecc_mean(curcond, i) = mean(R2_ecc{curcond}{i}, 'omitnan'); % Omit NaNs
    end
end
  

%% save mean ecc, size as masks
meanecc_lh = mean(Allecc_lh,2, "omitnan");
meanecc_rh = mean(Allecc_rh,2, "omitnan");
% meansize_lh = mean(Allsize_lh,2, "omitnan");
% meansize_rh = mean(Allsize_rh,2, "omitnan");

edges = [4.2, 3, 2, 1, 0]; % Define edges in descending order
values = [4, 3, 2, 1];     % Define corresponding values for the ranges
% image size 8.4° of visual angle. Cut off should be half of the size
meanecc_lh(meanecc_lh > 4.2) = 5;
meanecc_lh(meanecc_lh == edges(1)) = 4;
meanecc_rh(meanecc_rh > 4.2) = 5;
meanecc_rh(meanecc_rh == edges(1)) = 4;
for i = 1:length(values)
    meanecc_lh(meanecc_lh >= edges(i+1) & meanecc_lh < edges(i)) = values(i);
    meanecc_rh(meanecc_rh >= edges(i+1) & meanecc_rh < edges(i)) = values(i);
end


% edges = [8.4, 7, 6, 5, 4, 3, 2, 1, 0];
% values = [8, 7, 6, 5, 4, 3, 2, 1];
% meansize_lh(meansize_lh > 8.4) = 9;
% meansize_lh(meansize_lh == edges(1)) = 8;
% meansize_rh(meansize_rh > 8.4) = 9;
% meansize_rh(meansize_rh == edges(1)) = 8;
% 
%     for i = 1:length(values)
%         meansize_lh(meansize_lh >= edges(i+1) & meansize_lh < edges(i)) = values(i);
%         meansize_rh(meansize_rh >= edges(i+1) & meansize_rh < edges(i)) = values(i);
%     end


for i = 1:5
    % Eccentricity masks
    ecc_mask_lh = (meanecc_lh == i);
    ecc_mask_rh = (meanecc_rh == i);
    
    save_label(filedir, ecc_mask_lh, ['ecc_level', num2str(i), '_lh.label']);
    save_label(filedir, ecc_mask_rh, ['ecc_level', num2str(i), '_rh.label']);
end

% for i = 1:9
%     % Size masks
%     size_mask_lh = (meansize_lh == i);
%     size_mask_rh = (meansize_rh == i);
% 
%     % Save left hemisphere label
%     save_label(filedir, size_mask_lh, ['size_level', num2str(i), '_lh.label']);
%     % Save right hemisphere label
%     save_label(filedir, size_mask_rh, ['size_level', num2str(i), '_rh.label']);
% end

function save_label(output_dir, mask, filename)
    % Convert logical mask to vertex indices
    vertex_indices = find(mask) - 1; % FreeSurfer uses 0-based indexing
    
    % Write to .label file
    fileID = fopen([output_dir, filename], 'w');
    fprintf(fileID, '#!ascii label file\n');
    fprintf(fileID, '%d\n', numel(vertex_indices)); % Number of vertices
    for v = vertex_indices'
        fprintf(fileID, '%d 0.0 0.0 0.0 1.0\n', v); % Vertex index and dummy coordinates
    end
    fclose(fileID);
end



%% plot

R2diff_avg_ecc_allH = mean(R2diff_avg_ecc, 1);
figure;
plot(R2diff_avg_ecc_allH(1,1:4),'-o', 'Color','black','LineWidth',2);
% ylim([0 0.3])
xticks(1:4);
xticklabels({'1', '2', '3', '4'});
xlabel('Eccentricity', 'FontSize', 15, 'FontName', 'Helvetica');
% legend('left','right','FontSize',15,'FontName','Helvetica');
legend('hide')
ax = gca;
ax.YAxis.FontSize = 15;
ax.YAxis.FontName = 'Helvetica';
ax.XAxis.FontSize = 15; 
ax.XAxis.FontName = 'Helvetica';
box off;
legend boxoff
% saveas(gcf,['r2diff_ecc' '.pdf']);


figure;
plot(R2_ecc_mean(1,1:4),'-o', 'Color','#0072BD','LineWidth',2);
hold on;
plot(R2_ecc_mean(2,1:4),'-o', 'Color','#65B74A','LineWidth',2);
hold on;
plot(R2_ecc_mean(3,1:4),'-o', 'Color','#F35872','LineWidth',2);

% ylim([0 0.3])
xticks(1:4);
xticklabels({'≤1', '≤2', '≤3', '≤4.2'});
xlabel('Eccentricity', 'FontSize', 15, 'FontName', 'Helvetica');
% legend('contour','photo-steerable pyramid','FontSize',15,'FontName','Helvetica');
legend('hide')
ax = gca;
ax.YAxis.FontSize = 15;
ax.YAxis.FontName = 'Helvetica';
ax.XAxis.FontSize = 15; 
ax.XAxis.FontName = 'Helvetica';
box off;
legend boxoff

saveas(gcf,['r2_ecc' '.pdf']);