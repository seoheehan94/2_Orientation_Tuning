%%cd '/home/hanseohe/Documents/GitHub/nsdOtopy'
clear all;
cd '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/brainVolume';
corrList = [];
folderDir = '/bwdata/NSDData/Seohee/Orientation/';
condition = {'old', 'ori', 'control'};
subfolderName = {'prfsample/', 'prfsample_Ori/', 'prfsample_Ori_control/'};
saveName = {'oldBrain','oriBrain','controlBrain'};

%% voxel preference
for curcond = 1:3
    for isub = 1:8


        %% get volume data
        fileName = [folderDir, subfolderName{curcond}, 'voxModelPref_sub', num2str(isub), '.mat'];
        load(fileName);

        ourBrain = visRoiData;
        ourBrain(ourBrain == 2) = 1;
        ourBrain(ourBrain == 3 | ourBrain == 4) = 2;
        ourBrain(ourBrain == 5 | ourBrain == 6) = 3;
        ourBrain(ourBrain == 7) = 4;
        angleBrain = ourBrain;
        angleBrain(angleBrain > 0) = 0;
        angleBrain(ourBrain == 1) = roiOri{1,1}(3,:);
        angleBrain(ourBrain == 2) = roiOri{1,2}(3,:);
        angleBrain(ourBrain == 3) = roiOri{1,3}(3,:);
        angleBrain(ourBrain == 4) = roiOri{1,4}(3,:);
        angleBrain = angleBrain / pi *180;
        % angleBrain(angleBrain > 0) = 180 - angleBrain(angleBrain > 0);
        angleBrain(angleBrain < 0) = NaN;
        angleBrain(angleBrain == 0) = NaN;

        %hist(angleBrain);

        cursaveName = [saveName{curcond}, '_sub', num2str(isub), '.mat'];
        save(cursaveName, 'angleBrain');

        %% save afni
        %clear all;
        %cd '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/brainVolume'
        %load("angleBrain_new.mat");
        % size(angleBrain)

        %3dinfo betas_session01.nii.gz
        %3dcalc -a betas_session01.nii.gz[1] -expr a -prefix oneBeta
        % command = ['3dinfo betas_session01_sub', num2str(sub), '.nii.gz'];
        % system(command);
        % command = ['3dcalc -a betas_session01_sub', num2str(sub), '.nii.gz[1] -expr a -prefix oneBeta_sub', num2str(sub)];
        % system(command);

        currOneBeta = ['oneBeta_sub', num2str(isub), '+orig'];
        [err,V,Info] = BrikLoad(currOneBeta);
        Info.RootName = ['angleBrain_sub', num2str(isub), '+orig'];
        opt.Prefix = [saveName{curcond}, '_sub', num2str(isub)];
        WriteBrik(angleBrain,Info,opt);

        % compare original and new
        % roiOri_New = roiOri;
        % origfileName = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/prfsample/voxModelPref_sub', num2str(sub), '.mat'];
        % load(origfileName);
        %
        % original = roiOri{1,1}(3,:)';
        % new = roiOri_New{1,1}(3,:)';
        % [r,p] = corr(original, new);
        %
        % corrList = [corrList; r, p];

        %% save nifti
        % load(['angleBrain_sub', num2str(isub), '.mat']);

        niftiwrite(angleBrain,[saveName{curcond}, '_sub', num2str(isub),'.nii']);
        info_old = niftiinfo(['betas_session01_sub', num2str(isub),'.nii.gz']);
        info_new = niftiinfo([saveName{curcond}, '_sub', num2str(isub),'.nii']);

        info_new.PixelDimensions = [1.8, 1.8, 1.8];
        info_new.TransformName = info_old.TransformName;
        info_new.SpatialDimension = info_old.SpatialDimension;
        info_new.Transform = info_old.Transform;
        info_new.Qfactor = info_old.Qfactor;
        info_new.AuxiliaryFile = info_old.AuxiliaryFile;
        info_new.raw.pixdim = info_old.raw.pixdim;
        info_new.raw.aux_file = info_old.raw.aux_file;
        info_new.raw.sform_code = info_old.raw.sform_code;
        info_new.raw.srow_x = info_old.raw.srow_x;
        info_new.raw.srow_y = info_old.raw.srow_y;
        info_new.raw.srow_z = info_old.raw.srow_z;

        niftiwrite(angleBrain,[saveName{curcond}, '_sub', num2str(isub),'.nii'],info_new);
    end
end
%% R2
nrois = 4;
saveName = {'old','ori','control'};

% change condition!
for curcond = 1:3

    for isub = 1:8
        fileName = [folderDir, subfolderName{curcond}, 'voxModelPref_sub', num2str(isub), '.mat'];
        load(fileName,'roiNsdOriR2', 'visRoiData');
        allNsdOriR2 = cell(1,nrois);
        for iroi=1:nrois
            allNsdOriR2{iroi} = [allNsdOriR2{iroi} roiNsdOriR2{iroi}];
        end

        ourBrain = visRoiData;
        ourBrain(ourBrain == 2) = 1;
        ourBrain(ourBrain == 3 | ourBrain == 4) = 2;
        ourBrain(ourBrain == 5 | ourBrain == 6) = 3;
        ourBrain(ourBrain == 7) = 4;
        angleBrain = ourBrain;
        angleBrain(angleBrain > 0) = 0;
        angleBrain(ourBrain == 1) = allNsdOriR2{1}(3,:);
        angleBrain(ourBrain == 2) = allNsdOriR2{2}(3,:);
        angleBrain(ourBrain == 3) = allNsdOriR2{3}(3,:);
        angleBrain(ourBrain == 4) = allNsdOriR2{4}(3,:);
        angleBrain(angleBrain < 0) = NaN;
        angleBrain(angleBrain == 0) = NaN;

        niftiwrite(angleBrain,[saveName{curcond}, 'R2_sub', num2str(isub),'.nii']);
        info_old = niftiinfo(['betas_session01_sub', num2str(isub),'.nii.gz']);
        info_new = niftiinfo([saveName{curcond}, 'R2_sub', num2str(isub),'.nii']);

        info_new.PixelDimensions = [1.8, 1.8, 1.8];
        info_new.TransformName = info_old.TransformName;
        info_new.SpatialDimension = info_old.SpatialDimension;
        info_new.Transform = info_old.Transform;
        info_new.Qfactor = info_old.Qfactor;
        info_new.AuxiliaryFile = info_old.AuxiliaryFile;
        info_new.raw.pixdim = info_old.raw.pixdim;
        info_new.raw.aux_file = info_old.raw.aux_file;
        info_new.raw.sform_code = info_old.raw.sform_code;
        info_new.raw.srow_x = info_old.raw.srow_x;
        info_new.raw.srow_y = info_old.raw.srow_y;
        info_new.raw.srow_z = info_old.raw.srow_z;

        niftiwrite(angleBrain,[saveName{curcond}, 'R2_sub', num2str(isub),'.nii'],info_new);
    end
end