%analyze_orientation.m
% get voxel preference from model weights and create brainVolume
%   uses files created by: regressPrfSplit.m
%   creates files used by:
clear all;
savefolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/brainVolume_regress/maxmin';
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4'};
combinedRoiNames = {'V1','V2','V3','hV4'};
condition = {'old', 'ori', 'control'};
prffolder = {'prfsample_maxmin', 'prfsample_maxmin_Ori', 'prfsample_maxmin_control'};
subList = {'sub1', 'sub2','sub3', 'sub4','sub5', 'sub6', 'sub7', 'sub8'};
imgType = {'top20k', 'bottom20k'};
totalR2OriSplit = struct;
totalaicOriSplit = struct;
totalbicOriSplit = struct;


%% 1.  Coef voxel preference %%
for curtype =1:2
    for isub = 1:8
        for con = 1:3

            curPrf = ['/bwdata/NSDData/Seohee/Orientation/', prffolder{con}, '/'];
            if con == 1 || con == 3
                bandpassStr = '_bandpass1to7';
            elseif con == 2
                bandpassStr = '_bandpass1to1';
            end

            % cd('/home/hanseohe/Documents/GitHub/nsdOtopy');
            fprintf('Current loop - curtype: %d, isub: %d, con: %d\n', curtype, isub, con);
            clearvars -except isub con subList curtype curPrf imgType totalaicOriSplit totalbicOriSplit totalR2OriSplit roiNames combinedRoiNames curPrf prffolder savefolder bandpassStr condition
            allR2OriSplit = struct;

            for visualRegion = 1:4
                thisfield = combinedRoiNames{visualRegion};
                load(fullfile(curPrf,[imgType{curtype}, 'regressPrfSplit' bandpassStr '_v' num2str(visualRegion) '_sub' num2str(isub)  '.mat']), ...
                    'nsd','rois','roiPrf','nsplits');

                if length(rois)>1 %combine across ventral and dorsal ROIs
                    oldNsd = nsd;
                    nsd.voxOriCoef{1} = [];
                    nsd.roiInd{1} = [];
                    nsd.r2oriSplit{1} = [];
                    nsd.aicOriSplit{1} = [];
                    nsd.bicOriSplit{1} = [];

                    for iroi=1:length(rois)
                        nsd.voxOriCoef{1} = cat(2,nsd.voxOriCoef{1},oldNsd.voxOriCoef{iroi});
                        nsd.roiInd{1} = cat(1,nsd.roiInd{1}, oldNsd.roiInd{iroi});
                        nsd.r2oriSplit{1} = cat(2,nsd.r2oriSplit{1},oldNsd.r2oriSplit{iroi});
                        nsd.aicOriSplit{1} = cat(2,nsd.aicOriSplit{1},oldNsd.aicOriSplit{iroi});
                        nsd.bicOriSplit{1} = cat(2,nsd.bicOriSplit{1},oldNsd.bicOriSplit{iroi});
                    end

                    oldPrf = roiPrf; clear roiPrf;
                    roiPrf{1}.ecc=[];
                    roiPrf{1}.ang=[];
                    roiPrf{1}.sz=[];
                    %         roiPrf{1}.exponent=[];
                    %         roiPrf{1}.gain=[];
                    roiPrf{1}.r2=[];
                    roiPrf{1}.x=[];
                    roiPrf{1}.y=[];
                    for iroi=1:length(rois)
                        roiPrf{1}.ecc = cat(1,roiPrf{1}.ecc,oldPrf{iroi}.ecc);
                        roiPrf{1}.ang = cat(1,roiPrf{1}.ang,oldPrf{iroi}.ang);
                        roiPrf{1}.sz = cat(1,roiPrf{1}.sz,oldPrf{iroi}.sz);
                        %             roiPrf{1}.exponent = cat(1,roiPrf{1}.exponent,oldPrf{iroi}.exponent);
                        %             roiPrf{1}.gain = cat(1,roiPrf{1}.gain,oldPrf{iroi}.gain);
                        roiPrf{1}.r2 = cat(1,roiPrf{1}.r2,oldPrf{iroi}.r2);
                        roiPrf{1}.x = cat(1,roiPrf{1}.x,oldPrf{iroi}.x);
                        roiPrf{1}.y = cat(1,roiPrf{1}.y,oldPrf{iroi}.y);
                    end
                    rois = 1;
                end

                % AVERAGE SPLITS
                nsd.voxOriCoef{1}(nsplits+1,:,:) = mean(nsd.voxOriCoef{1},1);
                nsd.r2oriSplit{1}(nsplits+1,:) = mean(nsd.r2oriSplit{1},1);
                nsd.aicOriSplit{1}(nsplits+1,:) = mean(nsd.aicOriSplit{1},1);
                nsd.bicOriSplit{1}(nsplits+1,:) = mean(nsd.bicOriSplit{1},1);

                allR2OriSplit.(thisfield) = nsd.r2oriSplit{1}(3,:);
                allaicOriSplit.(thisfield) = nsd.aicOriSplit{1}(3,:);
                allbicOriSplit.(thisfield) = nsd.bicOriSplit{1}(3,:);
                allRoiPrf{visualRegion} = roiPrf{1};
                roiNsdOriR2{visualRegion} = nsd.r2oriSplit{1};

                %% get regress angle
                prefAngle =[];
                fullCoef = squeeze(nsd.voxOriCoef{1}(3,:,1:end-1));
                numvox = size(fullCoef,1);
                numLevels = 7;
                
                if con == 1 || con ==3
                    coefMat = reshape(fullCoef,numvox,numLevels,numOrientations);%vox x levels x orientations
                    coefMat_meansf = squeeze(mean(coefMat, 2));
                elseif con ==2
                    coefMat_meansf = fullCoef;
                end
                coefMat_meansf = coefMat_meansf - min(coefMat_meansf,[],2);
                theta = linspace(0,2*pi,numOrientations+1);%for circular calculation
                theta = theta(1:end-1);
                for ivox=1:numvox
                    prefAngle(ivox) = circ_mean(theta',coefMat_meansf(ivox,:)');
                end
                prefAngle = mod(prefAngle,2*pi);%from [-pi, pi] to [0 2pi]
                prefAngle = prefAngle./2;%range 0 to pi.
                prefAngle = prefAngle / pi *180;


                % allMeanCoef.(thisfield) = prefAngle;
                allMeanCoef{visualRegion} = prefAngle;

                % saveName = [curPrf, imgType{curtype}, 'voxOriCoef_sub', num2str(isub), '.mat'];
                % save(saveName, 'allMeanCoef', 'allRoiPrf','roiNsdOriR2');

                

            end
            
            %% total values of R2, aic, bic
            totalR2OriSplit.(condition{con}).(subList{isub}) = allR2OriSplit;
            totalaicOriSplit.(condition{con}).(subList{isub}) = allaicOriSplit;
            totalbicOriSplit.(condition{con}).(subList{isub}) = allbicOriSplit;

            %% make a brain volume
            % save all ROIs to create overlay
            roifolder = ['/bwdata/NSDData/nsddata/ppdata/subj0' num2str(isub) '/func1pt8mm/'];
            visualRoisFile = fullfile(roifolder,'roi/prf-visualrois.nii.gz');%V1v, V1d, V2v, V2d, V3v, V3d, and hV4
            visRoiData = niftiread(visualRoisFile);

            ourBrain = visRoiData;
            ourBrain(ourBrain == 2) = 1;
            ourBrain(ourBrain == 3 | ourBrain == 4) = 2;
            ourBrain(ourBrain == 5 | ourBrain == 6) = 3;
            ourBrain(ourBrain == 7) = 4;


            % make a brain volume
            newBrain = ourBrain;
            newBrain(newBrain > 0) = 0;
            for visualRegion = 1:4
                curOurBrain = ourBrain;
                % if visualRegion == 2
                %     curOurBrain(visRoiData == 3 | visRoiData == 4) = 2;
                % elseif visualRegion == 3
                %     curOurBrain(visRoiData == 5 | visRoiData == 6) = 3;
                % end
                thisfield = combinedRoiNames{visualRegion};
                newBrain(curOurBrain == visualRegion) = allMeanCoef{visualRegion};
            end
            newBrain(newBrain <= 0) = -1;

            for visualRegion = 1:4
                curOurBrain = ourBrain;
                % if visualRegion == 2
                %     curOurBrain(visRoiData == 3 | visRoiData == 4) = 2;
                % elseif visualRegion == 3
                %     curOurBrain(visRoiData == 5 | visRoiData == 6) = 3;
                % end
                curNewBrain = curOurBrain;
                curNewBrain(curOurBrain ~= visualRegion) = -1;
                thisfield = combinedRoiNames{visualRegion};

                curNewBrain(curOurBrain == visualRegion) = allMeanCoef{visualRegion};

                newBrainbyROI(:,:,:,visualRegion) =curNewBrain;
            end

            % save(fullfile(savefolder, [condition{con}, 'Brain_sub', num2str(isub), '.mat']), 'newBrain');
            % save(fullfile(savefolder, [condition{con}, 'BrainbyROI_sub', num2str(isub), '.mat']), 'newBrainbyROI');

            
            % R2
            r2Brain = ourBrain;
            r2Brain(ourBrain == 0 | ourBrain == -1) = NaN;
            r2Brain(ourBrain > 0) = 0;
            for visualRegion = 1:4
                thisfield = combinedRoiNames{visualRegion};
                r2Brain(ourBrain == visualRegion) = allR2OriSplit.(thisfield);
            end

            % %% save afni file
            % % size(newBrain)
            % %3dinfo betas_session01.nii.gz
            % %3dcalc -a betas_session01.nii.gz[1] -expr a -prefix oneBeta
            % % command = ['3dinfo betas_session01_sub', num2str(isub), '.nii.gz'];
            % % system(command);
            % % command = ['3dcalc -a betas_session01_sub', num2str(isub), '.nii.gz[1] -expr a -prefix oneBeta_sub', num2str(isub)];
            % % system(command);
            %
            % currOneBeta = ['oneBeta_sub', num2str(isub), '+orig'];
            % [err,V,Info] = BrikLoad(currOneBeta);
            %
            % Info.RootName = ['curvBrain_sub', num2str(isub), '+orig'];
            % opt.Prefix = ['curvBrain_sub', num2str(isub)];
            % WriteBrik(newBrain,Info,opt);
            % Info.RootName = ['curvBrainbyROI_sub', num2str(isub), '+orig'];
            % opt.Prefix = ['curvBrainbyROI_sub', num2str(isub)];
            % WriteBrik(newBrainbyROI,Info,opt);

            %% save nifti
            % cd('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature/brainVolume_regress');
            % load(['oriBrain_sub', num2str(isub), '.mat']);
            % info_old = niftiinfo(['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/brainVolume/betas_session01_sub', num2str(isub),'.nii.gz']);
            % 
            % niftiwrite(newBrain,[savefolder, '/', imgType{curtype}, '_', condition{con}, 'Brain_sub', num2str(isub),'.nii']);
            % info_new = niftiinfo([savefolder, '/', imgType{curtype}, '_', condition{con}, 'Brain_sub', num2str(isub),'.nii']);
            % info_new.PixelDimensions = [1.8, 1.8, 1.8];
            % info_new.TransformName = info_old.TransformName;
            % info_new.SpatialDimension = info_old.SpatialDimension;
            % info_new.Transform = info_old.Transform;
            % info_new.Qfactor = info_old.Qfactor;
            % info_new.AuxiliaryFile = info_old.AuxiliaryFile;
            % info_new.raw.pixdim = info_old.raw.pixdim;
            % info_new.raw.aux_file = info_old.raw.aux_file;
            % info_new.raw.sform_code = info_old.raw.sform_code;
            % info_new.raw.srow_x = info_old.raw.srow_x;
            % info_new.raw.srow_y = info_old.raw.srow_y;
            % info_new.raw.srow_z = info_old.raw.srow_z;
            % niftiwrite(newBrain,[savefolder, '/', imgType{curtype}, '_', condition{con}, 'Brain_sub', num2str(isub),'.nii'], info_new);
            % 
            % niftiwrite(newBrainbyROI,[savefolder, '/', imgType{curtype}, '_', condition{con}, 'BrainbyROI_sub', num2str(isub),'.nii']);
            % info_new = niftiinfo([savefolder, '/', imgType{curtype}, '_', condition{con}, 'BrainbyROI_sub', num2str(isub),'.nii']);
            % info_new.PixelDimensions = info_old.PixelDimensions;
            % info_new.TransformName = info_old.TransformName;
            % info_new.SpatialDimension = info_old.SpatialDimension;
            % info_new.Transform = info_old.Transform;
            % info_new.Qfactor = info_old.Qfactor;
            % info_new.AuxiliaryFile = info_old.AuxiliaryFile;
            % info_new.raw.pixdim = info_old.raw.pixdim;
            % info_new.raw.aux_file = info_old.raw.aux_file;
            % info_new.raw.sform_code = info_old.raw.sform_code;
            % info_new.raw.srow_x = info_old.raw.srow_x;
            % info_new.raw.srow_y = info_old.raw.srow_y;
            % info_new.raw.srow_z = info_old.raw.srow_z;
            % niftiwrite(newBrainbyROI,[savefolder, '/', imgType{curtype}, '_', condition{con}, 'BrainbyROI_sub', num2str(isub),'.nii'],info_new);
            % 
            % niftiwrite(r2Brain,[savefolder, '/', imgType{curtype}, '_', condition{con}, 'BrainR2_sub', num2str(isub),'.nii']);
            % info_new = niftiinfo([savefolder, '/', imgType{curtype}, '_', condition{con}, 'BrainR2_sub', num2str(isub),'.nii']);
            % info_new.PixelDimensions = [1.8, 1.8, 1.8];
            % info_new.TransformName = info_old.TransformName;
            % info_new.SpatialDimension = info_old.SpatialDimension;
            % info_new.Transform = info_old.Transform;
            % info_new.Qfactor = info_old.Qfactor;
            % info_new.AuxiliaryFile = info_old.AuxiliaryFile;
            % info_new.raw.pixdim = info_old.raw.pixdim;
            % info_new.raw.aux_file = info_old.raw.aux_file;
            % info_new.raw.sform_code = info_old.raw.sform_code;
            % info_new.raw.srow_x = info_old.raw.srow_x;
            % info_new.raw.srow_y = info_old.raw.srow_y;
            % info_new.raw.srow_z = info_old.raw.srow_z;
            % niftiwrite(r2Brain,[savefolder, '/', imgType{curtype}, '_', condition{con}, 'BrainR2_sub', num2str(isub),'.nii'],info_new);

        end
    end

%% mean R2
keyboard;
allroiR2OriSplit=[];
V1R2OriSplit=[];
fieldsCon = fieldnames(totalR2OriSplit);
fieldsSub = fieldnames(totalR2OriSplit.old);
fieldsRoi = fieldnames(totalR2OriSplit.old.sub1);
for i = 1:numel(fieldsCon)
    curRoiR2OriSplit = [];
    curV1R2OriSplit = [];
    for j = 1:numel(fieldsSub)
        for k = 1:numel(fieldsRoi)
            curRoiR2OriSplit = [curRoiR2OriSplit, totalR2OriSplit.(fieldsCon{i}).(fieldsSub{j}).(fieldsRoi{k})];
        end

       curV1R2OriSplit = [curV1R2OriSplit, totalR2OriSplit.(fieldsCon{i}).(fieldsSub{j}).V1];
    end
    writematrix(curRoiR2OriSplit', ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/MaxMin/allroiR2', imgType{curtype}, '_', (fieldsCon{i}), '.csv']);
    writematrix(curV1R2OriSplit', ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/MaxMin/V1R2', imgType{curtype}, '_', (fieldsCon{i}), '.csv']);

    allroiR2OriSplit(i) = mean(curRoiR2OriSplit);
    V1R2OriSplit(i) = mean(curV1R2OriSplit);
end

save(fullfile(['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/MaxMin/', imgType{curtype}, '_', 'meanR2.mat']), "allroiR2OriSplit", "V1R2OriSplit");


%% AIC/BIC
allroiaicOriSplit=[];
V1aicOriSplit=[];
allroibicOriSplit=[];
V1bicOriSplit=[];
fieldsCon = fieldnames(totalaicOriSplit);
fieldsSub = fieldnames(totalaicOriSplit.old);
fieldsRoi = fieldnames(totalaicOriSplit.old.sub1);
for i = 1:numel(fieldsCon)
    curRoiaicOriSplit = [];
    curV1aicOriSplit = [];
    curRoibicOriSplit = [];
    curV1bicOriSplit = [];
    for j = 1:numel(fieldsSub)
        for k = 1:numel(fieldsRoi)
            curRoiaicOriSplit = [curRoiaicOriSplit, totalaicOriSplit.(fieldsCon{i}).(fieldsSub{j}).(fieldsRoi{k})];
            curRoibicOriSplit = [curRoibicOriSplit, totalbicOriSplit.(fieldsCon{i}).(fieldsSub{j}).(fieldsRoi{k})];
        end

       curV1aicOriSplit = [curV1aicOriSplit, totalaicOriSplit.(fieldsCon{i}).(fieldsSub{j}).V1];
       curV1bicOriSplit = [curV1bicOriSplit, totalbicOriSplit.(fieldsCon{i}).(fieldsSub{j}).V1];
    end
    writematrix(curRoiaicOriSplit, ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/MaxMin/allroiR2', imgType{curtype}, '_', (fieldsCon{i}), '.csv']);
    writematrix(curRoiaicOriSplit, ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/MaxMin/V1R2', imgType{curtype}, '_', (fieldsCon{i}), '.csv']);

    allroiaicOriSplit(i) = mean(curRoiaicOriSplit,"omitnan");
    V1aicOriSplit(i) = mean(curV1aicOriSplit,"omitnan");
    allroibicOriSplit(i) = mean(curRoibicOriSplit,"omitnan");
    V1bicOriSplit(i) = mean(curV1bicOriSplit,"omitnan");
end

% allroiaicOriSplit
% V1aicOriSplit
% allroibicOriSplit
% V1bicOriSplit

save(fullfile(['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/MaxMin/', imgType{curtype}, '_', 'meanAICBIC.mat']), "allroiaicOriSplit", "V1aicOriSplit", "allroibicOriSplit", "V1bicOriSplit");
end


%%
imgNumTrue = [];
for isub = 1:8

    load(fullfile(['/bwdata/NSDData/Seohee/Orientation/prfsample_maxmin_Ori/top20kregressPrfSplit_bandpass1to1_v1_sub', num2str(isub) '.mat']), ...
        'nsd');
    imgNumTrue = [imgNumTrue; size(nsd.imgNumTrue)];

end