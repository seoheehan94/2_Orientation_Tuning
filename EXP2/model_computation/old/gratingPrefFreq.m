function prefFreqNum = gratingPrefFreq(fullCoef,modelOriEnergy);
        %         global prefAnalysis
        [numFreqs, numAngles, numLevels, numOrientations] = size(modelOriEnergy);
        numvox = size(fullCoef,1);
        keyboard;
        coefMat = reshape(fullCoef,numvox,numLevels,numOrientations);%vox x levels x orientations
        coefVec = reshape(fullCoef,numvox,numLevels*numOrientations);%vox x coefficients
        modelEnergy = reshape(modelOriEnergy,numFreqs*numAngles,numLevels*numOrientations);
        modelEnergy = modelEnergy';%(ilev*orientation, isf*iangle)
        
        voxGratingResp = coefVec*modelEnergy;%vox,isf*iangle
        switch prefAnalysis
            case 1
                % 1 - simple max
                [m maxInd] = max(voxGratingResp,[],2);%vox
                [prefFreqNum, prefAngleNum] = ind2sub([numFreqs, numAngles],maxInd);%vox
            case 2
                % 2 - average across gratings orientation, and then max
                voxFreqResp = squeeze(mean(reshape(voxGratingResp,numvox, numFreqs, numAngles),3));
                [m prefFreqNum] = max(voxFreqResp,[],2);%ivox
            case 3
                % 3 - average across gratings orientation, and then weighted mean
                voxFreqResp = squeeze(mean(reshape(voxGratingResp,numvox, numFreqs, numAngles),3));
                %subtract minimum response
                voxFreqResp = voxFreqResp - min(voxFreqResp,[],2);
                prefFreqNum = (voxFreqResp * [1:numFreqs]')./sum(voxFreqResp,2);
        end
    end