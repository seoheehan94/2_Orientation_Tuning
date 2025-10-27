%% prfSampleModel

addpath(genpath('/home/hanseohe/Documents/GitHub/3rdParty/stimulusVignetting'))

for sub = 1:8
    for region = 1:4
        fprintf('%s. %d. %d ...\n','prfSampleModel_photoSP',sub,region);
        prfSampleModel_photoSP(sub,region);
    end
end

for sub = 1:8
    for region = 1:4
        fprintf('%s. %d. %d ...\n','prfSampleModel_ldSP',sub,region);
        prfSampleModel_ldSP(sub,region);
   
    end
end
for sub = 1:8
    for region = 1:4
        fprintf('%s. %d. %d ...\n','prfSampleModel_contour',sub,region);
        prfSampleModel_contour(sub,region);
    end
end

%% regressPrfSplit
for sub = 1:8
    for region = 1:4
        fprintf('%s. %d. %d ...\n','regressPrfSplit_photoSP',sub, region);
        regressPrfSplit_photoSP(sub, region);
    end
end

for sub = 1:8
    for region = 1:4
        fprintf('%s. %d. %d ...\n','regressPrfSplit_ldSP',sub, region);
        regressPrfSplit_ldSP(sub, region);
    end
end

for sub = 1:8
    for region = 1:4
        fprintf('%s. %d. %d ...\n','regressPrfSplit_contour',sub, region);
        regressPrfSplit_contour(sub, region);
    end
end

%% getVoxPref
for sub = 1:8
    fprintf('%s. %d. %d ...\n','getVoxPref_photoSP',sub, 1);
    getVoxPref(sub,4, 1);
end
for sub = 1:8
    fprintf('%s. %d. %d ...\n','getVoxPref_ldSP',sub, 2);
    getVoxPref(sub,4, 2);
end
for sub = 1:8
    fprintf('%s. %d. %d ...\n','getVoxPref_contour',sub, 3);
    getVoxPref(sub,4, 3);
end