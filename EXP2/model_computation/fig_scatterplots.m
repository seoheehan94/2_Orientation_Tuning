% fig2.m
%
% associated with the following publication: Roth, ZN, Kay, K, and Merriam, EP (2022).
% Massive natural scene sampling reveals reliable coarse-scale orientation tuning in human V1
% DOI:
%
%   usage: fig2()
%   by: zvi roth
%   date: 7/29/2022
%   purpose: Plot scatter plots for Figure 2
%   uses files created by: getVoxPref.m

close all
clear all
tic

figRoi=1;

isplit=3;

toSavePdf = 0;
imgFormat = 'jpg';
%subjects = [1:8];
subjects = [1:8];

ifig=0;
nrois = 4;
%nrois = 1;

imgScaling = 0.5;
global interpSz; interpSz= 714*imgScaling;
global backgroundSz; backgroundSz= 1024*imgScaling;
global degPerPix; degPerPix = 8.4/(interpSz*imgScaling);

nhistbins = 30;
histAlpha = 0.5;

eccMin = 0.1;
eccMax = 10;
nbins = 20;
binBorders = logspace(log10(eccMin),log10(eccMax),nbins+1);
nbins = length(binBorders)-1;
for i=2:length(binBorders)
    binCenters(i-1) = (binBorders(i)+binBorders(i-1))/2;
end
binBorders;

%bins for pRF angle
angMin = 0;
angMax = 2*pi;
angBinBorders = linspace(angMin,angMax,nbins+1);
nbins = length(angBinBorders)-1;
for i=2:length(angBinBorders)
    angBinCenters(i-1) = (angBinBorders(i)+angBinBorders(i-1))/2;
end
angBinBorders;

%bins for R^2
rMin = 0;
rMax = 1;
rBinBorders = linspace(rMin,rMax,nbins+1);
nbins = length(rBinBorders)-1;
for i=2:length(rBinBorders)
    rBinCenters(i-1) = (rBinBorders(i)+rBinBorders(i-1))/2;
end
rBinBorders;

%scatter parameters
markersize = 1;
edgeAlpha = 0.3;%0.07
markerColor = [0 0 0];
prfThresh = 0;

prffolder_photoSP = '/bwdata/NSDData/Seohee/Orientation/prfsample_photoSP/';
prffolder_ldSP = '/bwdata/NSDData/Seohee/Orientation/prfsample_ldSP/';
prffolder_contour = '/bwdata/NSDData/Seohee/Orientation/prfsample_contour/';
folder_list= {prffolder_photoSP, prffolder_ldSP, prffolder_contour};
figFolder = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/figures/'];

allOri = cell(3,nrois);
allLevFull = cell(3,nrois);
allPrfR2 = cell(3,nrois);
allNsdOriCorr = cell(3,nrois);
allPrfX = cell(3,nrois);
allPrfY = cell(3,nrois);
allPrfEcc = cell(3,nrois);
allPrfAng = cell(3,nrois);
allNsdOriR2 = cell(3,nrois);

for condition = 1:3 %1=photoSP, 2=ldSP, 3=contour
    for isub=1:length(subjects)
        subnum = subjects(isub);

        load([folder_list{condition} 'voxModelPref_sub' num2str(isub) '.mat'],'allRoiPrf',...
            'roiOri','roiNsdOriCorr','roiNsdOriR2',...
            'visRoiData','roiNames','combinedRoiNames','roiInd','prefAnalysis');

        subAnalysis(isub) = prefAnalysis;
        %     subNsdSynthImprov_corr(isub,:,:) = nsdSynthImprov_corr;
        %     subNsdSynthImprov_pval(isub,:,:) = nsdSynthImprov_pval;
        for iroi=1:nrois
            allPrfX{condition, iroi} = [allPrfX{condition, iroi}; allRoiPrf{iroi}.x];
            allPrfY{condition, iroi} = [allPrfY{condition, iroi}; allRoiPrf{iroi}.y];
            allPrfEcc{condition, iroi} = [allPrfEcc{condition, iroi}; allRoiPrf{iroi}.ecc];
            allPrfAng{condition, iroi} = [allPrfAng{condition, iroi}; allRoiPrf{iroi}.ang];
            allPrfR2{condition, iroi} = [allPrfR2{condition, iroi}; allRoiPrf{iroi}.r2];
            allOri{condition, iroi} =  [allOri{condition, iroi} roiOri{iroi}];
            allNsdOriCorr{condition, iroi} = [allNsdOriCorr{condition, iroi} roiNsdOriCorr{iroi}];
            allNsdOriR2{condition, iroi} = [allNsdOriR2{condition, iroi} roiNsdOriR2{iroi}];
            %         allSynthCorr{iroi} = [allSynthCorr{iroi} roiSynthCorr{iroi}];
            %         allSynthOriCorr{iroi} = [allSynthOriCorr{iroi} roiSynthOriCorr{iroi}];
            %         allSynthOri{iroi} = [allSynthOri{iroi}; roiSynthOri{iroi}'];
            %         allSynthLevVig{iroi} = [allSynthLevVig{iroi} roiSynthLevVig{iroi}'];
            %         allSynthLevFull{iroi} = [allSynthLevFull{iroi} roiSynthLevFull{iroi}'];
            %         allImprovCorr{iroi} = [allImprovCorr{iroi} roiOriDeviation{iroi}];
        end
    end
end

% ifig=ifig+1; h=figure(ifig); clf;
% rows=2;
% cols = 4;
% isubplot=0;


% constColor = [133,149,225]/255;
% fullColor = [224,123,145]/255;%pinkish red
r2Color = [156,222,214]/255;
linewidth=2;

% Base color
color_photoSP = [133,149,225]/255;
color_ldSP = [225, 99, 116] / 255;
color_contour = [141,213,147]/255;

iroi=figRoi;
%%
%histogram of R2 for all three models
%isubplot=isubplot+1; subplot(rows,cols,isubplot);
histogram(allNsdOriR2{2, iroi}(isplit,:),nhistbins,'faceColor',color_lsSP,'faceAlpha',histAlpha,'Normalization','probability'); hold on;
histogram(allNsdOriR2{3, iroi}(isplit,:),nhistbins,'faceColor',color_contour,'faceAlpha',histAlpha,'Normalization','probability'); hold on;
histogram(allNsdOriR2{1, iroi}(isplit,:),nhistbins,'faceColor',color_photoSP,'faceAlpha',histAlpha,'Normalization','probability'); hold on;

xlabel('\itR^2');
ylabel('\itproportion of voxels');
legend('\itLine Drawing-Steerable Pyramid','\itContour', '\itPhoto-Steerable Pyramid');
axis square
mean(allNsdOriR2{1,iroi}(isplit,:))
mean(allNsdOriR2{2,iroi}(isplit,:))
mean(allNsdOriR2{3,iroi}(isplit,:))

%pR2improvMedian = ranksum(allNsdOriR2{1,iroi}(isplit,:),allNsdOriR2{2,iroi}(isplit,:))
[corrR2, pCorrR2] = corr(allNsdOriCorr{1,iroi}(isplit,:)',allNsdOriCorr{2,iroi}(isplit,:)')
[corrR2, pCorrR2] = corr(allNsdOriCorr{2,iroi}(isplit,:)',allNsdOriCorr{3,iroi}(isplit,:)')
[corrR2, pCorrR2] = corr(allNsdOriCorr{1,iroi}(isplit,:)',allNsdOriCorr{3,iroi}(isplit,:)')

%% scatter plot of full R2 vs. constrained R2
% isubplot=isubplot+1; subplot(rows,cols,isubplot);
scatter(allNsdOriR2{1, iroi}(isplit,:),allNsdOriR2{2, iroi}(isplit,:),markersize,markerColor,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',edgeAlpha);
xlim([0 0.2]); ylim([0 0.2]);
axis square
xlabel('\itPhoto-Steerable Pyramid R^2');
ylabel('\itLine Drawing-Steerable Pyramid R^2');
box on

scatter(allNsdOriR2{1, iroi}(isplit,:),allNsdOriR2{3, iroi}(isplit,:),markersize,markerColor,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',edgeAlpha);
xlim([0 0.2]); ylim([0 0.2]);
axis square
xlabel('\itPhoto-Steerable Pyramid R^2');
ylabel('\itContour R^2');
box on

scatter(allNsdOriR2{2, iroi}(isplit,:),allNsdOriR2{3, iroi}(isplit,:),markersize,markerColor,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',edgeAlpha);
xlim([0 0.2]); ylim([0 0.2]);
axis square
xlabel('\itLine Drawing-Steerable Pyramid R^2');
ylabel('\itContour R^2');
box on

%% pRF R2
% Photo-Steerable Pyramid R2  vs. pRF R2
% isubplot=isubplot+1; subplot(rows,cols,isubplot);
scatter(allPrfR2{1, iroi}./100,allNsdOriR2{1, iroi}(isplit,:),markersize,markerColor,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',edgeAlpha);
xlim([0 1]);  ylim([0 0.2]); %ylim([0 1]);
axis square
xlabel('\itpRF R^2');
ylabel('\itPhoto-Steerable Pyramid R^2');
box on
hold on
temp = allNsdOriR2{1, iroi}(isplit,:);
for ibin=1:nbins
    binData(ibin) = mean(temp(allPrfR2{1, iroi}./100>rBinBorders(ibin) & allPrfR2{1, iroi}./100<=rBinBorders(ibin+1)));
end
p=plot(rBinCenters,binData,'color',color_photoSP,'linestyle','-','linewidth',linewidth);
[corrConstR2PrfR2, pCorrConstR2PrfR2] = corr(allNsdOriR2{1, iroi}(isplit,:)',allPrfR2{1, iroi}./100,'type','Pearson')

% Line Drawing-Steerable Pyramid R2  vs. pRF R2
scatter(allPrfR2{2, iroi}./100,allNsdOriR2{2, iroi}(isplit,:),markersize,markerColor,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',edgeAlpha);
xlim([0 1]);  ylim([0 0.2]); %ylim([0 1]);
axis square
xlabel('\itpRF R^2');
ylabel('\itLine Drawing-Steerable Pyramid R^2');
box on
hold on
temp = allNsdOriR2{2, iroi}(isplit,:);
for ibin=1:nbins
    binData(ibin) = mean(temp(allPrfR2{2, iroi}./100>rBinBorders(ibin) & allPrfR2{2, iroi}./100<=rBinBorders(ibin+1)));
end
p=plot(rBinCenters,binData,'color',color_ldSP,'linestyle','-','linewidth',linewidth);
[corrConstR2PrfR2, pCorrConstR2PrfR2] = corr(allNsdOriR2{2, iroi}(isplit,:)',allPrfR2{2, iroi}./100,'type','Pearson')

% Contour R2  vs. pRF R2
scatter(allPrfR2{3, iroi}./100,allNsdOriR2{3, iroi}(isplit,:),markersize,markerColor,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',edgeAlpha);
xlim([0 1]);  ylim([0 0.2]); %ylim([0 1]);
axis square
xlabel('\itpRF R^2');
ylabel('\itLine Drawing-Steerable Pyramid R^2');
box on
hold on
temp = allNsdOriR2{3, iroi}(isplit,:);
for ibin=1:nbins
    binData(ibin) = mean(temp(allPrfR2{3, iroi}./100>rBinBorders(ibin) & allPrfR2{3, iroi}./100<=rBinBorders(ibin+1)));
end
p=plot(rBinCenters,binData,'color',color_contour,'linestyle','-','linewidth',linewidth);
[corrConstR2PrfR2, pCorrConstR2PrfR2] = corr(allNsdOriR2{3, iroi}(isplit,:)',allPrfR2{3, iroi}./100,'type','Pearson')

%% improvement  vs. pRF R2
% ldSP - photoSP vs. pRF R2
improv = allNsdOriR2{2, iroi}(isplit,:) - allNsdOriR2{1, iroi}(isplit,:);
scatter(allPrfR2{1, iroi}./100,allNsdOriR2{2, iroi}(isplit,:) - allNsdOriR2{1, iroi}(isplit,:),markersize,markerColor,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',edgeAlpha);
xlim([0 1]);
axis square
xlabel('\itpRF R^2');
ylabel('\itldSP R^2 - photoSP R^2');
box on
hold on
temp = improv;
for ibin=1:nbins
    binData(ibin) = mean(temp(allPrfR2{1, iroi}./100>rBinBorders(ibin) & allPrfR2{1, iroi}./100<=rBinBorders(ibin+1)));
end
p=plot(rBinCenters,binData,'color',r2Color,'linestyle','-','linewidth',linewidth);
ylim([-0.01 0.08]);
[corrImprovPrfR2, pCorrImprovPrfR2] = corr(temp',allPrfR2{1, iroi}./100,'type','Pearson')

% ldSP - contour vs. pRF R2
improv = allNsdOriR2{2, iroi}(isplit,:) - allNsdOriR2{3, iroi}(isplit,:);
scatter(allPrfR2{1, iroi}./100,allNsdOriR2{2, iroi}(isplit,:) - allNsdOriR2{3, iroi}(isplit,:),markersize,markerColor,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',edgeAlpha);
xlim([0 1]);
axis square
xlabel('\itpRF R^2');
ylabel('\itldSP R^2 - contour R^2');
box on
hold on
temp = improv;
for ibin=1:nbins
    binData(ibin) = mean(temp(allPrfR2{1, iroi}./100>rBinBorders(ibin) & allPrfR2{1, iroi}./100<=rBinBorders(ibin+1)));
end
p=plot(rBinCenters,binData,'color',r2Color,'linestyle','-','linewidth',linewidth);
ylim([-0.01 0.08]);
[corrImprovPrfR2, pCorrImprovPrfR2] = corr(temp',allPrfR2{1, iroi}./100,'type','Pearson')

% contour - photoSP vs. pRF R2
improv = allNsdOriR2{3, iroi}(isplit,:) - allNsdOriR2{1, iroi}(isplit,:);
scatter(allPrfR2{1, iroi}./100,allNsdOriR2{3, iroi}(isplit,:) - allNsdOriR2{1, iroi}(isplit,:),markersize,markerColor,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',edgeAlpha);
xlim([0 1]);
axis square
xlabel('\itpRF R^2');
ylabel('\itcontour R^2 - photoSP R^2');
box on
hold on
temp = improv;
for ibin=1:nbins
    binData(ibin) = mean(temp(allPrfR2{1, iroi}./100>rBinBorders(ibin) & allPrfR2{1, iroi}./100<=rBinBorders(ibin+1)));
end
p=plot(rBinCenters,binData,'color',r2Color,'linestyle','-','linewidth',linewidth);
ylim([-0.01 0.08]);
[corrImprovPrfR2, pCorrImprovPrfR2] = corr(temp',allPrfR2{1, iroi}./100,'type','Pearson')


%% improvement vs. pRF eccentricity
% ldSP - photoSP vs. pRF eccentricity
figure;
improv = allNsdOriR2{2,iroi}(isplit,:) - allNsdOriR2{1,iroi}(isplit,:);
temp = improv;
scatter(allPrfEcc{1,iroi},temp,markersize,markerColor,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',edgeAlpha);
axis square
xlabel('\itpRF eccentricity (deg)');
ylabel('\itldSP R^2 - photoSP R^2');
box on
hold on
for ibin=1:nbins
    binData(ibin) = mean(temp(allPrfR2{1,iroi}>prfThresh &  allPrfEcc{1,iroi}>=binBorders(ibin) & allPrfEcc{iroi}<binBorders(ibin+1)));
end
p=plot(binCenters,binData,'color',r2Color,'linestyle','-','linewidth',linewidth);
xlim([0 10]);
ylim([-0.01 0.08]);

% ldSP - contour vs. pRF eccentricity
figure;
improv = allNsdOriR2{2,iroi}(isplit,:) - allNsdOriR2{3,iroi}(isplit,:);
temp = improv;
scatter(allPrfEcc{1,iroi},temp,markersize,markerColor,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',edgeAlpha);
axis square
xlabel('\itpRF eccentricity (deg)');
ylabel('\itldSP R^2 - contour R^2');
box on
hold on
for ibin=1:nbins
    binData(ibin) = mean(temp(allPrfR2{1,iroi}>prfThresh &  allPrfEcc{1,iroi}>=binBorders(ibin) & allPrfEcc{iroi}<binBorders(ibin+1)));
end
p=plot(binCenters,binData,'color',r2Color,'linestyle','-','linewidth',linewidth);
xlim([0 10]);
ylim([-0.01 0.08]);

% contour - photoSP vs. pRF eccentricity
figure;
improv = allNsdOriR2{3,iroi}(isplit,:) - allNsdOriR2{1,iroi}(isplit,:);
temp = improv;
scatter(allPrfEcc{1,iroi},temp,markersize,markerColor,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',edgeAlpha);
axis square
xlabel('\itpRF eccentricity (deg)');
ylabel('\itldSP R^2 - photoSP R^2');
box on
hold on
for ibin=1:nbins
    binData(ibin) = mean(temp(allPrfR2{1,iroi}>prfThresh &  allPrfEcc{1,iroi}>=binBorders(ibin) & allPrfEcc{iroi}<binBorders(ibin+1)));
end
p=plot(binCenters,binData,'color',r2Color,'linestyle','-','linewidth',linewidth);
xlim([0 10]);
ylim([-0.01 0.08]);

%% scatter plot: pRF R^2 vs. pRF eccentricity
scatter(allPrfEcc{1, iroi},allPrfR2{1, iroi}/100,markersize,markerColor,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',edgeAlpha);
xlim([0 10]);
axis square
xlabel('\itpRF eccentricity (deg)');
ylabel('\itpRF R^2');
box on
hold on
for ibin=1:nbins
    binData(ibin) = mean(allPrfR2{1, iroi}(allPrfR2{1,iroi}>prfThresh & allPrfEcc{1,iroi}>=binBorders(ibin) & allPrfEcc{iroi}<binBorders(ibin+1))/100);
end
plot(binCenters,binData,'color',r2Color,'linestyle','-','linewidth',linewidth);
ylim([0 1]);
xlim([0 10]);

%%
set(gcf,'position',[150 180 940 420]);
if toSavePdf
    savepdf(h, [figFolder 'fig2.pdf']);
    saveas(h, [figFolder 'fig2.' imgFormat]);
end

toc
