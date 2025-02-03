% fig4.m
%
% associated with the following publication: Roth, ZN, Kay, K, and Merriam, EP (2022).
% Massive natural scene sampling reveals reliable coarse-scale orientation tuning in human V1
% DOI:
%
%   usage: fig4()
%   by: zvi roth
%   date: 7/29/2022
%   purpose: plot Figure 4
%   uses files created by: getVoxPref.m
addpath(genpath('/home/hanseohe/Documents/GitHub/mrTools'))

% close all
% clear all
tic

condition =2; %1=old, 2=ori, 3=control
curimgtype = 1;
curpairtype =2;
imgTypes = {'Top', 'Bottom'};
pairTypes = {'old_control', 'old_ori', 'control_ori'};

figRoi=1;
pvalThresh = 0.025;

toSavePdf = 1;
nbins = 20;
sigMarkersize = 15;

imgFormat = 'jpg';
subjects = [1:8];

ifig=0;
nrois = 4;
linewidth=2;
imgScaling = 0.5;
global interpSz; interpSz= 714*imgScaling;
global backgroundSz; backgroundSz= 1024*imgScaling;
global degPerPix; degPerPix = 8.4/interpSz;

nhistbins = 30;
histEdges = linspace(0,1,nhistbins);
histAlpha = 0.5;

%bins for pRF eccentricity
eccMin = 0.1;
eccMax = 10;

binBorders = logspace(log10(eccMin),log10(eccMax),nbins+1);
nbins = length(binBorders)-1;
for i=2:length(binBorders)
    binCenters(i-1) = (binBorders(i)+binBorders(i-1))/2;
end
binBorders;

binColorMean = 'm.-';
binColorMedian = 'c.-';
binColor = 'm.-';

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

% vertColor = [133,149,225]/255;
% cardColor = [141,213,147]/255;
% vertColor = [239,151,8]/255;
% cardColor = [239,151,8]/255;
% radColor = [239,151,8]/255;

% Base color
if condition == 1
    baseColor = [133,149,225]/255;
elseif condition == 2
    baseColor = [225, 99, 116] / 255;
elseif condition == 3
    baseColor = [141,213,147]/255;
end

% Variations of the base color
vertColor = baseColor * 1.2;  % Slightly lighter (scaled up to brighten)
cardColor = baseColor;         % Original color
radColor = baseColor * 0.8;    % Slightly darker (scaled down to darken)

% Ensure the values stay within [0, 1] range
vertColor = min(vertColor, 1);  % Cap values to 1
radColor = max(radColor, 0);    % Ensure values stay above 0

r2Color = [156,222,214]/255;%turquoise
surfaceAlpha = 0.1;

% Setting different line styles
lineStyleVert = ':';  % dotted line for binVert
lineStyleCard = '--'; % dashed line for binCard
lineStyleRad = '-';   % solid line for binRad

%scatter parameters
markersize = 1;
edgeAlpha = 0.3;
markerColor = [0 0 0];
prfThresh = 0;

if condition == 1
    prffolder = ['/bwdata/NSDData/Seohee/Orientation/prfsample_maxmin/'];
elseif condition == 2
    prffolder = ['/bwdata/NSDData/Seohee/Orientation/prfsample_maxmin_Ori/'];
elseif condition == 3
    prffolder = ['/bwdata/NSDData/Seohee/Orientation/prfsample_maxmin_control/'];
end
figFolder = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/figures/'];

allOri = cell(1,nrois);
% allLevVig = cell(1,nrois);
% allLevFull = cell(1,nrois);
allPrfR2 = cell(1,nrois);
% allNsdCorr = cell(1,nrois);
% allNsdOriCorr = cell(1,nrois);
% allSynthCorr = cell(1,nrois);
% allSynthOriCorr = cell(1,nrois);
allPrfX = cell(1,nrois);
allPrfY = cell(1,nrois);
allPrfEcc = cell(1,nrois);
allPrfAng = cell(1,nrois);
% allSynthOri = cell(1,nrois);
% allSynthLevVig = cell(1,nrois);
% allSynthLevFull = cell(1,nrois);
allOriDeviation = cell(1,nrois);
allVertDeviation = cell(1,nrois);
allCardDeviation = cell(1,nrois);
allNsdOriR2 = cell(1,nrois);
allNsdR2 = cell(1,nrois);
allSubInd = cell(1,nrois);

for isub=1:length(subjects)
    subnum = subjects(isub);
    load([prffolder 'indices' imgTypes{curimgtype} '_' pairTypes{curpairtype} 'voxModelPref_regress_sub' num2str(isub) '.mat']);
    % subAnalysis(isub) = prefAnalysis;
%     subNsdSynthImprov_corr(isub,:,:) = nsdSynthImprov_corr;
%     subNsdSynthImprov_pval(isub,:,:) = nsdSynthImprov_pval;
    for iroi=figRoi
        allPrfX{iroi} = [allPrfX{iroi}; allRoiPrf{iroi}.x];
        allPrfY{iroi} = [allPrfY{iroi}; allRoiPrf{iroi}.y];
        allPrfEcc{iroi} = [allPrfEcc{iroi}; allRoiPrf{iroi}.ecc];
        allPrfAng{iroi} = [allPrfAng{iroi}; allRoiPrf{iroi}.ang];
        allPrfR2{iroi} = [allPrfR2{iroi}; allRoiPrf{iroi}.r2];
        allOri{iroi} =  [allOri{iroi} roiOri{iroi}];
        % allLevVig{iroi} = [allLevVig{iroi} roiLevVig{iroi}];
        % allLevFull{iroi} = [allLevFull{iroi} roiLevFull{iroi}];
        % allNsdCorr{iroi} = [allNsdCorr{iroi} roiNsdCorr{iroi}];
        % allNsdOriCorr{iroi} = [allNsdOriCorr{iroi} roiNsdOriCorr{iroi}];
        % allNsdOriR2{iroi} = [allNsdOriR2{iroi} roiNsdOriR2{iroi}];
        % allNsdR2{iroi} = [allNsdR2{iroi} roiNsdR2{iroi}];
%         allSynthCorr{iroi} = [allSynthCorr{iroi} roiSynthCorr{iroi}];
%         allSynthOriCorr{iroi} = [allSynthOriCorr{iroi} roiSynthOriCorr{iroi}];
%         allSynthOri{iroi} = [allSynthOri{iroi}; roiSynthOri{iroi}'];
%         allSynthLevVig{iroi} = [allSynthLevVig{iroi} roiSynthLevVig{iroi}'];
%         allSynthLevFull{iroi} = [allSynthLevFull{iroi} roiSynthLevFull{iroi}'];
         allOriDeviation{iroi} = [allOriDeviation{iroi} roiOriDeviation{iroi}];
         allVertDeviation{iroi} = [allVertDeviation{iroi} roiVertDeviation{iroi}];
         allCardDeviation{iroi} = [allCardDeviation{iroi} roiCardDeviation{iroi}];
         allSubInd{iroi} = [allSubInd{iroi}; subnum*ones(size(roiOriDeviation{iroi},2),1)];
    end
end
%% SCHEMATIC
% ifig=ifig+1; h=figure(ifig); clf;
% rows=1;
% cols=3;
% isubplot=0;

% %lines
% nlines = 15;
% necc = 4;
% eccVals = linspace(1,9,necc);
% linesPerEcc = 16;
% ecc = [];
% ang = [];
% for iecc=1:necc
%     ecc = [ecc eccVals(iecc)*ones(1,nlines)./degPerPix];
%     tempang = linspace(0, 2*pi, nlines+1);
%     ang = [ang tempang(1:end-1)];
% end
% [x, y] = pol2cart(ang, ecc);
% r2 = 0.9*ones(size(x));
% 
% % vertical
% isubplot=isubplot+1;
% subplot(rows,cols,isubplot);
% vertOri = (0)*ones(size(x));
% plotOriLines(vertOri, x, y,(0.4*r2),vertColor);
% 
% % cardinal
% isubplot=isubplot+1;
% subplot(rows,cols,isubplot);
% cardOri = pi/2*ones(size(x));
% cardOri((ang>pi/4 & ang<3*pi/4) |  (ang<7*pi/4 & ang>5*pi/4)) = 0;
% plotOriLines(cardOri, x, y,(0.4*r2),cardColor);
% 
% % radial
% isubplot=isubplot+1;
% subplot(rows,cols,isubplot);
% radOri = pi/2-mod(ang,pi);
% plotOriLines(radOri, x, y,(0.4*r2),radColor);
% 
% %%
% set(gcf,'position',[150 180 420 200]);
% if toSavePdf
%     print('-painters','-dpdf',[figFolder 'ideal_schematic']);
% end
% toSavePdf = 0;
%%
ifig=ifig+1; h=figure(ifig); clf;
iroi=figRoi;
isplit=nsplits;
isubplot=0;
rows=3;
cols=2;
%% line plot: vertical, cardinal, and radial deviation vs. pRF eccentricity
isubplot=isubplot+1;
subplot(rows,cols, isubplot);
tempVert = abs(allVertDeviation{iroi}(isplit,:));
tempRad = abs(allOriDeviation{iroi}(isplit,:));
tempCard = abs(allCardDeviation{iroi}(isplit,:));
clear binSubVert binSubCard binSubRad binVert binCard binRad binSemVert binSemCard binSemRad
for ibin=1:nbins
    for isub=subjects
        goodInd = allSubInd{iroi}==isub & allPrfR2{iroi}>prfThresh & allPrfEcc{iroi}>=binBorders(ibin) & allPrfEcc{iroi}<binBorders(ibin+1);
        binSubVert(isub,ibin) = mean(tempVert(goodInd), "omitnan");
        binSubCard(isub,ibin) = mean(tempCard(goodInd), "omitnan");
        binSubRad(isub,ibin) = mean(tempRad(goodInd), "omitnan");
    end
end
binVert = mean(binSubVert, "omitnan");
binCard = mean(binSubCard, "omitnan");
binRad = mean(binSubRad, "omitnan");
binSemVert = std(binSubVert)/sqrt(length(subjects));
binSemCard = std(binSubCard)/sqrt(length(subjects));
binSemRad = std(binSubRad)/sqrt(length(subjects));

plot(binCenters,binVert,'linewidth',linewidth,'color',vertColor, 'LineStyle', lineStyleVert); hold all
dsErrorsurface(binCenters,binVert,binSemVert,vertColor,surfaceAlpha);
plot(binCenters,binCard,'linewidth',linewidth,'color',cardColor, 'LineStyle', lineStyleCard);
dsErrorsurface(binCenters,binCard,binSemCard,cardColor,surfaceAlpha);
plot(binCenters,binRad,'linewidth',linewidth,'color',radColor, 'LineStyle', lineStyleRad);
dsErrorsurface(binCenters,binRad,binSemRad,radColor,surfaceAlpha);

xlim([0 10]);
ylim([0.35 1.05]);
axis square
xlabel('\iteccentricity (deg)');
ylabel('\itdeviation (rad)');


%% POLAR plot: vertical, cardinal, and radial deviation vs. pRF angle
isubplot=isubplot+1;
subplot(rows,cols,isubplot);
% thresh = -inf;
% improv = allNsdOriCorr{iroi}(isplit,:) - allNsdCorr{iroi}(isplit,:);
% improv = improv';
tempVert = abs(allVertDeviation{iroi}(isplit,:));
tempRad = abs(allOriDeviation{iroi}(isplit,:));
tempCard = abs(allCardDeviation{iroi}(isplit,:));

clear binSubVert binSubCard binSubRad binVert binCard binRad binSemVert binSemCard binSemRad
for ibin=1:nbins
    for isub=subjects
        goodInd = allSubInd{iroi}==isub & allPrfR2{iroi}>prfThresh & allPrfAng{iroi}*pi/180>=angBinBorders(ibin) & allPrfAng{iroi}*pi/180<angBinBorders(ibin+1);
        binSubVert(isub,ibin) = mean(tempVert(goodInd), "omitnan");
        binSubCard(isub,ibin) = mean(tempCard(goodInd), "omitnan");
        binSubRad(isub,ibin) = mean(tempRad(goodInd), "omitnan");
    end
end
binVert = nanmean(binSubVert);
binCard = nanmean(binSubCard);
binRad = nanmean(binSubRad);
binSemVert = std(binSubVert)/sqrt(length(subjects));
binSemCard = std(binSubCard)/sqrt(length(subjects));
binSemRad = std(binSubRad)/sqrt(length(subjects));

polarplot([angBinCenters angBinCenters(1)],[binVert binVert(1)],'color',vertColor,'linewidth',linewidth, 'LineStyle', lineStyleVert); hold all
polarplot([angBinCenters angBinCenters(1)],[binCard binCard(1)],'color',cardColor, 'linewidth',linewidth, 'LineStyle', lineStyleCard);
polarplot([angBinCenters angBinCenters(1)],[binRad binRad(1)],'color',radColor,'linewidth',linewidth, 'LineStyle', lineStyleRad);
thetaticks(0:45:315); thetaticklabels({'0','','\pi/2','', '\pi','','3\pi/2',''});

%% line plot: per subject, orientation deviation from VERTICAL minus from RADIAL vs. pRF eccentricity
isubplot=isubplot+1;
% thresh = -inf;
subplot(rows,cols, isubplot);
% improv = allNsdOriCorr{iroi}(isplit,:) - allNsdCorr{iroi}(isplit,:);
% improv = improv';
temp = abs(allVertDeviation{iroi}(isplit,:)) - abs(allOriDeviation{iroi}(isplit,:));
xlim([0 10]);
plot([0:10],0*[0:10],'k');
hold on
clear binData
for isub=subjects
    for ibin=1:nbins
        binData(isub,ibin) = mean(temp(allSubInd{iroi}==isub & allPrfR2{iroi}>prfThresh & allPrfEcc{iroi}>=binBorders(ibin) & allPrfEcc{iroi}<binBorders(ibin+1)), "omitnan");
    end
end
plot(binCenters,binData,'color',0.7*[1 1 1]);

[~, pvalRadVert, ciRadVert, statsRadVert] = ttest(binData);
sigBins = pvalRadVert<pvalThresh;
plot(binCenters(sigBins),mean(binData(:,sigBins), "omitnan"),'LineStyle','none','color',cardColor,'linewidth',linewidth,'markersize',sigMarkersize,'marker','.');
box on
axis square
xlabel('\iteccentricity (deg)');
ylabel('\it\Delta deviation (rad)');

%% POLAR plot: per subject, orientation deviation from VERTICAL minus from RADIAL vs. pRF eccentricity

isubplot=isubplot+1;
subplot(rows,cols, isubplot);
% thresh = -inf;
% improv = allNsdOriCorr{iroi}(isplit,:) - allNsdCorr{iroi}(isplit,:);
% improv = improv';
temp = abs(allVertDeviation{iroi}(isplit,:)) - abs(allOriDeviation{iroi}(isplit,:));
for isub=subjects
    for ibin=1:nbins
        binData(isub,ibin) = mean(temp(allSubInd{iroi}==isub & allPrfR2{iroi}>prfThresh & allPrfAng{iroi}*pi/180>=angBinBorders(ibin) & allPrfAng{iroi}*pi/180<angBinBorders(ibin+1)), "omitnan");
    end
end
polarplot([angBinCenters angBinCenters(1)],[binData binData(:,1)],'color',0.7*[1 1 1]);
hold on

[~, pvalRadVert, ciRadVert, statsRadVert] = ttest(binData);
sigBins = pvalRadVert<pvalThresh;
polarplot(angBinCenters(sigBins),mean(binData(:,sigBins), "omitnan"),'LineStyle','none','color',cardColor,'linewidth',linewidth,'markersize',sigMarkersize,'marker','.');

rlim([min(binData(:))-0.1 max(binData(:))+0.1])
thetaticks(0:45:315); thetaticklabels({'0','','\pi/2','', '\pi','','3\pi/2',''});

%% line plot: per subject, orientation deviation from CARDINAL minus from RADIAL vs. pRF eccentricity
isubplot=isubplot+1;
subplot(rows,cols, isubplot);
% thresh = -inf;
% improv = allNsdOriCorr{iroi}(isplit,:) - allNsdCorr{iroi}(isplit,:);
% improv = improv';
temp = abs(allCardDeviation{iroi}(isplit,:)) - abs(allOriDeviation{iroi}(isplit,:));
xlim([0 10]);
plot([0:10],0*[0:10],'k');
hold on
clear binData
for isub=subjects
    for ibin=1:nbins
        binData(isub,ibin) = mean(temp(allSubInd{iroi}==isub & allPrfR2{iroi}>prfThresh & allPrfEcc{iroi}>=binBorders(ibin) & allPrfEcc{iroi}<binBorders(ibin+1)), "omitnan");
    end
end
plot(binCenters,binData,'color',0.7*[1 1 1]);
[~, pvalRadCard, ciRadCard, statsRadCard] = ttest(binData);
sigBins = pvalRadCard<pvalThresh;
plot(binCenters(sigBins),mean(binData(:,sigBins), "omitnan"),'LineStyle','none','color',cardColor,'linewidth',linewidth,'markersize',sigMarkersize,'marker','.');
box on
axis square
xlabel('\iteccentricity (deg)');
ylabel('\it\Delta deviation (rad)');


%% POLAR plot: per subject, orientation deviation from CARDINAL minus from RADIAL vs. pRF eccentricity
isubplot=isubplot+1;
subplot(rows,cols, isubplot);
% thresh = -inf;
% improv = allNsdOriCorr{iroi}(isplit,:) - allNsdCorr{iroi}(isplit,:);
% improv = improv';
temp = abs(allCardDeviation{iroi}(isplit,:)) - abs(allOriDeviation{iroi}(isplit,:));
for isub=subjects
    for ibin=1:nbins
        binData(isub,ibin) = mean(temp(allSubInd{iroi}==isub & allPrfR2{iroi}>prfThresh & allPrfAng{iroi}*pi/180>=angBinBorders(ibin) & allPrfAng{iroi}*pi/180<angBinBorders(ibin+1)), "omitnan");
    end
end
polarplot([angBinCenters angBinCenters(1)],[binData binData(:,1)],'color',0.7*[1 1 1]);
hold on

[~, pvalRadCard, ciRadCard, statsRadCard] = ttest(binData);
sigBins = pvalRadCard<pvalThresh;
polarplot(angBinCenters(sigBins),mean(binData(:,sigBins), "omitnan"),'LineStyle','none','color',cardColor,'linewidth',linewidth,'markersize',sigMarkersize,'marker','.');

rlim([min(binData(:))-0.1 max(binData(:))+0.1])
thetaticks(0:45:315); thetaticklabels({'0','','\pi/2','', '\pi','','3\pi/2',''});

%% line plot: orientation deviation from RADIAL vs. angular distance from MERIDIAN
% isubplot=isubplot+1;
% subplot(rows,cols, isubplot);
% thresh = -inf;
% prfThresh = 00;
% numMeridians = 4;
% 
% improv = allNsdOriCorr{iroi}(isplit,:) - allNsdCorr{iroi}(isplit,:);
% improv = improv';
% temp = mod(allPrfAng{iroi},(360/numMeridians));
% temp = temp*2*pi/(360/numMeridians);
% meridDist = circ_dist(temp,0)/numMeridians;
% meridDist = abs(meridDist);
% 
% [corrMeridRadDevi(iroi,isplit), pCorrMeridRadDevi(iroi,isplit)] = corr(meridDist(allPrfR2{iroi}>prfThresh),abs(allOriDeviation{iroi}(isplit,allPrfR2{iroi}>prfThresh))');
% goodInd = allPrfR2{iroi}>prfThresh & allPrfAng{iroi}<180;
% [corrMeridRadDeviUpper(iroi,isplit), pCorrMeridRadDeviUpper(iroi,isplit)] = corr(meridDist(goodInd),abs(allOriDeviation{iroi}(isplit,goodInd))');
% goodInd = allPrfR2{iroi}>prfThresh & allPrfAng{iroi}>=180;
% [corrMeridRadDeviLower(iroi,isplit), pCorrMeridRadDeviLower(iroi,isplit)] = corr(meridDist(goodInd),abs(allOriDeviation{iroi}(isplit,goodInd))');
% [corrMeridPrfR2(iroi), pCorrMeridPrfR2(iroi)] = corr(meridDist(allPrfR2{iroi}>prfThresh),allPrfR2{iroi}(allPrfR2{iroi}>prfThresh));
% 
% [partCorrMeridRad(iroi,isplit), pPartCorrMeridRad(iroi,isplit)] = partialcorr(meridDist(allPrfR2{iroi}>prfThresh),abs(allOriDeviation{iroi}(isplit,allPrfR2{iroi}>prfThresh))',allPrfR2{iroi}(allPrfR2{iroi}>prfThresh));
% 
% clear binR2 binData binSubR2 binSubData
% for ibin=1:nbins
%     for isub=subjects
%         binSubR2(isub,ibin) = mean(allPrfR2{iroi}(allSubInd{iroi}==isub & allPrfR2{iroi}>prfThresh &  meridDist>=angBinBorders(ibin)/(numMeridians*2) & meridDist<angBinBorders(ibin+1)/(numMeridians*2))/100);
%         binSubData(isub,ibin) = mean(abs(allOriDeviation{iroi}(isplit,allSubInd{iroi}==isub & allPrfR2{iroi}>prfThresh &  meridDist>=angBinBorders(ibin)/(numMeridians*2) & meridDist<angBinBorders(ibin+1)/(numMeridians*2))));
%     end
% end
% binR2 = mean(binSubR2);
% binData = mean(binSubData);
% binSemR2 = std(binSubR2)/sqrt(length(subjects));
% binSemData = std(binSubData)/sqrt(length(subjects));
% axis square
% xlabel('\itdist from meridian (rad)');
% box on
% hold on
% p=plot(angBinCenters/(numMeridians*2),binData,'color',radColor,'linewidth',linewidth);
% 
% dsErrorsurface(angBinCenters/(numMeridians*2),binData,binSemData,radColor,surfaceAlpha);
% ylabel('\itradial deviation (rad)','color','k');
% 
% xlim([0 pi/4]);

%% line plot: pRF R^2 vs. angular distance from MERIDIAN

% isubplot=isubplot+1;
% subplot(rows,cols, isubplot);
% 
% improv = allNsdOriCorr{iroi}(end,:) - allNsdCorr{iroi}(end,:);
% improv = improv';
% temp = mod(allPrfAng{iroi},(360/numMeridians));
% temp = temp*2*pi/(360/numMeridians);
% meridDist = circ_dist(temp,0)/numMeridians;
% meridDist = abs(meridDist);
% 
% clear binR2 binData binSubR2 binSubData
% for ibin=1:nbins
%     for isub=subjects
%         binSubR2(isub,ibin) = mean(allPrfR2{iroi}(allSubInd{iroi}==isub & allPrfR2{iroi}>prfThresh &  meridDist>=angBinBorders(ibin)/(numMeridians*2) & meridDist<angBinBorders(ibin+1)/(numMeridians*2))/100);
%     end
% end
% binR2 = mean(binSubR2);
% binSemR2 = std(binSubR2)/sqrt(length(subjects));
% plot(angBinCenters/(numMeridians*2),binR2,'color',r2Color,'linewidth',linewidth);
% dsErrorsurface(angBinCenters/(numMeridians*2),binR2,binSemR2,r2Color,surfaceAlpha*2);
% 
% axis square
% xlabel('\itdist from meridian (rad)');
% ylabel('\itpRF R^2');
% box on
% hold on

%%
set(gcf,'position',[150 180 420 740]);

if condition == 1
    figName = 'quantMap_maxmin_old';
elseif condition == 2
    figName = 'quantMap_maxmin_ori';
elseif condition == 3
    figName = 'quantMap_maxmin_control';
end

if toSavePdf
    print('-painters','-dpdf',[figFolder imgTypes{curimgtype} '_' pairTypes{curpairtype} '_' figName]);
end


%% MEAN DEVIATION DIRECTION VS PRF ANGLE
% ifig=ifig+1; h=figure(ifig); clf;
% 
% tempVert = allVertDeviation{iroi}(isplit,:);
% tempRad = allOriDeviation{iroi}(isplit,:);
% tempCard = allCardDeviation{iroi}(isplit,:);
% 
% clear binSubVert binSubCard binSubRad binVert binCard binRad binSemVert binSemCard binSemRad
% for ibin=1:nbins
%     for isub=subjects
%         goodInd = allSubInd{iroi}==isub & allPrfR2{iroi}>prfThresh & improv>thresh &  allPrfAng{iroi}*pi/180>=angBinBorders(ibin) & allPrfAng{iroi}*pi/180<angBinBorders(ibin+1);
%         binSubVert(isub,ibin) = circ_mean(tempVert(goodInd)');
%         binSubCard(isub,ibin) = circ_mean(tempCard(goodInd)');
%         binSubRad(isub,ibin) = circ_mean(tempRad(goodInd)');
%     end
% end
% binVert = nanmean(binSubVert);
% binCard = nanmean(binSubCard);
% binRad = nanmean(binSubRad);
% binSemVert = std(binSubVert)/sqrt(length(subjects));
% binSemCard = std(binSubCard)/sqrt(length(subjects));
% binSemRad = std(binSubRad)/sqrt(length(subjects));
% plot(angBinCenters,binRad,'color',radColor,'linewidth',linewidth);
% dsErrorsurface(angBinCenters,binRad,binSemRad,radColor,surfaceAlpha);
% 

%%
function prefOri = plotOriLines(prefOri, prfX, prfY,r2,lineColor)
r2 = (r2.^2)*600;
minWidth = 0.001;
r2(isnan(r2)) = minWidth;
r2(r2<minWidth) = minWidth;
numvox = length(prefOri);

lineWidth = 0.01*r2;
lineLength = 0.4*r2;

for ivox=1:numvox
    %if the coefficients are NaN, don't plot
    if ~isnan(prefOri(ivox))
        h=drawOriLine(prfX(ivox), prfY(ivox), pi/2-prefOri(ivox), lineLength(ivox), lineWidth(ivox), lineColor);
        hold on
    end
end

global interpSz;% = 714;
global backgroundSz;% = 1024;
global degPerPix;
set(gca,'XTick',[]);
set(gca,'YTick',[]);
linecolor = [0.5 0.5 0.5];
line([0 0], [-backgroundSz backgroundSz],'color',linecolor);
line([-backgroundSz backgroundSz],[0 0], 'color',linecolor);
line([-interpSz/2 -interpSz/2], [-interpSz/2 interpSz/2], 'color',linecolor);
line([interpSz/2 interpSz/2], [-interpSz/2 interpSz/2], 'color',linecolor);
line([-interpSz/2 interpSz/2],[interpSz/2 interpSz/2], 'color',linecolor);
line([-interpSz/2 interpSz/2],[-interpSz/2 -interpSz/2], 'color',linecolor);
xlim([-interpSz interpSz]); ylim([-interpSz interpSz]);
box on
axis square

end
