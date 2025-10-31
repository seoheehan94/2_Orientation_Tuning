% fig5.m
%
% associated with the following publication: Roth, ZN, Kay, K, and Merriam, EP (2022).
% Massive natural scene sampling reveals reliable coarse-scale orientation tuning in human V1
% DOI:
%
%   usage: fig5()
%   by: zvi roth
%   date: 7/29/2022
%   purpose: plot preferred orientations for Figure 5
%   uses files created by: getVoxPref.m


close all
clear all
tic

con = 3;

if con == 1
    condition = 'ldSPresidual_photoSP';
elseif con == 2
    condition = 'photoSPresidual_ldSP';
elseif con == 3
    condition = 'contourresidual_ldSP';
end




toSavePdf = 1;
imgFormat = 'jpg';
subjects = [1:8];

ifig=0;
nrois = 4;

imgScaling = 0.5;
global interpSz; interpSz= 714*imgScaling;
global backgroundSz; backgroundSz= 1024*imgScaling;
global degPerPix; degPerPix = 8.4/interpSz;

%scatter parameters
markersize = 1;
edgeAlpha = 0.3;%0.07
markerColor = [0 0 0];
prfThresh = 0;

prffolder = ['/bwdata/NSDData/Seohee/Orientation/'];
figFolder = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/figures/'];

allOri = cell(1,nrois);
allResidOri = cell(1,nrois);
allPrfR2 = cell(1,nrois);
allPrfX = cell(1,nrois);
allPrfY = cell(1,nrois);
allPrfEcc = cell(1,nrois);
allPrfAng = cell(1,nrois);
allNsdOriR2 = cell(1,nrois);

for isub=1:length(subjects)
    subnum = subjects(isub);
    load([prffolder condition '_sub' num2str(isub) '.mat'],'allRoiPrf',...
        'residOri','roiNsdOriR2', 'nsplits');

    % subNsdSynthImprov_corr(isub,:,:) = nsdSynthImprov_corr;
    % subNsdSynthImprov_pval(isub,:,:) = nsdSynthImprov_pval;
    for iroi=1:nrois
        allPrfX{iroi} = [allPrfX{iroi}; allRoiPrf{iroi}.x];
        allPrfY{iroi} = [allPrfY{iroi}; allRoiPrf{iroi}.y];
        allPrfEcc{iroi} = [allPrfEcc{iroi}; allRoiPrf{iroi}.ecc];
        allPrfAng{iroi} = [allPrfAng{iroi}; allRoiPrf{iroi}.ang];
       allPrfR2{iroi} = [allPrfR2{iroi}; allRoiPrf{iroi}.r2];
       allResidOri{iroi} =  [allResidOri{iroi} residOri{iroi}];
       allNsdOriR2{iroi} = [allNsdOriR2{iroi} roiNsdOriR2{iroi}];
      
    end
end


% ifig=ifig+1; h=figure(ifig); clf;
% rows=1;
% cols=3;
isplit = nsplits;
% isubplot=0;

iroi=1;

%% RESIDUAL preferred ORIENTATION
% isubplot=isubplot+1;
% subplot(rows,cols, isubplot);
%plotOriLines(allResidOri{iroi}(isplit,:), allPrfX{iroi}, allPrfY{iroi}, allPrfEcc{iroi},(0.4*allNsdOriR2{iroi}(isplit,:)));
fig=figure;
plotOriLines(allResidOri{iroi}(isplit,:), allPrfX{iroi}, allPrfY{iroi}, allPrfEcc{iroi},(3*allNsdOriR2{iroi}(isplit,:)));

xlabel('\itx position (deg)');
ylabel('\ity position (deg)');

exportgraphics(fig, ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/figures/' condition, '.png'])
exportgraphics(fig, ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/figures/' condition, '.pdf'])

toc
%%
function prefOri = plotOriLines(prefOri, prfX, prfY, prfEcc,r2)
r2 = (r2)*200;
minWidth = 0.001;
r2(isnan(r2)) = minWidth;
r2(r2<minWidth) = minWidth;
numvox = length(prefOri);

lineWidth = 0.01*r2;
lineLength = 0.4*r2;

cMap = turbo(256);

for ivox=1:numvox
    %if the coefficients are NaN, don't plot
    if ~isnan(prefOri(ivox))
        h=drawOriLine(prfX(ivox), prfY(ivox), pi/2-prefOri(ivox), lineLength(ivox), lineWidth(ivox), cMap(1+floor((prefOri(ivox))*255/(pi)),:));
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
set(gca,'xTick',[-interpSz/2 0 interpSz/2]);
set(gca,'xTicklabels',{-degPerPix*interpSz/2, 0, degPerPix*interpSz/2});
set(gca,'yTick',[-interpSz/2 0 interpSz/2]);
set(gca,'yTicklabels',{-degPerPix*interpSz/2, 0, degPerPix*interpSz/2});
box on
axis square

end
