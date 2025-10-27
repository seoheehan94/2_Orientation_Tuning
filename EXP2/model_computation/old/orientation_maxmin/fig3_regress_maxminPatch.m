% fig3.m
%
% associated with the following publication: Roth, ZN, Kay, K, and Merriam, EP (2022).
% Massive natural scene sampling reveals reliable coarse-scale orientation tuning in human V1
% DOI:
%
%   usage: fig3()
%   by: zvi roth
%   date: 7/29/2022
%   purpose: Plot preferred orientations for Figure 3
%   uses files created by: getVoxPref.m
% clear all
tic

condition =2; %1=old, 2=ori, 3=control
curpairtype =2;
imgTypes = {'Top', 'Bottom'};
pairTypes = {'old_control', 'old_ori', 'control_ori'};

toSavePdf = 1;

imgFormat = 'jpg';
subjects = [1:8];
% subjects = [7];
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

if condition == 1
    prffolder = ['/bwdata/NSDData/Seohee/Orientation/prfsample_maxmin/'];
elseif condition == 2
    prffolder = ['/bwdata/NSDData/Seohee/Orientation/prfsample_maxmin_Ori/'];
elseif condition == 3
    prffolder = ['/bwdata/NSDData/Seohee/Orientation/prfsample_maxmin_control/'];
end
figFolder = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/figures/'];

allOri = cell(length(imgTypes), nrois);
allPrfR2 = cell(length(imgTypes), nrois);
allPrfX = cell(length(imgTypes), nrois);
allPrfY = cell(length(imgTypes), nrois);
allPrfEcc = cell(length(imgTypes), nrois);
allPrfAng = cell(length(imgTypes), nrois);
allNsdOriR2 = cell(length(imgTypes), nrois);

for curimgtype = 1:length(imgTypes)
    for isub = 1:length(subjects)
        subnum = subjects(isub);
        % Load data for the current imgType
        dataFile = [prffolder 'indices' imgTypes{curimgtype} '_' pairTypes{curpairtype} 'voxModelPref_regress_sub' num2str(isub) '.mat'];
        load(dataFile);
        
        % Process each ROI for the current imgType
        for iroi = 1:nrois
            allPrfX{curimgtype, iroi} = [allPrfX{curimgtype, iroi}; allRoiPrf{iroi}.x];
            allPrfY{curimgtype, iroi} = [allPrfY{curimgtype, iroi}; allRoiPrf{iroi}.y];
            allPrfEcc{curimgtype, iroi} = [allPrfEcc{curimgtype, iroi}; allRoiPrf{iroi}.ecc];
            allPrfAng{curimgtype, iroi} = [allPrfAng{curimgtype, iroi}; allRoiPrf{iroi}.ang];
            allPrfR2{curimgtype, iroi} = [allPrfR2{curimgtype, iroi}; allRoiPrf{iroi}.r2];
            allOri{curimgtype, iroi} = [allOri{curimgtype, iroi} roiOri{iroi}];
            allNsdOriR2{curimgtype, iroi} = [allNsdOriR2{curimgtype, iroi} roiNsdOriR2{iroi}];
        end
    end
end

%%
figure;
ifig=1; h=figure(ifig); clf;
rows=1;
cols=2;
isplit = 3;

iroi=1;

%% preferred ORIENTATION

for curimgtype = 1:2
    subplot(1, 2, curimgtype);
    set(gca, 'FontName', 'Helvetica', 'FontSize', 18, 'FontAngle', 'normal');
    plotOriLines(allOri{curimgtype, iroi}(isplit, :), ...
                 allPrfX{curimgtype, iroi}, ...
                 allPrfY{curimgtype, iroi}, ...
                 allPrfEcc{curimgtype, iroi}, ...
                 (3 * allNsdOriR2{curimgtype, iroi}(isplit, :)));
end

% set(gcf,'position',[150 180 3*250 rows*210]);
h.Units = 'centimeters';
h.PaperSize=[10 5];
if toSavePdf
    print('-painters','-dpdf',[figFolder 'radialBias_maxminPatch_old_ori_ori']);
end
% 
% toc
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
set(gca,'xTick',[-interpSz/2 0 interpSz/2],  'FontSize', 18);
set(gca,'xTicklabels',{-degPerPix*interpSz/2, 0, degPerPix*interpSz/2});
set(gca,'yTick',[-interpSz/2 0 interpSz/2],  'FontSize', 18);
set(gca,'yTicklabels',{-degPerPix*interpSz/2, 0, degPerPix*interpSz/2});
box on
axis square
set(gca, 'LineWidth', 1.5);

end
