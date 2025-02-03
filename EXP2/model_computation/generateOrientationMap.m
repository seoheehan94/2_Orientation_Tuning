function oriMap = generateOrientationMap(vecLD,backgroundValue, backgroundSize,renderSize)

if nargin < 4
    renderSize = vecLD.imsize;
end
if nargin < 3
    backgroundSize = renderSize;
end
if nargin < 2
    backgroundValue = NaN;
end
oriMap = zeros(backgroundSize) + backgroundValue;

scaleY = renderSize(1)/vecLD.imsize(1);
scaleX = renderSize(2)/vecLD.imsize(2);
offsetY = round((backgroundSize(1) - renderSize(1))/2);
offsetX = round((backgroundSize(2) - renderSize(2))/2);

for c = 1:vecLD.numContours
     oris = mod(vecLD.orientations{c},180);
     % oris = mod((180-vecLD.orientations{c}),180);
     

    for s = 1:size(vecLD.contours{c},1)
        thisMap = zeros(backgroundSize(1),backgroundSize(2),3);
        newcord = vecLD.contours{c}(s,:);
        newcord(1,[1,3]) = newcord(1,[1,3])*scaleX + offsetX;
        newcord(1,[2,4]) = newcord(1,[2,4])*scaleY + offsetY;
        thisMap = insertShape(thisMap,'Line',newcord,'Color',[1,0,0],'LineWidth',1,'Opacity',1,'SmoothEdges',false);
        thisMap = thisMap(:,:,1);
        thisIdx = (thisMap > 0);
        oriMap(thisIdx) = oris(s);
    end
end
    
% %vecLD.orientationBins
% for binIdx = 1: length(vecLD.orientationBins)
%     oribinmap{1,binIdx} = (abs(oriMap-vecLD.orientationBins(binIdx)) < 11.25);
% 
% end
% modelOri=cat(3,oribinmap{:});
% modelOri = permute(modelOri,[3 1 2]);
% 
