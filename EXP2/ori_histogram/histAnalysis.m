%%
clear all;
rng(4228);

%% Load data
original_folder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/ori_hist';
cd(original_folder);

hist_vecLD = load("vecLD/vecLD.mat");
hist_photo = load("photo/photo.mat");
hist_LD = load("LD/LD.mat");

% quad_fields = {'LeftTop';'RightTop';'LeftBottom';'RightBottom'};
sq_fields = {'sq1';'sq2';'sq3';'sq4';'sq5';'sq6';'sq7'};
ori_fields = {'ori1';'ori2';'ori3';'ori4';'ori5';'ori6';'ori7';'ori8'};

imageSize = [512 512];
% totalSum = imageSize(1)*imageSize(2); 

% normalize
% for sqIdx = 1:7
%     for oriIdx = 1:8
%         hist_photo.modelOriSum.(sq_fields{sqIdx}).(ori_fields{oriIdx}) = hist_photo.modelOriSum.(sq_fields{sqIdx}).(ori_fields{oriIdx})/totalSum;
%         hist_LD.modelOriSum.(sq_fields{sqIdx}).(ori_fields{oriIdx}) = hist_LD.modelOriSum.(sq_fields{sqIdx}).(ori_fields{oriIdx})/totalSum;
%     end
% end

for oriIdx = 1:8
    hist_photo_allsq.modelOriSum.(ori_fields{oriIdx}) = hist_photo.modelOriSum.sq1.(ori_fields{oriIdx}) +...
        hist_photo.modelOriSum.sq2.(ori_fields{oriIdx}) + ...
        hist_photo.modelOriSum.sq3.(ori_fields{oriIdx}) + ...
        hist_photo.modelOriSum.sq4.(ori_fields{oriIdx}) + ...
        hist_photo.modelOriSum.sq5.(ori_fields{oriIdx}) + ...
        hist_photo.modelOriSum.sq6.(ori_fields{oriIdx}) + ...
        hist_photo.modelOriSum.sq7.(ori_fields{oriIdx});
end


for oriIdx = 1:8
    hist_LD_allsq.modelOriSum.(ori_fields{oriIdx}) = hist_LD.modelOriSum.sq1.(ori_fields{oriIdx}) +...
        hist_LD.modelOriSum.sq2.(ori_fields{oriIdx}) + ...
        hist_LD.modelOriSum.sq3.(ori_fields{oriIdx}) + ...
        hist_LD.modelOriSum.sq4.(ori_fields{oriIdx}) + ...
        hist_LD.modelOriSum.sq5.(ori_fields{oriIdx}) + ...
        hist_LD.modelOriSum.sq6.(ori_fields{oriIdx}) + ...
        hist_LD.modelOriSum.sq7.(ori_fields{oriIdx});
end

%% get sum for each orientation bin
%vecLD
for oriIdx = 1:8
    meanOri_vecLD(oriIdx) = mean(hist_vecLD.modelOriSum.(ori_fields{oriIdx}));
end
% figure;
% plot(meanOri_vecLD)

%LD
for oriIdx = 1:8
    meanOri_LD(oriIdx) = mean(hist_LD_allsq.modelOriSum.(ori_fields{oriIdx}));
end

for sqIdx = 1:7
for oriIdx = 1:8
    sumOri_LD.(sq_fields{sqIdx})(oriIdx) = sum(hist_LD.modelOriSum.(sq_fields{sqIdx}).(ori_fields{oriIdx}));
end
end

for sqIdx = 1:7
    thissqSum = sum(sumOri_LD.(sq_fields{sqIdx}));
    for oriIdx = 1:8
        normMeanOri_LD.(sq_fields{sqIdx})(oriIdx) = (sumOri_LD.(sq_fields{sqIdx})(oriIdx))/thissqSum;
    end
end

% for sqIdx = 1:7
%     for oriIdx = 1:8
%         meanOri_LD.(sq_fields{sqIdx})(oriIdx) = mean(hist_LD.modelOriSum.(sq_fields{sqIdx}).(ori_fields{oriIdx}));
%     end
% end
% 
% for sqIdx = 1:7
%     thissqSum = sum(meanOri_LD.(sq_fields{sqIdx}));
%     for oriIdx = 1:8
%         normMeanOri_LD.(sq_fields{sqIdx})(oriIdx) = (meanOri_LD.(sq_fields{sqIdx})(oriIdx))/thissqSum;
%     end
% end
% 
% figure;
% for sqIdx = 1:7
%     hold on;
% plot(meanOri_LD.(sq_fields{sqIdx}));
% end
% legend('1','2','3','4','5','6','7','FontSize',15,'FontName','Helvetica');
% title('all sq LD');
% 
figure;
for sqIdx = 1:7
    hold on;
plot(normMeanOri_LD.(sq_fields{sqIdx}));
end
legend('1','2','3','4','5','6','7','FontSize',15,'FontName','Helvetica');
title('proportion all sq LD');

% 
% sumLD = struct;
% for sfIdx = 1:7
%     for oriIdx = 1:8
%         sumLD.(sq_fields{sfIdx})(oriIdx) = sum(hist_LD.modelOriSum.(sq_fields{sfIdx}).(ori_fields{oriIdx}));
%     end
% end
% 
% totalSumLD = struct;
% for sfIdx = 1:7
%     totalSumLD.(sq_fields{sfIdx}) = sum(sumLD.(sq_fields{sfIdx}), 'All');
% end
% 
% normSumLD = struct;
% for sfIdx = 1:7
%     for oriIdx = 1:8
%         normSumLD.(sq_fields{sfIdx})(oriIdx) = sumLD.(sq_fields{sfIdx})(oriIdx)/totalSumLD.(sq_fields{sfIdx});
%     end
% end
% 
% for oriIdx = 1:8
%     sumLD_allsq(oriIdx) = sum(hist_LD_allsq.modelOriSum.(ori_fields{oriIdx}));
% end
% 
% totalSumLD_allsq = sum(sumLD_allsq, 'All');
% 
% for oriIdx = 1:8
%     normSumLD_allsq(oriIdx) = sumLD_allsq(oriIdx)/totalSumLD_allsq;
% end

%Photo
for oriIdx = 1:8
    meanOri_photo(oriIdx) = mean(hist_photo_allsq.modelOriSum.(ori_fields{oriIdx}));
end

for sqIdx = 1:7
for oriIdx = 1:8
    sumOri_photo.(sq_fields{sqIdx})(oriIdx) = sum(hist_photo.modelOriSum.(sq_fields{sqIdx}).(ori_fields{oriIdx}));
end
end

for sqIdx = 1:7
    thissqSum = sum(sumOri_photo.(sq_fields{sqIdx}));
    for oriIdx = 1:8
        normMeanOri_photo.(sq_fields{sqIdx})(oriIdx) = (sumOri_photo.(sq_fields{sqIdx})(oriIdx))/thissqSum;
    end
end
% 
% 
% figure;
% plot(meanOri_vecLD)
% for sqIdx = 1:7
%     hold on;
% plot(normMeanOri_photo.(sq_fields{sqIdx}));
% end
% legend('vecLD','1','2','3','4','5','6','7','FontSize',15,'FontName','Helvetica');
% title('proportion: vecLD vs. all sq Photo');
% 
% figure;
% for sqIdx = 1:7
%     hold on;
% plot(meanOri_photo.(sq_fields{sqIdx}));
% end
% legend('1','2','3','4','5','6','7','FontSize',15,'FontName','Helvetica');
% title('all sq Photo');
% 
figure;
for sqIdx = 1:7
    hold on;
plot(normMeanOri_photo.(sq_fields{sqIdx}));
end
legend('1','2','3','4','5','6','7','FontSize',15,'FontName','Helvetica');
title('proportion all sq Photo');


%% Allsq photo vs vecLD vs LD
figure;
plot(meanOri_photo,'-o', 'Color','#0072BD','LineWidth',2);
hold on;
plot(meanOri_vecLD,'-o', 'Color','#F35872','LineWidth',2);
hold on;
plot(meanOri_LD,'-o', 'LineWidth',2, 'Color', '#65B74A');

ylim([0 0.3])
set(gca,'xticklabel',[])
legend('photograph-filter','contour','line drawing-filter','FontSize',15,'FontName','Helvetica');
ax = gca;
ax.YAxis.FontSize = 15;
ax.YAxis.FontName = 'Helvetica';
title('Proportion of Orientations', 'FontSize',20,'FontName','Helvetica');
box off;
legend boxoff

% saveas(gcf,['ori_hist' '.pdf']);

%% sq1 photo vs vecLD vs sq1 LD

figure;
plot(normMeanOri_photo.sq1,'-o', 'Color','#0072BD','LineWidth',2);
hold on;
plot(meanOri_vecLD,'-o', 'Color','#F35872','LineWidth',2);
hold on;
plot(normMeanOri_LD.sq1,'-o', 'LineWidth',2, 'Color', '#65B74A');

ylim([0 0.3])
set(gca,'xticklabel',[])
legend('photograph-filter1','contour','line drawing-filter1','FontSize',15,'FontName','Helvetica');
ax = gca;
ax.YAxis.FontSize = 15;
ax.YAxis.FontName = 'Helvetica';
title('Proportion of Orientations', 'FontSize',20,'FontName','Helvetica');
box off;
legend boxoff

%%
figure;
plot(normMeanOri_photo.sq1,'-o', 'Color','#0072BD','LineWidth',2);
hold on;
plot(normMeanOri_LD.sq1,'-o', 'LineWidth',2, 'Color', '#65B74A');

ylim([0 0.3])
set(gca,'xticklabel',[])
legend('photograph-filter1','contour','line drawing-filter1','FontSize',15,'FontName','Helvetica');
ax = gca;
ax.YAxis.FontSize = 15;
ax.YAxis.FontName = 'Helvetica';
title('Proportion of Orientations', 'FontSize',20,'FontName','Helvetica');
box off;
legend boxoff
%% Photo vs LD
% for sfIdx = 1:7
% 
%         figure;
%         subplot(2,2,1)
%         bar(normSumPhoto.LeftTop.(sq_fields{sfIdx}),'facealpha',.7,'edgecolor','none');
%         hold on;
%         bar(normSumLD.LeftTop.(sq_fields{sfIdx}),'facealpha',.3,'FaceColor','red', 'edgecolor','none');
%         legend('photo','LD');
%         title('LeftTop');
% 
% 
%         subplot(2,2,2)
%         bar(normSumPhoto.RightTop.(sq_fields{sfIdx}),'facealpha',.7,'edgecolor','none');
%         hold on;
%         bar(normSumLD.RightTop.(sq_fields{sfIdx}),'facealpha',.3,'FaceColor','red','edgecolor','none');
%         legend('photo','LD');
%         title('RightTop');
% 
%         subplot(2,2,3)
%         bar(normSumPhoto.LeftBottom.(sq_fields{sfIdx}),'facealpha',.7,'edgecolor','none');
%         hold on;
%         bar(normSumLD.LeftBottom.(sq_fields{sfIdx}),'facealpha',.3,'FaceColor','red','edgecolor','none');
%         legend('photo','LD');
%         title('LeftBottom');
% 
%         subplot(2,2,4)
%         bar(normSumPhoto.RightBottom.(sq_fields{sfIdx}),'facealpha',.7,'edgecolor','none');
%         hold on;
%         bar(normSumLD.RightBottom.(sq_fields{sfIdx}),'facealpha',.3,'FaceColor','red','edgecolor','none');
%         legend('photo','LD');
%         title('RightBottom');
% 
%         sgtitle(sq_fields{sfIdx});
% 
%         saveas(gcf,['photoLD_' (sq_fields{sfIdx}) '.png']);
% 
% 
% end
% 
%   figure;
%     subplot(2,2,1)
%     bar(normSumLD_allsq.LeftTop,'facealpha',.7,'edgecolor','none');
%     hold on;
%     bar(normSumPhoto_allsq.LeftTop,'facealpha',.3,'edgecolor','none','FaceColor','red');
%     legend('LD','Photo');
%     title('LeftTop');
% 
%     subplot(2,2,2)
%     bar(normSumLD_allsq.RightTop,'facealpha',.7,'edgecolor','none');
%     hold on;
%     bar(normSumPhoto_allsq.RightTop,'facealpha',.3,'edgecolor','none','FaceColor','red');
%     legend('LD','Photo');
%     title('RightTop');
% 
%     subplot(2,2,3)
%     bar(normSumLD_allsq.LeftBottom,'facealpha',.7,'edgecolor','none');
%     hold on;
%     bar(normSumPhoto_allsq.LeftBottom,'facealpha',.3,'edgecolor','none','FaceColor','red');
%     legend('LD','Photo');
%     title('LeftBottom');
% 
%     subplot(2,2,4)
%     bar(normSumLD_allsq.RightBottom,'facealpha',.7,'edgecolor','none');
%     hold on;
%     bar(normSumPhoto_allsq.RightBottom,'facealpha',.3,'edgecolor','none','FaceColor','red');
%     legend('LD','Photo');
%     title('RightBottom');
% 
%     sgtitle('LDPhoto');
% 
%     saveas(gcf,['LDPhoto_allsq' '.png']);
% 
% 
% %% Photo vs vecLD
% 
%     figure;
%     subplot(2,2,1)
%     bar(normSumPhoto_allsq.LeftTop,'facealpha',.7,'edgecolor','none');
%     hold on;
%     bar(normSumVecLD.LeftTop,'facealpha',.3,'edgecolor','none','FaceColor','red');
%     legend('photo','vecLD');
%     title('LeftTop');
% 
%     subplot(2,2,2)
%     bar(normSumPhoto_allsq.RightTop,'facealpha',.7,'edgecolor','none');
%     hold on;
%     bar(normSumVecLD.RightTop,'facealpha',.3,'edgecolor','none','FaceColor','red');
%     legend('photo','vecLD');
%     title('RightTop');
% 
%     subplot(2,2,3)
%     bar(normSumPhoto_allsq.LeftBottom,'facealpha',.7,'edgecolor','none');
%     hold on;
%     bar(normSumVecLD.LeftBottom,'facealpha',.3,'edgecolor','none','FaceColor','red');
%     legend('photo','vecLD');
%     title('LeftBottom');
% 
%     subplot(2,2,4)
%     bar(normSumPhoto_allsq.RightBottom,'facealpha',.7,'edgecolor','none');
%     hold on;
%     bar(normSumVecLD.RightBottom,'facealpha',.3,'edgecolor','none','FaceColor','red');
%     legend('photo','vecLD');
%     title('RightBottom');
% 
%     sgtitle('photovecLD');
% 
%     saveas(gcf,['photovecLD_allsq' '.png']);
% 
% 
% 
% 
% %% LD vs vecLD
% 
%     figure;
%     subplot(2,2,1)
%     bar(normSumLD_allsq.LeftTop,'facealpha',.7,'edgecolor','none');
%     hold on;
%     bar(normSumVecLD.LeftTop,'facealpha',.3,'edgecolor','none','FaceColor','red');
%     legend('LD','vecLD');
%     title('LeftTop');
% 
%     subplot(2,2,2)
%     bar(normSumLD_allsq.RightTop,'facealpha',.7,'edgecolor','none');
%     hold on;
%     bar(normSumVecLD.RightTop,'facealpha',.3,'edgecolor','none','FaceColor','red');
%     legend('LD','vecLD');
%     title('RightTop');
% 
%     subplot(2,2,3)
%     bar(normSumLD_allsq.LeftBottom,'facealpha',.7,'edgecolor','none');
%     hold on;
%     bar(normSumVecLD.LeftBottom,'facealpha',.3,'edgecolor','none','FaceColor','red');
%     legend('LD','vecLD');
%     title('LeftBottom');
% 
%     subplot(2,2,4)
%     bar(normSumLD_allsq.RightBottom,'facealpha',.7,'edgecolor','none');
%     hold on;
%     bar(normSumVecLD.RightBottom,'facealpha',.3,'edgecolor','none','FaceColor','red');
%     legend('LD','vecLD');
%     title('RightBottom');
% 
%     sgtitle('LDvecLD');
% 
%     saveas(gcf,['LDvecLD_allsq' '.png']);
% 
% 
%  %% LD each Freq vs vecLD
% 
% for sfIdx = 1:7
% 
%         figure;
%         subplot(2,2,1)
%         bar(normSumVecLD.LeftTop,'facealpha',.7,'edgecolor','none');
%         hold on;
%         bar(normSumLD.LeftTop.(sq_fields{sfIdx}),'facealpha',.3,'FaceColor','red', 'edgecolor','none');
%         legend('vecLD','LD');
%         title('LeftTop');
% 
% 
%         subplot(2,2,2)
%         bar(normSumVecLD.RightTop,'facealpha',.7,'edgecolor','none');
%         hold on;
%         bar(normSumLD.RightTop.(sq_fields{sfIdx}),'facealpha',.3,'FaceColor','red', 'edgecolor','none');
%         legend('vecLD','LD');
%         title('RightTop');
% 
%         subplot(2,2,3)
%         bar(normSumVecLD.LeftBottom,'facealpha',.7,'edgecolor','none');
%         hold on;
%         bar(normSumLD.LeftBottom.(sq_fields{sfIdx}),'facealpha',.3,'FaceColor','red', 'edgecolor','none');
%         legend('vecLD','LD');
%         title('LeftBottom');
% 
%         subplot(2,2,4)
%         bar(normSumVecLD.RightBottom,'facealpha',.7,'edgecolor','none');
%         hold on;
%         bar(normSumLD.RightBottom.(sq_fields{sfIdx}),'facealpha',.3,'FaceColor','red', 'edgecolor','none');
%         legend('vecLD','LD');
%         title('RightBottom');
% 
%         sgtitle(sq_fields{sfIdx});
% 
%         saveas(gcf,['LDvecLD_' (sq_fields{sfIdx}) '.png']);
% 
% 
% end