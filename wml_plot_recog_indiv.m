clear all; close all; clc
%%%%%%%%%%% NOTE YET WORKING %%%%%%%%%%%%

% Set working directories.
rootDir = '/Volumes/Seagate/wml/';
remove = []; %[27 32 35 40];
testortestgen = 'test'; %'test' or 'test_gen'

% Create date-specific file name that indicates how many subjects.
datestring = '20211015';
filename = sprintf('wml_beh_data_recog_%s_%s', testortestgen, datestring);

% Load data.
load(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-supportFiles', filename), 'data_recog', 'data_recog_rt_mean', 'data_recog_acc_mean');

% Remove outliers.
keep = find(~ismember(data_recog_rt_mean.Var1, remove));
subjectlist = data_recog_rt_mean.Var1(keep);

rt = table2array(data_recog_rt_mean(keep, 2:end));
acc = table2array(data_recog_acc_mean(keep, 2:end));

alphastat = 0.66; % to return 1 SD, for 95% CI use .05

color_DI = [0.8500 0.3250 0.0980]; %orange 
coloralpha = .05;

% Set up plot and measure-specific details.
capsize = 0;
marker = 'o';
linewidth = 1.5;
linestyle = 'none';
markersize = 100;
xtickvalues = [1 2 3 4];
xlim_lo = min(xtickvalues)-0.5; xlim_hi = max(xtickvalues)+0.5;
fontname = 'Arial';
fontsize = 20;
fontangle = 'italic';
yticklength = 0;
xticklength = 0.05;

ylimlo = 0.5; ylimhi = 0.75;

% % Get individual subject means for each day.
% subjectlist = unique(data_recog.subID(find(~ismember(data_recog.subID, remove))));
% for sub = 1:length(subjectlist)
%     
%     for day = 1:length(unique(data_recog.day))
%     
%         clear idx;
%         idx = find(data_recog.subID == subjectlist(sub) & data_recog.day == day);
%         
%         if isempty(idx)
%             
%             rt(sub, day) = NaN;
%             acc(sub, day) = NaN;
% 
%         else
%             
%             rt(sub, day) = nanmean(data_recog.RT(idx));
%             acc(sub, day) = sum(data_recog.acc(idx))/length(data_recog.acc(idx));
% 
%         end
%     
%     end
%     
% end

%% Plot means and standard deviations, but including individual data.
clr = [randi(255, [10, 1]) randi(255, [10, 1]) randi(255, [10, 1])]./255;
figure(1)
hold on;

% Calculate mean and standard deviation across-subject, within-day.
DI_mean_rt = [nanmean(rt(:, 1)) nanmean(rt(:, 2)) nanmean(rt(:, 3)) nanmean(rt(:, 4))] ;
% DI_std_rt = [nanstd(rt(:, 1)) nanstd(rt(:, 2)) nanstd(rt(:, 3)) nanstd(rt(:, 4))] ;

% Plot means (do this first for legend and then second to keep it on top layer).
xval = linspace(1, length(DI_mean_rt), length(DI_mean_rt));
% scatter(xval, DI_mean, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_DI, 'MarkerEdgeColor', color_DI)

% Individual data points for DI.
gscatter(repmat(1, [size(rt, 1) 1]), rt(:, 1), subjectlist, clr)
gscatter(repmat(2, [size(rt, 1) 1]), rt(:, 2), subjectlist, clr)
gscatter(repmat(3, [size(rt, 1) 1]), rt(:, 3), subjectlist, clr)
gscatter(repmat(4, [size(rt, 1) 1]), rt(:, 4), subjectlist, clr)

% Means (second time to put it on top layer).
scatter(xval, DI_mean_rt, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_DI, 'MarkerEdgeColor', color_DI)
errorbar(xval, DI_mean_rt, DI_std_rt, 'Color', color_DI, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = xtickvalues;
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xlabels = {'Day 1', 'Day 2', 'Day 3', 'Day 4'};
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylimlo ylimhi];
yax.TickValues = [ylimlo (ylimlo+ylimhi)/2 ylimhi];
yax.TickDirection = 'out';
yax.TickLength = [yticklength yticklength];
yax.TickLabels = {num2str(ylimlo, '%2.2f'), num2str((ylimlo+ylimhi)/2, '%2.2f'), num2str(ylimhi, '%2.2f')};
yax.FontName = fontname;
yax.FontSize = fontsize;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

if strcmp(testortestgen, 'test')
    title({'Test'; '(typed symbols)'})
else
    title({'Generalization'; '(variable symbols)'})
end

legend('Location', 'eastoutside')
legend('boxoff');

a.YLabel.String = 'Reaction Time (seconds)';
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
pbaspect([1 1 1])

%     pos=get(gca,'Position');
%     pos1=pos-[0 .02 0 0];
%     set(gca,'Position', pos1);

% Write.
% if strcmp(save_figures, 'yes')
    
    print(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-plots', ['plot_recog_' testortestgen '_rt']), '-dpng')
    print(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-plots', 'eps', ['plot_recog_' testortestgen '_rt']), '-depsc')
    
% end

hold off;

%% Accuracy, letters: Plot means, but including individual data.
figure(2)
hold on;

ylimlo = 0.5; ylimhi = 1;

% Calculate mean and standard deviation across-subject, within-day.
DI_mean_acc = [nanmean(acc(:, 1)) nanmean(acc(:, 2)) nanmean(acc(:, 3)) nanmean(acc(:, 4))] ;
DI_std_acc = [nanstd(acc(:, 1)) nanstd(acc(:, 2)) nanstd(acc(:, 3)) nanstd(acc(:, 4))] ;

% Plot means (do this first for legend and then second to keep it on top layer).
xval = linspace(1, length(DI_mean_acc), length(DI_mean_acc));
% scatter(xval, DI_mean, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_DI, 'MarkerEdgeColor', color_DI)

% Individual data points for DI.
gscatter(repmat(1, [size(acc, 1) 1]), acc(:, 1), subjectlist, clr)
gscatter(repmat(2, [size(acc, 1) 1]), acc(:, 2), subjectlist, clr)
gscatter(repmat(3, [size(acc, 1) 1]), acc(:, 3), subjectlist, clr)
gscatter(repmat(4, [size(acc, 1) 1]), acc(:, 4), subjectlist, clr)

% Means (second time to put it on top layer).
scatter(xval, DI_mean_acc, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_DI, 'MarkerEdgeColor', color_DI)
errorbar(xval, DI_mean_acc, DI_std_acc, 'Color', color_DI, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = xtickvalues;
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xlabels = {'Day 1', 'Day 2', 'Day 3', 'Day 4'};
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylimlo ylimhi];
yax.TickValues = [ylimlo (ylimlo+ylimhi)/2 ylimhi];
yax.TickDirection = 'out';
yax.TickLength = [yticklength yticklength];
yax.TickLabels = {num2str(ylimlo, '%2.2f'), num2str((ylimlo+ylimhi)/2, '%2.2f'), num2str(ylimhi, '%2.2f')};
yax.FontName = fontname;
yax.FontSize = fontsize;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

if strcmp(testortestgen, 'test')
    title({'Test'; '(typed symbols)'})
else
    title({'Generalization'; '(variable symbols)'})
end

% legend({'Draw Ink', 'Draw No Ink', 'Watch Dynamic'}, 'Location', 'southeast', 'Orientation', 'vertical', 'FontSize', fontsize);
% legend('boxoff');

a.YLabel.String = 'Accuracy (percentage)';
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
pbaspect([1 1 1])

%     pos=get(gca,'Position');
%     pos1=pos-[0 .02 0 0];
%     set(gca,'Position', pos1);

% Write.
% if strcmp(save_figures, 'yes')
    
    print(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-plots', ['plot_recog_' testortestgen '_acc']), '-dpng')
    print(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-plots', 'eps', ['plot_recog_' testortestgen '_acc']), '-depsc')
    
% end

hold off;

% clear data_recog

