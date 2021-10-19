clear all; close all; clc

% Set working directories.
rootDir = '/Volumes/Seagate/wml/';
remove = [];%[32];

% Create date-specific file name that indicates how many subjects.
datestring = '20211018';

% Load test data.
testortestgen = 'test'; %'test' or 'test_gen'
filename = sprintf('wml_beh_data_recog_%s_%s', testortestgen, datestring);
load(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-supportFiles', filename), 'data_recog', 'data_recog_rt_mean', 'data_recog_acc_mean');
test = data_recog; clear data_recog;

% Load test gen data.
testortestgen = 'testgen'; %'test' or 'testgen'
filename = sprintf('wml_beh_data_recog_%s_%s', testortestgen, datestring);
load(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-supportFiles', filename), 'data_recog', 'data_recog_rt_mean', 'data_recog_acc_mean');
testgen = data_recog; clear data_recog;
clear testortestgen;

alphastat = 0.66; % to return 1 SD, for 95% CI use .05

color_DI = [0.8500 0.3250 0.0980]; %orange 
coloralpha = .05;

% Set up plot and measure-specific details.
capsize = 0;
marker = 'o';
linewidth = 1.5;
linestyle = 'none';
markersize = 36;
xtickvalues = [1 2 3 4];
xlim_lo = min(xtickvalues)-0.5; xlim_hi = max(xtickvalues)+0.5;
fontname = 'Arial';
fontsize = 20;
fontangle = 'italic';
yticklength = 0;
xticklength = 0.05;

ylimlo = -0.1; ylimhi = 0.12;

% Get individual subject means for each day.
subjectlist = unique(test.subID(find(test.subID ~= remove)));
for sub = 1:length(subjectlist)
    
    for day = 1:length(unique(test.day))
    
        clear idx;
        idx_test = find(test.subID == subjectlist(sub) & test.day == day);
        idx_testgen = find(testgen.subID == subjectlist(sub) & testgen.day == day);

        if isempty(idx_test) || isempty(idx_testgen)
            
            rt(sub, day) = NaN;
            acc(sub, day) = NaN;

        else
            
            rt(sub, day) = (nanmean(test.RT(idx_test)) - nanmean(testgen.RT(idx_testgen)));
            acc(sub, day) = (sum(test.acc(idx_test))/length(test.acc(idx_test)) - sum(testgen.acc(idx_testgen))/length(testgen.acc(idx_testgen)));

        end
    
    end
    
end

%% Plot means and standard deviations, but including individual data.
clr = [randi(255, [10, 1]) randi(255, [10, 1]) randi(255, [10, 1])]./255;
figure(1)
hold on;

% Calculate mean and standard deviation across-subject, within-day.
DI_mean_rt = [nanmean(rt(:, 1)) nanmean(rt(:, 2)) nanmean(rt(:, 3)) nanmean(rt(:, 4))] ;
DI_std_rt = [nanstd(rt(:, 1)) nanstd(rt(:, 2)) nanstd(rt(:, 3)) nanstd(rt(:, 4))] ;

% Plot means (do this first for legend and then second to keep it on top layer).
xval = linspace(1, length(DI_mean_rt), length(DI_mean_rt));
% scatter(xval, DI_mean, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_DI, 'MarkerEdgeColor', color_DI)

% Individual data points for DI.
gscatter(repmat(1, [size(rt, 1) 1]), rt(:, 1), subjectlist, clr, '.', markersize)
gscatter(repmat(2, [size(rt, 1) 1]), rt(:, 2), subjectlist, clr, '.', markersize)
gscatter(repmat(3, [size(rt, 1) 1]), rt(:, 3), subjectlist, clr, '.', markersize)
gscatter(repmat(4, [size(rt, 1) 1]), rt(:, 4), subjectlist, clr, '.', markersize)

% % Means (second time to put it on top layer).
% scatter(xval, DI_mean_rt, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_DI, 'MarkerEdgeColor', color_DI)
% errorbar(xval, DI_mean_rt, DI_std_rt, 'Color', color_DI, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);

plot([0.5 4.5], [0 0], 'k')

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

title({'Generalization Ability'; '(typed - variable symbols)'})

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
    
    print(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-plots', 'plot_recog_genability_rt'), '-dpng')
    print(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-plots', 'eps', 'plot_recog_genability_rt'), '-depsc')
    
% end

hold off;

%% Accuracy, letters: Plot means, but including individual data.
figure(2)
hold on;

ylimlo = -0.1; ylimhi = 0.25;

% Calculate mean and standard deviation across-subject, within-day.
DI_mean_acc = [nanmean(acc(:, 1)) nanmean(acc(:, 2)) nanmean(acc(:, 3)) nanmean(acc(:, 4))] ;
DI_std_acc = [nanstd(acc(:, 1)) nanstd(acc(:, 2)) nanstd(acc(:, 3)) nanstd(acc(:, 4))] ;

% Plot means (do this first for legend and then second to keep it on top layer).
xval = linspace(1, length(DI_mean_acc), length(DI_mean_acc));
% scatter(xval, DI_mean, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_DI, 'MarkerEdgeColor', color_DI)

% Individual data points for DI.
gscatter(repmat(1, [size(acc, 1) 1]), acc(:, 1), subjectlist, clr, '.', markersize)
gscatter(repmat(2, [size(acc, 1) 1]), acc(:, 2), subjectlist, clr, '.', markersize)
gscatter(repmat(3, [size(acc, 1) 1]), acc(:, 3), subjectlist, clr, '.', markersize)
gscatter(repmat(4, [size(acc, 1) 1]), acc(:, 4), subjectlist, clr, '.', markersize)

% Means (second time to put it on top layer).
% scatter(xval, DI_mean_acc, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_DI, 'MarkerEdgeColor', color_DI)
% errorbar(xval, DI_mean_acc, DI_std_acc, 'Color', color_DI, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);

plot([0.5 4.5], [0 0], 'k')

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

title({'Generalization Ability'; '(typed - variable symbols)'})

legend('Location', 'eastoutside')
legend('boxoff');

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
    
    print(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-plots', 'plot_recog_genability_acc'), '-dpng')
    print(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-plots', 'eps', 'plot_recog_genability_acc'), '-depsc')
    
% end

hold off;

% clear data_recog

