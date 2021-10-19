clear all; close all; clc

% Set working directories.
rootDir = '/Volumes/Seagate/wml/';

% Create date-specific file name that indicates how many subjects.
datestring = '20211018';
remove = [];%[22 27 35 41 60 62];
filename = sprintf('wml_beh_data_write_%s', datestring);

% Load data.
load(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-supportFiles', filename), 'data_write', 'data_write_mean');

% Remove outliers.
keep = find(~ismember(data_write_mean.Var1, remove));
subjectlist = data_write_mean.Var1(keep);

drawduration = table2array(data_write_mean(keep, 2:end));

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

ylimlo = 2.25; ylimhi = 3.5;

%% Plot means and standard deviations, but including individual data.
clr = [randi(255, [10, 1]) randi(255, [10, 1]) randi(255, [10, 1])]./255;
figure(1)
hold on;

% Calculate mean and standard deviation across-subject, within-day.
DI_mean = [nanmean(drawduration(:, 1)) nanmean(drawduration(:, 2)) nanmean(drawduration(:, 3)) nanmean(drawduration(:, 4))];
DI_median = [nanmedian(drawduration(:, 1)) nanmedian(drawduration(:, 2)) nanmedian(drawduration(:, 3)) nanmedian(drawduration(:, 4))];
DI_ci = 1.96.*[nanstd(drawduration(:, 1)) nanstd(drawduration(:, 2)) nanstd(drawduration(:, 3)) nanstd(drawduration(:, 4))]./sqrt(size(drawduration, 1));

% Plot means (do this first for legend and then second to keep it on top layer).
xval = linspace(1, length(DI_mean), length(DI_mean));
% scatter(xval, DI_mean, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_DI, 'MarkerEdgeColor', color_DI)

% Individual data points for DI.
gscatter(repmat(1, [size(drawduration, 1) 1]), drawduration(:, 1), subjectlist, clr)
gscatter(repmat(2, [size(drawduration, 1) 1]), drawduration(:, 2), subjectlist, clr)
gscatter(repmat(3, [size(drawduration, 1) 1]), drawduration(:, 3), subjectlist, clr)
gscatter(repmat(4, [size(drawduration, 1) 1]), drawduration(:, 4), subjectlist, clr)

% Means (second time to put it on top layer).
s(1)=scatter(xval, DI_mean, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_DI, 'MarkerEdgeColor', color_DI);
s(2)=scatter(xval, DI_median, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_DI, 'MarkerEdgeColor', color_DI, ...
    'MarkerFaceAlpha', 0.2);
errorbar(xval, DI_mean, DI_ci, 'Color', color_DI, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);
% errorbar(xval, DnI_mean, DnI_std, 'Color', color_DnI, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);

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
%

legend(s([1, 2]), {'Mean +/- 95% CI', 'Median'}, 'Location', 'southeast', 'Orientation', 'vertical', 'FontSize', fontsize);
legend('boxoff');
% legend('off')

n = size(drawduration, 1);
title(['Draw Duration, n = ' num2str(n)])

a.YLabel.String = 'Draw Duration (seconds)';
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
pbaspect([1 1 1])

%     pos=get(gca,'Position');
%     pos1=pos-[0 .02 0 0];
%     set(gca,'Position', pos1);

% Write.
% if strcmp(save_figures, 'yes')

print(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-plots', ['plot_write_n=' num2str(n)]), '-dpng')
print(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-plots', 'eps', ['plot_write_n=' num2str(n)]), '-depsc')

% end

hold off;

%% Plot rate of change, but including individual data.
