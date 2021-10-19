clear all; clc;

% Set working directories.
rootDir = '/Volumes/Seagate/wml/wml-wmpredictslearning';

% Identify outliers for removal.
remove = [49 58]; % 49, 58 are processing now on brainlife

% Load recog data and remove outliers and sort by subID.
datestring = '20211018';
n = 25;
fittype = 'poly1';

filename = sprintf('WML_beh_data_recog_test_%s', datestring);
load(fullfile(rootDir, 'wml-wmpredictslearning-supportFiles', filename), 'data_recog', 'data_recog_rt_mean', 'data_recog_acc_mean');

keep = find(~ismember(data_recog_rt_mean.Var1, remove));

temp = table2array(data_recog_rt_mean(keep, :));
rt = array2table(temp, 'VariableNames', {'subID', 'rt_day1', 'rt_day2', 'rt_day3', 'rt_day4', 'rt_day5'}); clear temp;
rt = sortrows(rt);

temp = table2array(data_recog_acc_mean(keep, :));
acc = array2table(temp, 'VariableNames', {'subID', 'acc_day1', 'acc_day2', 'acc_day3', 'acc_day4', 'acc_day5'}); clear temp;
acc = sortrows(acc);

% Load writing data and remove outliers and sort by subID.
datestring = '20211018';
filename = sprintf('WML_beh_data_write_wlearningrate_%s_n=%d_%s', fittype, n, datestring);
load(fullfile(rootDir, 'wml-wmpredictslearning-supportFiles', filename), 'data_write', 'data_write_mean');

keep = find(~ismember(data_write_mean.Var1, remove));

temp = table2array(data_write_mean(keep, :));
mt = array2table(temp, 'VariableNames', {'subID', 'mt_day1', 'mt_day2', 'mt_day3', 'mt_day4', 'mt_day5', 'learningrate'}); clear temp;
mt = sortrows(mt, 1);
mt.learningrate = -mt.learningrate; % make faster learners positive (easier on my brain)

% Load tractprofiles data and remove outliers and sort by subID.
datestring = '20211018';
filename = sprintf('WML_mri_data_tractprofiles');%_%s', datestring);
load(fullfile(rootDir, 'wml-wmpredictslearning-supportFiles', filename), 'data_tractprofiles_mean');

keep = find(~ismember(data_tractprofiles_mean.subID, remove));

temp = table2array(data_tractprofiles_mean(keep, :));
mri = array2table(temp, 'VariableNames', data_tractprofiles_mean.Properties.VariableNames); clear temp;
mri = sortrows(mri);

% QA: make histograms of each mri measure.
fa = table2array(mri(:, 2:6:end)); figure(1); histogram(fa(:), 100); xlim([0, 1]); title('fa'); hold off;
md = table2array(mri(:, 3:6:end)); figure(2); histogram(md(:), 100); xlim([min(md(:)), max(md(:))]); title('md'); hold off;
t1t2 = table2array(mri(:, 4:6:end)); figure(3); histogram(t1t2(:), 100); xlim([min(t1t2(:)), max(t1t2(:))]); title('t1t2'); hold off; 
ndi = table2array(mri(:, 5:6:end)); figure(4); histogram(ndi(:), 100); xlim([min(ndi(:)), max(ndi(:))]); title('ndi'); hold off;
odi = table2array(mri(:, 6:6:end)); figure(5); histogram(odi(:), 100); xlim([min(odi(:)), max(odi(:))]); title('odi'); hold off;
isovf = table2array(mri(:, 7:6:end)); figure(6); histogram(isovf(:), 100); xlim([min(isovf(:)), max(isovf(:))]); title('isovf'); hold off;

% Concatenate all variables into a table and zscore each column separately.
t = [mt rt(:, 2:end) acc(:, 2:end) mri(:, 2:end)];
t_temp = cat(2, table2array(t(:, 1)), (table2array(t(:, 2:end))-nanmean(table2array(t(:, 2:end)), 1))./nanstd(table2array(t(:, 2:end)), 1));
t = array2table(t_temp, 'VariableNames', t.Properties.VariableNames);

% Delete all rows that contain NaN, for now. Columns first.
idxc = [6 12 17]; %6, 12, and 17 correspond to day 5 test measurements
t(:, idxc) = [];
idxr = [1];% 2 3 8 9 13]; 
% 1 is subID
% 2 and 8 are subjs that did not complete all days of training
% 3 and 8 have NaN for rightmdlfspl
% 3, NaN for right mdlfang
% 9 leftifof
% 13 right ifof
t(idxr, :) = [];

% Correlations.
% figure(7);
% mat = corr(table2array(t(:, 2:end)));% % m = cat(2, table2array(mt(:, 2:end)), table2array(mri(:, 2:end)));
% toplot = mat(1:12, :);
% imagesc(toplot); colorbar; caxis([0.4 .9]);
% xlabels = t.Properties.VariableNames(2:end);
% ylabels = t.Properties.VariableNames(2:13);

% FA
figure(7);
mat = corr(table2array(t(:, [6, 15:6:end])));% % m = cat(2, table2array(mt(:, 2:end)), table2array(mri(:, 2:end)));
toplot = mat;
imagesc(toplot); colorbar; caxis([-1 1]);
xlabels = t.Properties.VariableNames([6, 15:6:end]);
ylabels = t.Properties.VariableNames([6, 15:6:end]);

% Set up plot and measure-specific details.
capsize = 0;
marker = 'o';
linewidth = 1.5;
linestyle = 'none';
markersize = 100;
xtickvalues = [1 2 3 4];
xlim_lo = min(xtickvalues)-0.5; xlim_hi = max(xtickvalues)+0.5;
fontname = 'Arial';
fontsize = 14;
fontangle = 'italic';
yticklength = 0;
xticklength = 0.05;

% xaxis
xax = get(gca, 'xaxis');
% xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = 1:size(toplot, 2);
xax.TickDirection = 'out';
% xax.TickLength = [xticklength xticklength];
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
% xax.FontAngle = fontangle;
xax.TickLabelRotation = 90;

% yaxis
yax = get(gca,'yaxis');
% yax.Limits = [ylimlo ylimhi];
yax.TickValues = 1:size(toplot, 1);
% yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
yax.TickLabels = ylabels;
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.TickLabelRotation = 0;

title('FA')


% mat = corr(m);
% imagesc(mat);
% colorbar;

% t1t2
figure(8);
mat = corr(table2array(t(:, [6, 17:6:end])));% % m = cat(2, table2array(mt(:, 2:end)), table2array(mri(:, 2:end)));
toplot = mat;
imagesc(toplot); colorbar; caxis([-1 1]);
xlabels = t.Properties.VariableNames([6, 17:6:end]);
ylabels = t.Properties.VariableNames([6, 17:6:end]);

% Set up plot and measure-specific details.
capsize = 0;
marker = 'o';
linewidth = 1.5;
linestyle = 'none';
markersize = 100;
xtickvalues = [1 2 3 4];
xlim_lo = min(xtickvalues)-0.5; xlim_hi = max(xtickvalues)+0.5;
fontname = 'Arial';
fontsize = 14;
fontangle = 'italic';
yticklength = 0;
xticklength = 0.05;

% xaxis
xax = get(gca, 'xaxis');
% xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = 1:size(toplot, 2);
xax.TickDirection = 'out';
% xax.TickLength = [xticklength xticklength];
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
% xax.FontAngle = fontangle;
xax.TickLabelRotation = 90;

% yaxis
yax = get(gca,'yaxis');
% yax.Limits = [ylimlo ylimhi];
yax.TickValues = 1:size(toplot, 1);
% yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
yax.TickLabels = ylabels;
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.TickLabelRotation = 0;

title('t1/t2 ratio')


% mat = corr(m);
% imagesc(mat);
% colorbar;

% ndi
figure(9);
mat = corr(table2array(t(:, [6, 18:6:end])));% % m = cat(2, table2array(mt(:, 2:end)), table2array(mri(:, 2:end)));
toplot = mat;
imagesc(toplot); colorbar; caxis([-1 1]);
xlabels = t.Properties.VariableNames([6, 18:6:end]);
ylabels = t.Properties.VariableNames([6, 18:6:end]);

% Set up plot and measure-specific details.
capsize = 0;
marker = 'o';
linewidth = 1.5;
linestyle = 'none';
markersize = 100;
xtickvalues = [1 2 3 4];
xlim_lo = min(xtickvalues)-0.5; xlim_hi = max(xtickvalues)+0.5;
fontname = 'Arial';
fontsize = 14;
fontangle = 'italic';
yticklength = 0;
xticklength = 0.05;

% xaxis
xax = get(gca, 'xaxis');
% xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = 1:size(toplot, 2);
xax.TickDirection = 'out';
% xax.TickLength = [xticklength xticklength];
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
% xax.FontAngle = fontangle;
xax.TickLabelRotation = 90;

% yaxis
yax = get(gca,'yaxis');
% yax.Limits = [ylimlo ylimhi];
yax.TickValues = 1:size(toplot, 1);
% yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
yax.TickLabels = ylabels;
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.TickLabelRotation = 0;

title('NDI')


% mat = corr(m);
% imagesc(mat);
% colorbar;


% Fit a linear model: Does wm predict mt at day 1?
modelspec = 'learningrate~rightparc_fa';
mdl = fitlm(t, modelspec);
mdl.Coefficients
anova(mdl, 'summary')




