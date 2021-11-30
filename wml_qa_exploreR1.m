clear all; close all; clc;

% Set working directories.
rootdir = '/Volumes/Seagate/wml/wml-wmpredictslearning';
mp2rage = 'yes';
niter = 10000;
alphastat = 0.01;
% wmmeasure = 'r1';
% hemisphere = 'both'; %left, right, both

% Identify outliers for removal.
if strcmp(mp2rage, 'yes')
    remove = [27 32 22 24 25 26 28 30 31 33 35 42];
elseif strcmp(mp2rage, 'no')
    remove = [27 32];
end

% 51 is still processing on brainlife (with mp2rage data included)

% MRI
% 27 has severe motion
% 24, 31 withdrew mid-training, so multi-day learningrate is not accurate
% (but day 1 learning rate is accurate)
% 25, missing rightmdlfspl, rightmdlfang, rightvof
% 31, missing rightmdlfspl 
% 42, missing rightifof, leftuncinate
% 71, missing leftvof
% nomp2rage: 22 24 25 26 28 30 31 33 35 42
% non-perfect data quality: [27 28 30 31 32 33 35 41 42 45 46 47 50 51 54 56 57 59 60 64]; 

% Training
% 24 have only day 1 data
% 31 have only day 1 and day 2 data
% 32 had 0% accuracy by day 4
% 50 missing day 3, experimenter error

% Load recog data and remove outliers and sort by subID.
datestring = '20211119';
n = 36;
fittype = 'poly1';

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
xticklength = 0.02;

% Load recog acc data and remove outliers and sort by subID.
filename = sprintf('wml_beh_data_recog_test_%s', datestring);
load(fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', filename), 'data_recog', 'data_recog_acc_mean');

keep = find(~ismember(data_recog_acc_mean.Var1, remove));

temp = table2array(data_recog_acc_mean(keep, :));
acc = array2table(temp, 'VariableNames', {'subID', 'acc_day1', 'acc_day2', 'acc_day3', 'acc_day4', 'acc_day5'}); clear temp;
acc = sortrows(acc);

% Load writing data and remove outliers and sort by subID.
filename = sprintf('wml_beh_data_write_wlearningrate_%s_n=%d_%s', fittype, n, datestring);
load(fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', filename), 'data_write', 'data_write_mean');

keep = find(~ismember(data_write_mean.Var1, remove));

temp = table2array(data_write_mean(keep, :));
mt = array2table(temp, 'VariableNames', {'subID', 'mt_day1', 'mt_day2', 'mt_day3', 'mt_day4', 'mt_day5', 'mt_learningrate'}); clear temp;
mt = sortrows(mt, 1);
mt.mt_learningrate = -mt.mt_learningrate; % make faster learners positive (easier on my brain)

% Load recog rt data with learning rate and remove outliers and sort by subID.
filename = sprintf('wml_beh_data_recog_rt_wlearningrate_%s_n=%d_%s', fittype, n, datestring);
load(fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', filename), 'data_recog', 'data_recog_rt_mean');

keep = find(~ismember(data_recog_rt_mean.Var1, remove));

temp = table2array(data_recog_rt_mean(keep, :));
rt = array2table(temp, 'VariableNames', {'subID', 'rt_day1', 'rt_day2', 'rt_day3', 'rt_day4', 'rt_day5', 'rt_learningrate'}); clear temp;
rt = sortrows(rt, 1);
rt.rt_learningrate = -rt.rt_learningrate; % make faster learners positive (easier on my brain)

% Load tractprofiles data and remove outliers and sort by subID.
if strcmp(mp2rage, 'yes')
    filename = sprintf('wml_mri_data_tractprofiles_wmp2rage');%_%s', datestring);
elseif strcmp(mp2rage, 'no')
    filename = sprintf('wml_mri_data_tractprofiles_nomp2rage');%_%s', datestring);
end
load(fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', filename), 'data_tractprofiles_mean');

keep = find(~ismember(data_tractprofiles_mean.subID, remove));

temp = table2array(data_tractprofiles_mean(keep, :));
mri = array2table(temp, 'VariableNames', data_tractprofiles_mean.Properties.VariableNames); clear temp;
mri = sortrows(mri);

% Load streamline count data and remove outliers and sort by subID.
filename = sprintf('wml_wmlpredictslearning_data_streamlinecount');%_%s', datestring);
load(fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', filename), 'streamlinecounts');

keep = find(~ismember(streamlinecounts.subID, remove));
temp = table2array(streamlinecounts(keep, :));
ns = strcat('sc_', streamlinecounts.Properties.VariableNames);
tck = array2table(temp, 'VariableNames', ns); clear temp ns;
tck = sortrows(tck);

% Concatenate all variables into a table and zscore each column separately.
t = [mt rt(:, 2:end) acc(:, 2:end) mri(:, 2:end) tck(:, 2:end)];

%% Figures

figure(1);

tractname = 'leftparc';

% Fit a linear model: Does wm predict mt at day 1?
modelspec = 'leftparc_r1~sc_leftpArc';
mdl = fitlm(t, modelspec);
mdl.Coefficients
anova(mdl, 'summary')

% Another way, easier to plot with: Extract data for this fit.
y = t.leftparc_r1;
x = t.sc_leftpArc;

% Fit.
[f, f2] = fit(x, y, 'poly1'); %, 'Robust', 'Lar');

p=plot(f, x, y);
hold on;

% Customize scatter plot of data.
p(1).Color = [255 165 0]/255; %orange
p(1).MarkerSize = 40;
% Customize line plot of data.
p(2).Color = [255 165 0]/255; %orange
p(2).LineWidth = linewidth;
p(2).LineStyle = '--';

fontsize = 20;

xlim_lo = 0; xlim_hi = 7000;
ylim_lo = .5; ylim_hi = 1.5;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [0 1000 2000 3000 4000 5000 6000 7000];
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
% xax.FontAngle = fontangle;
% xax.TickLabelRotation = 90;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [.5 .75 1 1.25 1.5];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.TickLabels = ylabels;
yax.FontName = fontname;
yax.FontSize = fontsize;
% yax.TickLabelRotation = 0;

a = gca;
a.YLabel.String = {'R1'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Streamline Count'};
a.XLabel.FontSize = fontsize;
title(tractname)

box off;
legend('box', 'off');
legend('location', 'southeast');
pbaspect([1 1 1])

n = size(x, 1);
text(100, 1.4, ['adjr2 = ' num2str(f2.adjrsquare)])
text(100, 1.35, ['rmse = ' num2str(f2.rmse)])
text(100, 1.3, ['beta = ' num2str(f.p1)])
text(100, 1.25, ['p = ' num2str(mdl.Coefficients.pValue(2))])
text(100, 1.2, ['n = ' num2str(n)])

print(fullfile(rootdir, 'wml-wmpredictslearning-plots', ['plot_corr_r1_sc_' tractname '_n=' num2str(n)]), '-dpng')
print(fullfile(rootdir, 'wml-wmpredictslearning-plots', 'eps', ['plot_corr_r1_sc_' tractname '_n=' num2str(n)]), '-depsc')

figure(2);

tractname = 'rightparc';

% Fit a linear model: Does wm predict mt at day 1?
modelspec = 'rightparc_r1~sc_rightpArc';
mdl = fitlm(t, modelspec);
mdl.Coefficients
anova(mdl, 'summary')

% Another way, easier to plot with: Extract data for this fit.
y = t.rightparc_r1;
x = t.sc_rightpArc;

% Fit.
[f, f2] = fit(x, y, 'poly1'); %, 'Robust', 'Lar');

p=plot(f, x, y);
hold on;

% Customize scatter plot of data.
p(1).Color = [255 165 0]/255; %orange
p(1).MarkerSize = 40;
% Customize line plot of data.
p(2).Color = [255 165 0]/255; %orange
p(2).LineWidth = linewidth;
p(2).LineStyle = '--';

fontsize = 20;

xlim_lo = 0; xlim_hi = 7000;
ylim_lo = .5; ylim_hi = 1.5;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [0 1000 2000 3000 4000 5000 6000 7000];
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
% xax.FontAngle = fontangle;
% xax.TickLabelRotation = 90;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [.5 .75 1 1.25 1.5];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.TickLabels = ylabels;
yax.FontName = fontname;
yax.FontSize = fontsize;
% yax.TickLabelRotation = 0;

a = gca;
a.YLabel.String = {'R1'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Streamline Count'};
a.XLabel.FontSize = fontsize;
title(tractname)

box off;
legend('box', 'off');
legend('location', 'southeast');
pbaspect([1 1 1])

n = size(x, 1);
text(100, 1.4, ['adjr2 = ' num2str(f2.adjrsquare)])
text(100, 1.35, ['rmse = ' num2str(f2.rmse)])
text(100, 1.3, ['beta = ' num2str(f.p1)])
text(100, 1.25, ['p = ' num2str(mdl.Coefficients.pValue(2))])
text(100, 1.2, ['n = ' num2str(n)])

print(fullfile(rootdir, 'wml-wmpredictslearning-plots', ['plot_corr_r1_sc_' tractname '_n=' num2str(n)]), '-dpng')
print(fullfile(rootdir, 'wml-wmpredictslearning-plots', 'eps', ['plot_corr_r1_sc_' tractname '_n=' num2str(n)]), '-depsc')

figure(3);

tractname = 'lefttpc';

% Fit a linear model: Does wm predict mt at day 1?
modelspec = 'lefttpc_r1~sc_leftTPC';
mdl = fitlm(t, modelspec);
mdl.Coefficients
anova(mdl, 'summary')

% Another way, easier to plot with: Extract data for this fit.
y = t.lefttpc_r1;
x = t.sc_leftTPC;

% Fit.
[f, f2] = fit(x, y, 'poly1'); %, 'Robust', 'Lar');

p=plot(f, x, y);
hold on;

% Customize scatter plot of data.
p(1).Color = [255 165 0]/255; %orange
p(1).MarkerSize = 40;
% Customize line plot of data.
p(2).Color = [255 165 0]/255; %orange
p(2).LineWidth = linewidth;
p(2).LineStyle = '--';

fontsize = 20;

xlim_lo = 0; xlim_hi = 7000;
ylim_lo = .5; ylim_hi = 1.5;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [0 1000 2000 3000 4000 5000 6000 7000];
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
% xax.FontAngle = fontangle;
% xax.TickLabelRotation = 90;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [.5 .75 1 1.25 1.5];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.TickLabels = ylabels;
yax.FontName = fontname;
yax.FontSize = fontsize;
% yax.TickLabelRotation = 0;

a = gca;
a.YLabel.String = {'R1'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Streamline Count'};
a.XLabel.FontSize = fontsize;
title(tractname)

box off;
legend('box', 'off');
legend('location', 'southeast');
pbaspect([1 1 1])

n = size(x, 1);
text(100, 1.4, ['adjr2 = ' num2str(f2.adjrsquare)])
text(100, 1.35, ['rmse = ' num2str(f2.rmse)])
text(100, 1.3, ['beta = ' num2str(f.p1)])
text(100, 1.25, ['p = ' num2str(mdl.Coefficients.pValue(2))])
text(100, 1.2, ['n = ' num2str(n)])

print(fullfile(rootdir, 'wml-wmpredictslearning-plots', ['plot_corr_r1_sc_' tractname '_n=' num2str(n)]), '-dpng')
print(fullfile(rootdir, 'wml-wmpredictslearning-plots', 'eps', ['plot_corr_r1_sc_' tractname '_n=' num2str(n)]), '-depsc')

figure(4);

tractname = 'righttpc';

% Fit a linear model: Does wm predict mt at day 1?
modelspec = 'righttpc_r1~sc_rightTPC';
mdl = fitlm(t, modelspec);
mdl.Coefficients
anova(mdl, 'summary')

% Another way, easier to plot with: Extract data for this fit.
y = t.rightparc_r1;
x = t.sc_rightTPC;

% Fit.
[f, f2] = fit(x, y, 'poly1'); %, 'Robust', 'Lar');

p=plot(f, x, y);
hold on;

% Customize scatter plot of data.
p(1).Color = [255 165 0]/255; %orange
p(1).MarkerSize = 40;
% Customize line plot of data.
p(2).Color = [255 165 0]/255; %orange
p(2).LineWidth = linewidth;
p(2).LineStyle = '--';

fontsize = 20;

xlim_lo = 0; xlim_hi = 7000;
ylim_lo = .5; ylim_hi = 1.5;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [0 1000 2000 3000 4000 5000 6000 7000];
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
% xax.FontAngle = fontangle;
% xax.TickLabelRotation = 90;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [.5 .75 1 1.25 1.5];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.TickLabels = ylabels;
yax.FontName = fontname;
yax.FontSize = fontsize;
% yax.TickLabelRotation = 0;

a = gca;
a.YLabel.String = {'R1'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Streamline Count'};
a.XLabel.FontSize = fontsize;
title(tractname)

box off;
legend('box', 'off');
legend('location', 'southeast');
pbaspect([1 1 1])

n = size(x, 1);
text(100, 1.4, ['adjr2 = ' num2str(f2.adjrsquare)])
text(100, 1.35, ['rmse = ' num2str(f2.rmse)])
text(100, 1.3, ['beta = ' num2str(f.p1)])
text(100, 1.25, ['p = ' num2str(mdl.Coefficients.pValue(2))])
text(100, 1.2, ['n = ' num2str(n)])

print(fullfile(rootdir, 'wml-wmpredictslearning-plots', ['plot_corr_r1_sc_' tractname '_n=' num2str(n)]), '-dpng')
print(fullfile(rootdir, 'wml-wmpredictslearning-plots', 'eps', ['plot_corr_r1_sc_' tractname '_n=' num2str(n)]), '-depsc')

figure(5);

tractname = 'leftslf1and2';

% Fit a linear model: Does wm predict mt at day 1?
modelspec = 'leftslf1and2_r1~sc_leftSLF1And2';
mdl = fitlm(t, modelspec);
mdl.Coefficients
anova(mdl, 'summary')

% Another way, easier to plot with: Extract data for this fit.
y = t.leftslf1and2_r1;
x = t.sc_leftSLF1And2;

% Fit.
[f, f2] = fit(x, y, 'poly1'); %, 'Robust', 'Lar');

p=plot(f, x, y);
hold on;

% Customize scatter plot of data.
p(1).Color = [255 165 0]/255; %orange
p(1).MarkerSize = 40;
% Customize line plot of data.
p(2).Color = [255 165 0]/255; %orange
p(2).LineWidth = linewidth;
p(2).LineStyle = '--';

fontsize = 20;

xlim_lo = 8000; xlim_hi = 18000;
ylim_lo = .5; ylim_hi = 1.5;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [8000 10000 11000 12000 14000 16000 18000];
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
% xax.FontAngle = fontangle;
% xax.TickLabelRotation = 90;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [.5 .75 1 1.25 1.5];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.TickLabels = ylabels;
yax.FontName = fontname;
yax.FontSize = fontsize;
% yax.TickLabelRotation = 0;

a = gca;
a.YLabel.String = {'R1'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Streamline Count'};
a.XLabel.FontSize = fontsize;
title(tractname)

box off;
legend('box', 'off');
legend('location', 'southeast');
pbaspect([1 1 1])

n = size(x, 1);
text(100, 1.4, ['adjr2 = ' num2str(f2.adjrsquare)])
text(100, 1.35, ['rmse = ' num2str(f2.rmse)])
text(100, 1.3, ['beta = ' num2str(f.p1)])
text(100, 1.25, ['p = ' num2str(mdl.Coefficients.pValue(2))])
text(100, 1.2, ['n = ' num2str(n)])

print(fullfile(rootdir, 'wml-wmpredictslearning-plots', ['plot_corr_r1_sc_' tractname '_n=' num2str(n)]), '-dpng')
print(fullfile(rootdir, 'wml-wmpredictslearning-plots', 'eps', ['plot_corr_r1_sc_' tractname '_n=' num2str(n)]), '-depsc')

figure(6);

tractname = 'rightslf1and2';

% Fit a linear model: Does wm predict mt at day 1?
modelspec = 'rightslf1and2_r1~sc_rightSLF1And2';
mdl = fitlm(t, modelspec);
mdl.Coefficients
anova(mdl, 'summary')

% Another way, easier to plot with: Extract data for this fit.
y = t.rightslf1and2_r1;
x = t.sc_rightSLF1And2;

% Fit.
[f, f2] = fit(x, y, 'poly1'); %, 'Robust', 'Lar');

p=plot(f, x, y);
hold on;

% Customize scatter plot of data.
p(1).Color = [255 165 0]/255; %orange
p(1).MarkerSize = 40;
% Customize line plot of data.
p(2).Color = [255 165 0]/255; %orange
p(2).LineWidth = linewidth;
p(2).LineStyle = '--';

fontsize = 20;

xlim_lo = 8000; xlim_hi = 18000;
ylim_lo = .5; ylim_hi = 1.5;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [8000 10000 11000 12000 14000 16000 18000];
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
% xax.FontAngle = fontangle;
% xax.TickLabelRotation = 90;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [.5 .75 1 1.25 1.5];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.TickLabels = ylabels;
yax.FontName = fontname;
yax.FontSize = fontsize;
% yax.TickLabelRotation = 0;

a = gca;
a.YLabel.String = {'R1'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Streamline Count'};
a.XLabel.FontSize = fontsize;
title(tractname)

box off;
legend('box', 'off');
legend('location', 'southeast');
pbaspect([1 1 1])

n = size(x, 1);
text(100, 1.4, ['adjr2 = ' num2str(f2.adjrsquare)])
text(100, 1.35, ['rmse = ' num2str(f2.rmse)])
text(100, 1.3, ['beta = ' num2str(f.p1)])
text(100, 1.25, ['p = ' num2str(mdl.Coefficients.pValue(2))])
text(100, 1.2, ['n = ' num2str(n)])

print(fullfile(rootdir, 'wml-wmpredictslearning-plots', ['plot_corr_r1_sc_' tractname '_n=' num2str(n)]), '-dpng')
print(fullfile(rootdir, 'wml-wmpredictslearning-plots', 'eps', ['plot_corr_r1_sc_' tractname '_n=' num2str(n)]), '-depsc')

figure(7);

tractname = 'leftilf';

% Fit a linear model: Does wm predict mt at day 1?
modelspec = 'leftilf_r1~sc_leftILF';
mdl = fitlm(t, modelspec);
mdl.Coefficients
anova(mdl, 'summary')

% Another way, easier to plot with: Extract data for this fit.
y = t.leftilf_r1;
x = t.sc_leftILF;

% Fit.
[f, f2] = fit(x, y, 'poly1'); %, 'Robust', 'Lar');

p=plot(f, x, y);
hold on;

% Customize scatter plot of data.
p(1).Color = [255 165 0]/255; %orange
p(1).MarkerSize = 40;
% Customize line plot of data.
p(2).Color = [255 165 0]/255; %orange
p(2).LineWidth = linewidth;
p(2).LineStyle = '--';

fontsize = 20;

xlim_lo = 0; xlim_hi = 7000;
ylim_lo = .5; ylim_hi = 1.5;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [0 1000 2000 3000 4000 5000 6000 7000];
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
% xax.FontAngle = fontangle;
% xax.TickLabelRotation = 90;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [.5 .75 1 1.25 1.5];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.TickLabels = ylabels;
yax.FontName = fontname;
yax.FontSize = fontsize;
% yax.TickLabelRotation = 0;

a = gca;
a.YLabel.String = {'R1'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Streamline Count'};
a.XLabel.FontSize = fontsize;
title(tractname)

box off;
legend('box', 'off');
legend('location', 'southeast');
pbaspect([1 1 1])

n = size(x, 1);
text(100, 1.4, ['adjr2 = ' num2str(f2.adjrsquare)])
text(100, 1.35, ['rmse = ' num2str(f2.rmse)])
text(100, 1.3, ['beta = ' num2str(f.p1)])
text(100, 1.25, ['p = ' num2str(mdl.Coefficients.pValue(2))])
text(100, 1.2, ['n = ' num2str(n)])

print(fullfile(rootdir, 'wml-wmpredictslearning-plots', ['plot_corr_r1_sc_' tractname '_n=' num2str(n)]), '-dpng')
print(fullfile(rootdir, 'wml-wmpredictslearning-plots', 'eps', ['plot_corr_r1_sc_' tractname '_n=' num2str(n)]), '-depsc')

figure(8);

tractname = 'rightilf';

% Fit a linear model: Does wm predict mt at day 1?
modelspec = 'rightilf_r1~sc_rightILF';
mdl = fitlm(t, modelspec);
mdl.Coefficients
anova(mdl, 'summary')

% Another way, easier to plot with: Extract data for this fit.
y = t.rightilf_r1;
x = t.sc_rightILF;

% Fit.
[f, f2] = fit(x, y, 'poly1'); %, 'Robust', 'Lar');

p=plot(f, x, y);
hold on;

% Customize scatter plot of data.
p(1).Color = [255 165 0]/255; %orange
p(1).MarkerSize = 40;
% Customize line plot of data.
p(2).Color = [255 165 0]/255; %orange
p(2).LineWidth = linewidth;
p(2).LineStyle = '--';

fontsize = 20;

xlim_lo = 0; xlim_hi = 7000;
ylim_lo = .5; ylim_hi = 1.5;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [0 1000 2000 3000 4000 5000 6000 7000];
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
% xax.FontAngle = fontangle;
% xax.TickLabelRotation = 90;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [.5 .75 1 1.25 1.5];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.TickLabels = ylabels;
yax.FontName = fontname;
yax.FontSize = fontsize;
% yax.TickLabelRotation = 0;

a = gca;
a.YLabel.String = {'R1'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Streamline Count'};
a.XLabel.FontSize = fontsize;
title(tractname)

box off;
legend('box', 'off');
legend('location', 'southeast');
pbaspect([1 1 1])

n = size(x, 1);
text(100, 1.4, ['adjr2 = ' num2str(f2.adjrsquare)])
text(100, 1.35, ['rmse = ' num2str(f2.rmse)])
text(100, 1.3, ['beta = ' num2str(f.p1)])
text(100, 1.25, ['p = ' num2str(mdl.Coefficients.pValue(2))])
text(100, 1.2, ['n = ' num2str(n)])

print(fullfile(rootdir, 'wml-wmpredictslearning-plots', ['plot_corr_r1_sc_' tractname '_n=' num2str(n)]), '-dpng')
print(fullfile(rootdir, 'wml-wmpredictslearning-plots', 'eps', ['plot_corr_r1_sc_' tractname '_n=' num2str(n)]), '-depsc')

figure(9);

tractname = 'leftfat';

% Fit a linear model: Does wm predict mt at day 1?
modelspec = 'leftfat_r1~sc_leftAslant';
mdl = fitlm(t, modelspec);
mdl.Coefficients
anova(mdl, 'summary')

% Another way, easier to plot with: Extract data for this fit.
y = t.leftfat_r1;
x = t.sc_leftAslant;

% Fit.
[f, f2] = fit(x, y, 'poly1'); %, 'Robust', 'Lar');

p=plot(f, x, y);
hold on;

% Customize scatter plot of data.
p(1).Color = [255 165 0]/255; %orange
p(1).MarkerSize = 40;
% Customize line plot of data.
p(2).Color = [255 165 0]/255; %orange
p(2).LineWidth = linewidth;
p(2).LineStyle = '--';

fontsize = 20;

xlim_lo = 0; xlim_hi = 7000;
ylim_lo = .5; ylim_hi = 1.5;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [0 1000 2000 3000 4000 5000 6000 7000];
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
% xax.FontAngle = fontangle;
% xax.TickLabelRotation = 90;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [.5 .75 1 1.25 1.5];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.TickLabels = ylabels;
yax.FontName = fontname;
yax.FontSize = fontsize;
% yax.TickLabelRotation = 0;

a = gca;
a.YLabel.String = {'R1'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Streamline Count'};
a.XLabel.FontSize = fontsize;
title(tractname)

box off;
legend('box', 'off');
legend('location', 'southeast');
pbaspect([1 1 1])

n = size(x, 1);
text(100, 1.4, ['adjr2 = ' num2str(f2.adjrsquare)])
text(100, 1.35, ['rmse = ' num2str(f2.rmse)])
text(100, 1.3, ['beta = ' num2str(f.p1)])
text(100, 1.25, ['p = ' num2str(mdl.Coefficients.pValue(2))])
text(100, 1.2, ['n = ' num2str(n)])

print(fullfile(rootdir, 'wml-wmpredictslearning-plots', ['plot_corr_r1_sc_' tractname '_n=' num2str(n)]), '-dpng')
print(fullfile(rootdir, 'wml-wmpredictslearning-plots', 'eps', ['plot_corr_r1_sc_' tractname '_n=' num2str(n)]), '-depsc')

figure(10);

tractname = 'rightfat';

% Fit a linear model: Does wm predict mt at day 1?
modelspec = 'rightfat_r1~sc_rightAslant';
mdl = fitlm(t, modelspec);
mdl.Coefficients
anova(mdl, 'summary')

% Another way, easier to plot with: Extract data for this fit.
y = t.rightfat_r1;
x = t.sc_rightAslant;

% Fit.
[f, f2] = fit(x, y, 'poly1'); %, 'Robust', 'Lar');

p=plot(f, x, y);
hold on;

% Customize scatter plot of data.
p(1).Color = [255 165 0]/255; %orange
p(1).MarkerSize = 40;
% Customize line plot of data.
p(2).Color = [255 165 0]/255; %orange
p(2).LineWidth = linewidth;
p(2).LineStyle = '--';

fontsize = 20;

xlim_lo = 0; xlim_hi = 7000;
ylim_lo = .5; ylim_hi = 1.5;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [0 1000 2000 3000 4000 5000 6000 7000];
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
% xax.FontAngle = fontangle;
% xax.TickLabelRotation = 90;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [.5 .75 1 1.25 1.5];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.TickLabels = ylabels;
yax.FontName = fontname;
yax.FontSize = fontsize;
% yax.TickLabelRotation = 0;

a = gca;
a.YLabel.String = {'R1'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Streamline Count'};
a.XLabel.FontSize = fontsize;
title(tractname)

box off;
legend('box', 'off');
legend('location', 'southeast');
pbaspect([1 1 1])

n = size(x, 1);
text(100, 1.4, ['adjr2 = ' num2str(f2.adjrsquare)])
text(100, 1.35, ['rmse = ' num2str(f2.rmse)])
text(100, 1.3, ['beta = ' num2str(f.p1)])
text(100, 1.25, ['p = ' num2str(mdl.Coefficients.pValue(2))])
text(100, 1.2, ['n = ' num2str(n)])

print(fullfile(rootdir, 'wml-wmpredictslearning-plots', ['plot_corr_r1_sc_' tractname '_n=' num2str(n)]), '-dpng')
print(fullfile(rootdir, 'wml-wmpredictslearning-plots', 'eps', ['plot_corr_r1_sc_' tractname '_n=' num2str(n)]), '-depsc')

figure(11);

tractname = 'leftmdlfang';

% Fit a linear model: Does wm predict mt at day 1?
modelspec = 'leftmdlfang_r1~sc_leftMDLFang';
mdl = fitlm(t, modelspec);
mdl.Coefficients
anova(mdl, 'summary')

% Another way, easier to plot with: Extract data for this fit.
y = t.leftmdlfang_r1;
x = t.sc_leftMDLFang;

% Fit.
[f, f2] = fit(x, y, 'poly1'); %, 'Robust', 'Lar');

p=plot(f, x, y);
hold on;

% Customize scatter plot of data.
p(1).Color = [255 165 0]/255; %orange
p(1).MarkerSize = 40;
% Customize line plot of data.
p(2).Color = [255 165 0]/255; %orange
p(2).LineWidth = linewidth;
p(2).LineStyle = '--';

fontsize = 20;

xlim_lo = 0; xlim_hi = 7000;
ylim_lo = .5; ylim_hi = 1.5;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [0 1000 2000 3000 4000 5000 6000 7000];
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
% xax.FontAngle = fontangle;
% xax.TickLabelRotation = 90;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [.5 .75 1 1.25 1.5];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.TickLabels = ylabels;
yax.FontName = fontname;
yax.FontSize = fontsize;
% yax.TickLabelRotation = 0;

a = gca;
a.YLabel.String = {'R1'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Streamline Count'};
a.XLabel.FontSize = fontsize;
title(tractname)

box off;
legend('box', 'off');
legend('location', 'southeast');
pbaspect([1 1 1])

n = size(x, 1);
text(100, 1.4, ['adjr2 = ' num2str(f2.adjrsquare)])
text(100, 1.35, ['rmse = ' num2str(f2.rmse)])
text(100, 1.3, ['beta = ' num2str(f.p1)])
text(100, 1.25, ['p = ' num2str(mdl.Coefficients.pValue(2))])
text(100, 1.2, ['n = ' num2str(n)])

print(fullfile(rootdir, 'wml-wmpredictslearning-plots', ['plot_corr_r1_sc_' tractname '_n=' num2str(n)]), '-dpng')
print(fullfile(rootdir, 'wml-wmpredictslearning-plots', 'eps', ['plot_corr_r1_sc_' tractname '_n=' num2str(n)]), '-depsc')

figure(12);

tractname = 'rightmdlgang';

% Fit a linear model: Does wm predict mt at day 1?
modelspec = 'rightmdlfang_r1~sc_rightMDLFang';
mdl = fitlm(t, modelspec);
mdl.Coefficients
anova(mdl, 'summary')

% Another way, easier to plot with: Extract data for this fit.
y = t.rightmdlfang_r1;
x = t.sc_rightMDLFang;

% Fit.
[f, f2] = fit(x, y, 'poly1'); %, 'Robust', 'Lar');

p=plot(f, x, y);
hold on;

% Customize scatter plot of data.
p(1).Color = [255 165 0]/255; %orange
p(1).MarkerSize = 40;
% Customize line plot of data.
p(2).Color = [255 165 0]/255; %orange
p(2).LineWidth = linewidth;
p(2).LineStyle = '--';

fontsize = 20;

xlim_lo = 0; xlim_hi = 7000;
ylim_lo = .5; ylim_hi = 1.5;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [0 1000 2000 3000 4000 5000 6000 7000];
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
% xax.FontAngle = fontangle;
% xax.TickLabelRotation = 90;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [.5 .75 1 1.25 1.5];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.TickLabels = ylabels;
yax.FontName = fontname;
yax.FontSize = fontsize;
% yax.TickLabelRotation = 0;

a = gca;
a.YLabel.String = {'R1'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Streamline Count'};
a.XLabel.FontSize = fontsize;
title(tractname)

box off;
legend('box', 'off');
legend('location', 'southeast');
pbaspect([1 1 1])

n = size(x, 1);
text(100, 1.4, ['adjr2 = ' num2str(f2.adjrsquare)])
text(100, 1.35, ['rmse = ' num2str(f2.rmse)])
text(100, 1.3, ['beta = ' num2str(f.p1)])
text(100, 1.25, ['p = ' num2str(mdl.Coefficients.pValue(2))])
text(100, 1.2, ['n = ' num2str(n)])

print(fullfile(rootdir, 'wml-wmpredictslearning-plots', ['plot_corr_r1_sc_' tractname '_n=' num2str(n)]), '-dpng')
print(fullfile(rootdir, 'wml-wmpredictslearning-plots', 'eps', ['plot_corr_r1_sc_' tractname '_n=' num2str(n)]), '-depsc')