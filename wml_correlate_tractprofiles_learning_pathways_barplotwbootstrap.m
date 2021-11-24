clear all; close all; clc;

% Set working directories.
rootDir = '/Volumes/Seagate/wml/wml-wmpredictslearning';
mp2rage = 'no';
niter = 10000;
alphastat = 0.01;

% Identify outliers for removal.
if strcmp(mp2rage, 'yes')
    remove = [27 32 51 22 24 25 26 28 30 31 33 35 42];
elseif strcmp(mp2rage, 'no')
    remove = [27 32 51];
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
remove = [51];
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
load(fullfile(rootDir, 'wml-wmpredictslearning-supportFiles', filename), 'data_recog', 'data_recog_acc_mean');

keep = find(~ismember(data_recog_acc_mean.Var1, remove));

temp = table2array(data_recog_acc_mean(keep, :));
acc = array2table(temp, 'VariableNames', {'subID', 'acc_day1', 'acc_day2', 'acc_day3', 'acc_day4', 'acc_day5'}); clear temp;
acc = sortrows(acc);

% Load writing data and remove outliers and sort by subID.
filename = sprintf('wml_beh_data_write_wlearningrate_%s_n=%d_%s', fittype, n, datestring);
load(fullfile(rootDir, 'wml-wmpredictslearning-supportFiles', filename), 'data_write', 'data_write_mean');

keep = find(~ismember(data_write_mean.Var1, remove));

temp = table2array(data_write_mean(keep, :));
mt = array2table(temp, 'VariableNames', {'subID', 'mt_day1', 'mt_day2', 'mt_day3', 'mt_day4', 'mt_day5', 'mt_learningrate'}); clear temp;
mt = sortrows(mt, 1);
mt.mt_learningrate = -mt.mt_learningrate; % make faster learners positive (easier on my brain)

% Load recog rt data with learning rate and remove outliers and sort by subID.
filename = sprintf('wml_beh_data_recog_rt_wlearningrate_%s_n=%d_%s', fittype, n, datestring);
load(fullfile(rootDir, 'wml-wmpredictslearning-supportFiles', filename), 'data_recog', 'data_recog_rt_mean');

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
load(fullfile(rootDir, 'wml-wmpredictslearning-supportFiles', filename), 'data_tractprofiles_mean');

keep = find(~ismember(data_tractprofiles_mean.subID, remove));

temp = table2array(data_tractprofiles_mean(keep, :));
mri = array2table(temp, 'VariableNames', data_tractprofiles_mean.Properties.VariableNames); clear temp;
mri = sortrows(mri);

% Concatenate all variables into a table and zscore each column separately.
t = [mt rt(:, 2:end) acc(:, 2:end) mri(:, 2:end)];
t_temp = cat(2, table2array(t(:, 1)), (table2array(t(:, 2:end))-nanmean(table2array(t(:, 2:end)), 1))./nanstd(table2array(t(:, 2:end)), 1));
t = array2table(t_temp, 'VariableNames', t.Properties.VariableNames);

% Delete all columns that contain NaN, for now... essentially removing day 5 behavioral measures.
idxc = [1 6 12 18]; % 1 is subID, 6, 12, and 18 correspond to day 5 test measurements
t(:, idxc) = [];
% idxr = [3 11 34]; %[3 16]; %[2 9 15];
% t(idxr, :) = [];

%% Posterior Vertical Pathway -- Sensorimotor
lr = t.mt_learningrate;
wm = mean([t.leftparc_fa t.leftmdlfang_fa t.leftmdlfspl_fa t.lefttpc_fa], 2);
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_pvp_ss = mean(beta_resampled);
ci_pvp_ss = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled

%% Posterior Vertical Pathway -- Recognition
lr = t.rt_learningrate;
wm = mean([t.leftparc_fa t.leftmdlfang_fa t.leftmdlfspl_fa t.lefttpc_fa], 2);
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_pvp_recog = mean(beta_resampled);
ci_pvp_recog = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled

%% Ventral Horizontal Pathway -- Sensorimotor
lr = t.mt_learningrate;
wm = mean([t.leftilf_fa t.leftifof_fa], 2);
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_vhp_ss = mean(beta_resampled);
ci_vhp_ss = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled

%% Ventral Horizontal Pathway -- Recognition
lr = t.rt_learningrate;
wm = mean([t.leftilf_fa t.leftifof_fa], 2);
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_vhp_recog = mean(beta_resampled);
ci_vhp_recog = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled

%% Dorsal Horizontal Pathway -- Sensorimotor
lr = t.mt_learningrate;
wm = mean([t.leftslf1and2_fa t.leftslf3_fa], 2);
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_dhp_ss = mean(beta_resampled);
ci_dhp_ss = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled

%% Ventral Horizontal Pathway -- Recognition
lr = t.rt_learningrate;
wm = mean([t.leftslf1and2_fa t.leftslf3_fa], 2);
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_dhp_recog = mean(beta_resampled);
ci_dhp_recog = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled

%% Plot
% Specify colors.
vhp = [0 127 255]/255; % dark blue
dhp = [204 190 0]/255; % burnt yellow
pvp = [2 129 129]/255; % dark turquoise

figure(1); hold on;
b = bar([1 5], [beta_dhp_ss beta_dhp_recog], 0.15);
color = dhp;
b.FaceColor = color; b.EdgeColor = color; 
p = plot([1 1], ci_dhp_ss); p(1).Color = color/2;
p = plot([5 5], ci_dhp_recog); p(1).Color = color/2;

b = bar([2 6], [beta_pvp_ss beta_pvp_recog], 0.15);
color = pvp;
b.FaceColor = color; b.EdgeColor = color; 
p = plot([2 2], ci_pvp_ss); p(1).Color = color/2;
p = plot([6 6], ci_pvp_recog); p(1).Color = color/2;

b = bar([3 7], [beta_vhp_ss beta_vhp_recog], 0.15);
color = vhp;
b.FaceColor = color; b.EdgeColor = color; 
p = plot([3 3], ci_vhp_ss); p(1).Color = color/2;
p = plot([7 7], ci_vhp_recog); p(1).Color = color/2;

xlim([0.5 7.5]);

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0.5 7.5];
xax.TickValues = [2 6];
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = {'Sensorimotor Learning', 'Recognition Learning'};
xax.FontName = fontname;
xax.FontSize = fontsize;
% xax.FontAngle = fontangle;
% xax.TickLabelRotation = 90;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [-.5 1];
yax.TickValues = [-1 -.5 0 .5 1];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.TickLabels = ylabels;
yax.FontName = fontname;
yax.FontSize = fontsize;
% yax.TickLabelRotation = 0;

% title('Left pArc')
a = gca;
a.YLabel.String = {'Predicton Strength (beta)'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {' '};
% a.XLabel.FontSize = fontsize;

box off;
legend('box', 'off');
legend('location', 'northeast');
legend({'DHP', '', '', 'PVP', '', '', 'VHP'})
pbaspect([1 1 1])

% n = size(x, 1);
% text(-1.8, 1.8, ['adjr2 = ' num2str(f2.adjrsquare)])
% text(-1.8, 1.6, ['rmse = ' num2str(f2.rmse)])
% text(-1.8, 1.4, ['beta = ' num2str(f.p1)])
% text(-1.8, 1.2, ['p = ' num2str(mdl.Coefficients.pValue(2))])
% text(-1.8, 1.0, ['n = ' num2str(n)])

print(fullfile(rootDir, 'wml-wmpredictslearning-plots', ['plot_pred_learning_DHPPVPVHP_n=' num2str(n)]), '-dpng')
print(fullfile(rootDir, 'wml-wmpredictslearning-plots', 'eps', ['plot_pred_learning_DHPPVPVHP_n=' num2str(n)]), '-depsc')



