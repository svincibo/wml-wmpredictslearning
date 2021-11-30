clear all; close all; clc;

% Set working directories.
rootDir = '/Volumes/Seagate/wml/wml-wmpredictslearning';
mp2rage = 'yes';
niter = 10000;
alphastat = 0.01;
wmmeasure = 'r1';
hemisphere = 'both'; %left, right, both

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

%% Frontal motor -- vof -- Sensorimotor
lr = t.mt_learningrate;
if strcmp(hemisphere, 'left')
    wm = t.leftvof_r1;
elseif strcmp(hemisphere, 'right')
        wm = t.rightvof_r1;
elseif strcmp(hemisphere, 'both')
        wm = mean([t.leftvof_r1 t.rightvof_r1], 2);
end
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_vof_ss = mean(beta_resampled);
ci_vof_ss = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled

%% Frontal Motor -- vof -- Recognition
lr = t.rt_learningrate;
if strcmp(hemisphere, 'left')
    wm = t.leftvof_r1;
elseif strcmp(hemisphere, 'right')
        wm = t.rightvof_r1;
elseif strcmp(hemisphere, 'both')
            wm = mean([t.leftvof_r1 t.rightvof_r1], 2);
end
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_vof_recog = mean(beta_resampled);
ci_vof_recog = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled

%% Posterior Vertical Pathway -- parc -- Sensorimotor
lr = t.mt_learningrate;
if strcmp(hemisphere, 'left')
    wm = t.leftparc_r1;
elseif strcmp(hemisphere, 'right')
        wm = t.rightparc_r1;
elseif strcmp(hemisphere, 'both')
            wm = mean([t.leftparc_r1 t.rightparc_r1], 2);
end
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_parc_ss = mean(beta_resampled);
ci_parc_ss = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled

%% Posterior Vertical Pathway -- parc -- Recognition
lr = t.rt_learningrate;
if strcmp(hemisphere, 'left')
wm = t.leftparc_r1;
elseif strcmp(hemisphere, 'right')
    wm = t.rightparc_r1;
elseif strcmp(hemisphere, 'both')
            wm = mean([t.leftparc_r1 t.rightparc_r1], 2);
end
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_parc_recog = mean(beta_resampled);
ci_parc_recog = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled

%% Posterior Vertical Pathway -- tpc -- Sensorimotor
lr = t.mt_learningrate;
if strcmp(hemisphere, 'left')
wm = t.lefttpc_r1;
elseif strcmp(hemisphere, 'right')
    wm = t.righttpc_r1;
elseif strcmp(hemisphere, 'both')
            wm = mean([t.lefttpc_r1 t.righttpc_r1], 2);
end
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_tpc_ss = mean(beta_resampled);
ci_tpc_ss = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled

%% Posterior Vertical Pathway -- tpc -- Recognition
lr = t.rt_learningrate;
if strcmp(hemisphere, 'left')
wm = t.lefttpc_r1;
elseif strcmp(hemisphere, 'right')
    wm = t.righttpc_r1;
elseif strcmp(hemisphere, 'both')
            wm = mean([t.lefttpc_r1 t.righttpc_r1], 2);
end
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_tpc_recog = mean(beta_resampled);
ci_tpc_recog = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled

%% Posterior Vertical Pathway -- mdlfang -- Sensorimotor
lr = t.mt_learningrate;
if strcmp(hemisphere, 'left')
wm = t.leftmdlfang_r1;
elseif strcmp(hemisphere, 'right')
    wm = t.rightmdlfang_r1;
elseif strcmp(hemisphere, 'both')
            wm = mean([t.leftmdlfang_r1 t.rightmdlfang_r1], 2);
end
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_mdlfang_ss = mean(beta_resampled);
ci_mdlfang_ss = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled

%% Posterior Vertical Pathway -- mdlfang -- Recognition
lr = t.rt_learningrate;
if strcmp(hemisphere, 'left')
wm = t.leftmdlfang_r1;
elseif strcmp(hemisphere, 'right')
    wm = t.rightmdlfang_r1;
elseif strcmp(hemisphere, 'both')
            wm = mean([t.leftmdlfang_r1 t.rightmdlfang_r1], 2);
end
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_mdlfang_recog = mean(beta_resampled);
ci_mdlfang_recog = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled

%% Posterior Vertical Pathway -- mdlfspl -- Sensorimotor
lr = t.mt_learningrate;
if strcmp(hemisphere, 'left')
wm = t.leftmdlfspl_r1;
elseif strcmp(hemisphere, 'right')
    wm = t.rightmdlfspl_r1;
elseif strcmp(hemisphere, 'both')
            wm = mean([t.leftmdlfspl_r1 t.rightmdlfspl_r1], 2);
end
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_mdlfspl_ss = mean(beta_resampled);
ci_mdlfspl_ss = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled

%% Posterior Vertical Pathway -- mdlfspl -- Recognition
lr = t.rt_learningrate;
if strcmp(hemisphere, 'left')
wm = t.leftmdlfspl_r1;
elseif strcmp(hemisphere, 'right')
    wm = t.rightmdlfspl_r1;
elseif strcmp(hemisphere, 'both')
            wm = mean([t.leftmdlfspl_r1 t.rightmdlfspl_r1], 2);
end
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_mdlfspl_recog = mean(beta_resampled);
ci_mdlfspl_recog = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled

%% Ventral Horizontal Pathway -- ilf -- Sensorimotor
lr = t.mt_learningrate;
if strcmp(hemisphere, 'left')
wm = t.leftilf_r1;
elseif strcmp(hemisphere, 'right')
    wm = t.rightilf_r1;
elseif strcmp(hemisphere, 'both')
            wm = mean([t.leftilf_r1 t.rightilf_r1], 2);
end
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_ilf_ss = mean(beta_resampled);
ci_ilf_ss = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled

%% Ventral Horizontal Pathway -- ilf -- Recognition
lr = t.rt_learningrate;
if strcmp(hemisphere, 'left')
wm = t.leftilf_r1;
elseif strcmp(hemisphere, 'right')
    wm = t.rightilf_r1;
elseif strcmp(hemisphere, 'both')
            wm = mean([t.leftilf_r1 t.rightilf_r1], 2);
end
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_ilf_recog = mean(beta_resampled);
ci_ilf_recog = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled

%% Ventral Horizontal Pathway -- ifof -- Sensorimotor
lr = t.mt_learningrate;
if strcmp(hemisphere, 'left')
wm = t.leftifof_r1;
elseif strcmp(hemisphere, 'right')
    wm = t.rightifof_r1;
elseif strcmp(hemisphere, 'both')
            wm = mean([t.leftifof_r1 t.rightifof_r1], 2);
end
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_ifof_ss = mean(beta_resampled);
ci_ifof_ss = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled

%% Ventral Horizontal Pathway -- ifof -- Recognition
lr = t.rt_learningrate;
if strcmp(hemisphere, 'left')
wm = t.leftifof_r1;
elseif strcmp(hemisphere, 'right')
    wm = t.rightifof_r1;
elseif strcmp(hemisphere, 'both')
            wm = mean([t.leftifof_r1 t.rightifof_r1], 2);
end
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_ifof_recog = mean(beta_resampled);
ci_ifof_recog = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled


%% Dorsal Horizontal Pathway -- slf1and2 -- Sensorimotor
lr = t.mt_learningrate;
if strcmp(hemisphere, 'left')
wm = t.leftslf1and2_r1;
elseif strcmp(hemisphere, 'right')
    wm = t.rightslf1and2_r1;
elseif strcmp(hemisphere, 'both')
            wm = mean([t.leftslf1and2_r1 t.rightslf1and2_r1], 2);
end
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_slf1and2_ss = mean(beta_resampled);
ci_slf1and2_ss = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled

%% Dorsal Horizontal Pathway -- slf1and2 -- Recognition
lr = t.rt_learningrate;
if strcmp(hemisphere, 'left')
wm = t.leftslf1and2_r1;
elseif strcmp(hemisphere, 'right')
    wm = t.rightslf1and2_r1;
elseif strcmp(hemisphere, 'both')
            wm = mean([t.leftslf1and2_r1 t.rightslf1and2_r1], 2);
end
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_slf1and2_recog = mean(beta_resampled);
ci_slf1and2_recog = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled

%% Dorsal Horizontal Pathway -- slf3 -- Sensorimotor
lr = t.mt_learningrate;
if strcmp(hemisphere, 'left')
wm = t.leftslf3_r1;
elseif strcmp(hemisphere, 'right')
    wm = t.rightslf3_r1;
elseif strcmp(hemisphere, 'both')
            wm = mean([t.leftslf3_r1 t.rightslf3_r1], 2);
end
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_slf3_ss = mean(beta_resampled);
ci_slf3_ss = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled

%% Dorsal Horizontal Pathway -- slf3 -- Recognition
lr = t.rt_learningrate;
if strcmp(hemisphere, 'left')
wm = t.leftslf3_r1;
elseif strcmp(hemisphere, 'right')
    wm = t.rightslf3_r1;
elseif strcmp(hemisphere, 'both')
            wm = mean([t.leftslf3_r1 t.rightslf3_r1], 2);
end
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_slf3_recog = mean(beta_resampled);
ci_slf3_recog = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled

%% Frontal motor -- fat -- Sensorimotor
lr = t.mt_learningrate;
if strcmp(hemisphere, 'left')
wm = t.leftfat_r1;
elseif strcmp(hemisphere, 'right')
    wm = t.rightfat_r1;
elseif strcmp(hemisphere, 'both')
            wm = mean([t.leftfat_r1 t.rightfat_r1], 2);
end
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_fat_ss = mean(beta_resampled);
ci_fat_ss = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled

%% Frontal Motor -- fat -- Recognition
lr = t.rt_learningrate;
if strcmp(hemisphere, 'left')
    wm = t.leftfat_r1;
elseif strcmp(hemisphere, 'right')
    wm = t.rightfat_r1;
elseif strcmp(hemisphere, 'both')
            wm = mean([t.leftfat_r1 t.leftfat_r1], 2);
end
for r = 1:niter
    
    idx = randsample(1:size(t, 1), size(t, 1), true);
    this_lr = lr(idx);
    this_wm = wm(idx);
    
    mdl = fitlm(this_wm, this_lr, 'poly1');
    
    beta_resampled(r) = mdl.Coefficients.Estimate(end);
    
    clear idx this_lr this_wm mdl
      
end

beta_fat_recog = mean(beta_resampled);
ci_fat_recog = prctile(beta_resampled, [100*alphastat/2, 100*(1-alphastat/2)]);
clear beta_resampled

%% Plot
% Specify colors.
vhp = [0 127 255]/255; % dark blue
dhp = [204 190 0]/255; % burnt yellow
pvp = [2 129 129]/255; % dark turquoise

% Specify colors.
ifof = [142 198 255]/255; % light blue
ilf = [0 127 255]/255; % dark blue

slf1and2 = [237 177 32]/255; % burnt yellow
slf3 = [204 204 0]/255; % light burnt yellow

parc = [64 224 208]/255; % turquoise
tpc =  [27 102 87]/255; % dark turquoise
mdlfspl = [42, 102, 0]/255; % green
mdlfang =  [60 179 113]/255; % medium sea green

vof = [147 112 219]/255; % medium purple

fat = [240 128 128]/255; % light coral

figure(1); hold on;
b = bar([1 12], [beta_fat_ss beta_fat_recog], 0.075);
color = fat;
b.FaceColor = color; b.EdgeColor = color; 
p = plot([1 1], ci_fat_ss); p(1).Color = color/2;
p = plot([12 12], ci_fat_recog); p(1).Color = color/2;

b = bar([2 13], [beta_slf1and2_ss beta_slf1and2_recog], 0.075);
color = slf1and2;
b.FaceColor = color; b.EdgeColor = color; 
p = plot([2 2], ci_slf1and2_ss); p(1).Color = color/2;
p = plot([13 13], ci_slf1and2_recog); p(1).Color = color/2;

b = bar([3 14], [beta_slf3_ss beta_slf3_recog], 0.075);
color = slf3;
b.FaceColor = color; b.EdgeColor = color; 
p = plot([3 3], ci_slf3_ss); p(1).Color = color/2;
p = plot([14 14], ci_slf3_recog); p(1).Color = color/2;

b = bar([4 15], [beta_parc_ss beta_parc_recog], 0.075);
color = parc;
b.FaceColor = color; b.EdgeColor = color; 
p = plot([4 4], ci_parc_ss); p(1).Color = color/2;
p = plot([15 15], ci_parc_recog); p(1).Color = color/2;

b = bar([5 16], [beta_tpc_ss beta_tpc_recog], 0.075);
color = tpc;
b.FaceColor = color; b.EdgeColor = color; 
p = plot([5 5], ci_tpc_ss); p(1).Color = color/2;
p = plot([16 16], ci_tpc_recog); p(1).Color = color/2;

b = bar([6 17], [beta_mdlfang_ss beta_mdlfang_recog], 0.075);
color = mdlfang;
b.FaceColor = color; b.EdgeColor = color; 
p = plot([6 6], ci_mdlfang_ss); p(1).Color = color/2;
p = plot([17 17], ci_mdlfang_recog); p(1).Color = color/2;

b = bar([7 18], [beta_mdlfspl_ss beta_mdlfspl_recog], 0.075);
color = mdlfspl;
b.FaceColor = color; b.EdgeColor = color; 
p = plot([7 7], ci_mdlfspl_ss); p(1).Color = color/2;
p = plot([18 18], ci_mdlfspl_recog); p(1).Color = color/2;

b = bar([8 19], [beta_ilf_ss beta_ilf_recog], 0.075);
color = ilf;
b.FaceColor = color; b.EdgeColor = color; 
p = plot([8 8], ci_ilf_ss); p(1).Color = color/2;
p = plot([19 19], ci_ilf_recog); p(1).Color = color/2;

b = bar([9 20], [beta_ifof_ss beta_ifof_recog], 0.075);
color = ifof;
b.FaceColor = color; b.EdgeColor = color; 
p = plot([9 9], ci_ifof_ss); p(1).Color = color/2;
p = plot([20 20], ci_ifof_recog); p(1).Color = color/2;

b = bar([10 21], [beta_vof_ss beta_vof_recog], 0.075);
color = vof;
b.FaceColor = color; b.EdgeColor = color; 
p = plot([10 10], ci_vof_ss); p(1).Color = color/2;
p = plot([21 21], ci_vof_recog); p(1).Color = color/2;
xlim([0.5 17.5]);

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0.5 21.5];
xax.TickValues = [5.5 15.5];
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
yax.TickValues = [-.5 0 .5 1];
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
legend('location', 'eastoutside');
if strcmp(hemisphere, 'left')
legend({'left-fat', '', '', 'left-slf1and2', '', '', 'left-slf3', '', '', 'left-parc', '', '', 'left-tpc', '', '', ...
    'left-mdlfang', '', '', 'left-mdlfspl', '', '', 'left-ilf', '', '', 'left-ifof', '', '', 'left-vof'})
elseif strcmp(hemisphere, 'right')
    legend({'right-fat', '', '', 'right-slf1and2', '', '', 'right-slf3', '', '', 'right-parc', '', '', 'right-tpc', '', '', ...
    'right-mdlfang', '', '', 'right-mdlfspl', '', '', 'right-ilf', '', '', 'right-ifof', '', '', 'right-vof'})
elseif strcmp(hemisphere, 'both')
        legend({'fat', '', '', 'slf1and2', '', '', 'slf3', '', '', 'parc', '', '', 'tpc', '', '', ...
    'mdlfang', '', '', 'mdlfspl', '', '', 'ilf', '', '', 'ifof', '', '', 'vof'})
end
title(wmmeasure)

n = size(t, 1);
% text(-1.8, 1.8, ['adjr2 = ' num2str(f2.adjrsquare)])
% text(-1.8, 1.6, ['rmse = ' num2str(f2.rmse)])
% text(-1.8, 1.4, ['beta = ' num2str(f.p1)])
% text(-1.8, 1.2, ['p = ' num2str(mdl.Coefficients.pValue(2))])
% text(-1.8, 1.0, ['n = ' num2str(n)])

print(fullfile(rootDir, 'wml-wmpredictslearning-plots', ['plot_pred_learning_n=' num2str(n) '_' wmmeasure '_' hemisphere]), '-dpng')
print(fullfile(rootDir, 'wml-wmpredictslearning-plots', 'eps', ['plot_pred_learning_n=' num2str(n) '_' wmmeasure '_' hemisphere]), '-depsc')



