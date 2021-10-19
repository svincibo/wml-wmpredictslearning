clear all; clc;

% Set working directories.
rootDir = '/Volumes/Seagate/wml/wml-wmpredictslearning';

% Identify outliers for removal.
remove = [40 47 56 60 62];

% Load recog data and remove outliers and sort by subID.
datestring = '20211012';
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
datestring = '20211012';
filename = sprintf('WML_beh_data_write_%s', datestring);
load(fullfile(rootDir, 'wml-wmpredictslearning-supportFiles', filename), 'data_write', 'data_write_mean');

keep = find(~ismember(data_write_mean.Var1, remove));

temp = table2array(data_write_mean(keep, :));
mt = array2table(temp, 'VariableNames', {'subID', 'mt_day1', 'mt_day2', 'mt_day3', 'mt_day4', 'mt_day5'}); clear temp;
mt = sortrows(mt);

% Load writing data and remove outliers and sort by subID.
datestring = '20211012';
filename = sprintf('WML_mri_data_tractprofiles');%_%s', datestring);
load(fullfile(rootDir, 'wml-wmpredictslearning-supportFiles', filename), 'data_tractprofiles_mean');

keep = find(~ismember(data_tractprofiles_mean.subID, remove));

temp = table2array(data_tractprofiles_mean(keep, :));
mri = array2table(temp, 'VariableNames', data_tractprofiles_mean.Properties.VariableNames); clear temp;
mri = sortrows(mri);
fa = table2array(mri(:, 2:7:end)); figure(1); histogram(fa(:), 100); xlim([0, 1]); title('fa'); hold off;
md = table2array(mri(:, 3:7:end)); figure(2); histogram(md(:), 100); xlim([min(md(:)), max(md(:))]); title('md'); hold off;
t1t2 = table2array(mri(:, 4:7:end)); figure(3); histogram(t1t2(:), 100); xlim([min(t1t2(:)), max(t1t2(:))]); title('t1t2'); hold off; 
ndi = table2array(mri(:, 5:7:end)); figure(4); histogram(ndi(:), 100); xlim([min(ndi(:)), max(ndi(:))]); title('ndi'); hold off;
odi = table2array(mri(:, 6:7:end)); figure(5); histogram(odi(:), 100); xlim([min(odi(:)), max(odi(:))]); title('odi'); hold off;
isovf = table2array(mri(:, 7:7:end)); figure(6); histogram(isovf(:), 100); xlim([min(isovf(:)), max(isovf(:))]); title('isovf'); hold off;

% Concatenate all variables into a table.
t = [mt rt(:, 2:end) acc(:, 2:end) mri(:, 2:end)];

% Delete all rows that contain NaN, for now. Columns first.
idxc = [6 11 16];
t(:, idxc) = [];
idxr = [2 3 7 8 12 14 15 16];
t(idxr, :) = [];





