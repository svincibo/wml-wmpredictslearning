clear all; close all; clc;

% Identify outliers for removal.
remove = [27 32]; 
% MRI
% 32, 33, 35, 41 have spike artifacts
% 27 has severe motion
% 24, 31 withdrew mid-training, so multi-day learningrate is not accurate
% (but day 1 learning rate is accurate) 
% 22 missing left CST
% 3, NaN for rightmdlfspl, rightmdlfang, rightvof ~ sub 25
% 7, NaN for rightmdlfspl ~ sub 31
% 8, NaN for leftifof ~ sub 32
% 12, rightifof, leftuncinate ~ sub 42

% Training
% 24 have only day 1 data
% 31 have only day 1 and day 2 data
% 32 had 0% accuracy by day 4 
% 50 missing day 3, experimenter error
% 66 in progress, as of 10/26/21

% Set working directories.
rootdir = '/Volumes/Seagate/wml/wml-wmpredictslearning';

% Find all subject folders in the tractprofiles directory (generated with wml_datacat_tractprofiles_pca.m).
subfolders = dir(fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', 'tractprofiles'));

% Remove the '.' and '..' files.
subfolders = subfolders(arrayfun(@(x) x.name(1), subfolders) ~= '.');

% Keep only names that are csv files.
subfolders = subfolders(arrayfun(@(x) x.name(end), subfolders) == 'v');

% Perform PCA on each subjects white matter data.
for s = 2:size(subfolders, 1)
        
    % Grab demographics data.
    subID = str2num(subfolders(s).name(end-5:end-4));
    
    % Read in white matter data for this subject.
    d = readtable(fullfile(subfolders(s).folder, subfolders(s).name));
        
    % Prepare data for PCA by converting to matrix and removing NaN columns.
    dmat = table2array(d);
    keep = find(all(~isnan(dmat)));
    d = d(:, keep);
    dmat = dmat(:, keep);
    
    % Get tract names.
    tractnames = d.Properties.VariableNames;
    
    % Check that the difference in variance of different columns is not very large.
    bar(var(dmat));
    
    % PCA: Most of the time it seems large, so will default to scaling the data (i.e, weighted).
    w = 1./var(dmat);
%     [wcoeff, score, latent, tsquared, explained] = pca(dmat, 'NumComponents', 10, 'VariableWeights', w, 'Rows', 'pairwise', 'Economy', true);
        [wcoeff, score, latent, tsquared, explained] = pca(dmat, 'VariableWeights', w, 'Rows', 'pairwise', 'Economy', true);

    % Scree plot: latent contains the eigenvalues
    figure(1)
    plot(latent(1:10), 'o', 'LineStyle', '-', 'Color', 'b');
    ylabel('eigenvalue')
    xlabel('principle component')

    % Just use first three principle components for now.
    
    % Compute coefficients.
    c3 = wcoeff(:, 1:3);
    
    % Transform coefficients so that they are orthonormal.
    coefforth = diag(sqrt(w))*wcoeff;
    
%     % Check that they are, indeed, orthonormal.
%     I = coefforth'*coefforth;
%     I(1:3,1:3)

% Component scores
cscores = zscore(dmat)*coefforth;


figure(2)
pareto(explained)
xlabel('Principal Component')
ylabel('Variance Explained (%)')

figure(3)
plot(cscores(:,1),cscores(:,2),'+')
xlabel('1st Principal Component')
ylabel('2nd Principal Component')   

figure(4)
plot3(cscores(:,1),cscores(:,2),cscores(:, 3),'+')
xlabel('1st Principal Component')
ylabel('2nd Principal Component')   
zlabel('3rd Principal Component')   
    
end



% Load recog data and remove outliers and sort by subID.
datestring = '20211026';
n = 29;
fittype = 'poly1';

% Load recog acc data and remove outliers and sort by subID.
filename = sprintf('WML_beh_data_recog_test_%s', datestring);
load(fullfile(rootDir, 'wml-wmpredictslearning-supportFiles', filename), 'data_recog', 'data_recog_acc_mean');

keep = find(~ismember(data_recog_acc_mean.Var1, remove));

temp = table2array(data_recog_acc_mean(keep, :));
acc = array2table(temp, 'VariableNames', {'subID', 'acc_day1', 'acc_day2', 'acc_day3', 'acc_day4', 'acc_day5'}); clear temp;
acc = sortrows(acc);

% Load writing data and remove outliers and sort by subID.
filename = sprintf('WML_beh_data_write_wlearningrate_%s_n=%d_%s', fittype, n, datestring);
load(fullfile(rootDir, 'wml-wmpredictslearning-supportFiles', filename), 'data_write', 'data_write_mean');

keep = find(~ismember(data_write_mean.Var1, remove));

temp = table2array(data_write_mean(keep, :));
mt = array2table(temp, 'VariableNames', {'subID', 'mt_day1', 'mt_day2', 'mt_day3', 'mt_day4', 'mt_day5', 'mt_learningrate'}); clear temp;
mt = sortrows(mt, 1);
mt.mt_learningrate = -mt.mt_learningrate; % make faster learners positive (easier on my brain)

% Load recog rt data with learning rate and remove outliers and sort by subID.
filename = sprintf('WML_beh_data_recog_rt_wlearningrate_%s_n=%d_%s', fittype, n, datestring);
load(fullfile(rootDir, 'wml-wmpredictslearning-supportFiles', filename), 'data_recog', 'data_recog_rt_mean');

keep = find(~ismember(data_recog_rt_mean.Var1, remove));

temp = table2array(data_recog_rt_mean(keep, :));
rt = array2table(temp, 'VariableNames', {'subID', 'rt_day1', 'rt_day2', 'rt_day3', 'rt_day4', 'rt_day5', 'rt_learningrate'}); clear temp;
rt = sortrows(rt, 1);
rt.rt_learningrate = -rt.rt_learningrate; % make faster learners positive (easier on my brain)

% Load tractprofiles data and remove outliers and sort by subID.
filename = sprintf('WML_mri_data_tractprofiles');%_%s', datestring);
load(fullfile(rootDir, 'wml-wmpredictslearning-supportFiles', filename), 'data_tractprofiles_mean');

keep = find(~ismember(data_tractprofiles_mean.subID, remove));

temp = table2array(data_tractprofiles_mean(keep, :));
mri = array2table(temp, 'VariableNames', data_tractprofiles_mean.Properties.VariableNames); clear temp;
mri = sortrows(mri);

% % QA: make histograms of each mri measure.
% fa = table2array(mri(:, 2:6:end)); figure(1); histogram(fa(:), 100); xlim([0, 1]); title('fa'); hold off;
% md = table2array(mri(:, 3:6:end)); figure(2); histogram(md(:), 100); xlim([min(md(:)), max(md(:))]); title('md'); hold off;
% t1t2 = table2array(mri(:, 4:6:end)); figure(3); histogram(t1t2(:), 100); xlim([min(t1t2(:)), max(t1t2(:))]); title('t1t2'); hold off; 
% ndi = table2array(mri(:, 5:6:end)); figure(4); histogram(ndi(:), 100); xlim([min(ndi(:)), max(ndi(:))]); title('ndi'); hold off;
% odi = table2array(mri(:, 6:6:end)); figure(5); histogram(odi(:), 100); xlim([min(odi(:)), max(odi(:))]); title('odi'); hold off;
% isovf = table2array(mri(:, 7:6:end)); figure(6); histogram(isovf(:), 100); xlim([min(isovf(:)), max(isovf(:))]); title('isovf'); hold off;

% Concatenate all variables into a table and zscore each column separately.
t = [mt rt(:, 2:end) acc(:, 2:end) mri(:, 2:end)];
t_temp = cat(2, table2array(t(:, 1)), (table2array(t(:, 2:end))-nanmean(table2array(t(:, 2:end)), 1))./nanstd(table2array(t(:, 2:end)), 1));
t = array2table(t_temp, 'VariableNames', t.Properties.VariableNames);

% Delete all columns and rows that contain NaN, for now. Columns first.
idxc = [1 6 12 18]; % 1 is subID, 6, 12, and 18 correspond to day 5 test measurements
t(:, idxc) = [];
idxr = []; 
t(idxr, :) = [];

% Correlations.
% figure(7);
% mat = corr(table2array(t(:, 2:end)));% % m = cat(2, table2array(mt(:, 2:end)), table2array(mri(:, 2:end)));
% toplot = mat(1:12, :);
% imagesc(toplot); colorbar; caxis([0.4 .9]);
% xlabels = t.Properties.VariableNames(2:end);
% ylabels = t.Properties.VariableNames(2:13);

mt_lr_idx = 5;
rt_lr_idx = 10;
mt1_idx = 1;
rt1_idx = 6;
acc1_idx = 11;
first_fa_idx = 15; %15=fa, 16=md
% first_t1t2_idx = 17;
% first_ndi_idx = 18;
step = 6;
corr_thresh = .3;
p_thresh = .05;

% FA
figure(1);
[mat, p] = corr(table2array(t(:, [mt_lr_idx, mt1_idx, rt_lr_idx, rt1_idx, acc1_idx, first_fa_idx:step:end])));% % m = cat(2, table2array(mt(:, 2:end)), table2array(mri(:, 2:end)));
% id=find(abs(mat(:))<corr_thresh); mat(id) = NaN;
id=find(p>p_thresh); mat(id) = NaN;
toplot = mat;
imagesc(toplot); colorbar; caxis([-1 1]);
xlabels = t.Properties.VariableNames([mt_lr_idx, mt1_idx, rt_lr_idx, rt1_idx, acc1_idx, first_fa_idx:step:end]);
ylabels = t.Properties.VariableNames([mt_lr_idx, mt1_idx, rt_lr_idx, rt1_idx, acc1_idx, first_fa_idx:step:end]);

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

hold off;

% 
% % mat = corr(m);
% % imagesc(mat);
% % colorbar;
% 
% % t1t2
% figure(8);
% mat = corr(table2array(t(:, [lr_idx, mt1_idx, rt1_idx, acc1_idx, first_t1t2_idx:step:end])));% % m = cat(2, table2array(mt(:, 2:end)), table2array(mri(:, 2:end)));
% id=find(abs(mat(:))<corr_thresh); mat(id) = NaN;
% toplot = mat;
% imagesc(toplot); colorbar; caxis([-1 1]);
% xlabels = t.Properties.VariableNames([lr_idx, mt1_idx, rt1_idx, acc1_idx, first_t1t2_idx:step:end]);
% ylabels = t.Properties.VariableNames([lr_idx, mt1_idx, rt1_idx, acc1_idx, first_t1t2_idx:step:end]);
% 
% % Set up plot and measure-specific details.
% capsize = 0;
% marker = 'o';
% linewidth = 1.5;
% linestyle = 'none';
% markersize = 100;
% xtickvalues = [1 2 3 4];
% xlim_lo = min(xtickvalues)-0.5; xlim_hi = max(xtickvalues)+0.5;
% fontname = 'Arial';
% fontsize = 14;
% fontangle = 'italic';
% yticklength = 0;
% xticklength = 0.05;
% 
% % xaxis
% xax = get(gca, 'xaxis');
% % xax.Limits = [xlim_lo xlim_hi];
% xax.TickValues = 1:size(toplot, 2);
% xax.TickDirection = 'out';
% % xax.TickLength = [xticklength xticklength];
% % xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
% xax.FontName = fontname;
% xax.FontSize = fontsize;
% % xax.FontAngle = fontangle;
% xax.TickLabelRotation = 90;
% 
% % yaxis
% yax = get(gca,'yaxis');
% % yax.Limits = [ylimlo ylimhi];
% yax.TickValues = 1:size(toplot, 1);
% % yax.TickDirection = 'out';
% % yax.TickLength = [yticklength yticklength];
% yax.TickLabels = ylabels;
% yax.FontName = fontname;
% yax.FontSize = fontsize;
% yax.TickLabelRotation = 0;
% 
% title('t1/t2 ratio')
% 
% 
% % mat = corr(m);
% % imagesc(mat);
% % colorbar;
% 
% % ndi
% figure(9);
% mat = corr(table2array(t(:, [lr_idx, mt1_idx, rt1_idx, acc1_idx, first_ndi_idx:step:end])));% % m = cat(2, table2array(mt(:, 2:end)), table2array(mri(:, 2:end)));
% id=find(abs(mat(:))<corr_thresh); mat(id) = NaN;
% toplot = mat;
% imagesc(toplot); colorbar; caxis([-1 1]);
% xlabels = t.Properties.VariableNames([lr_idx, mt1_idx, rt1_idx, acc1_idx, first_ndi_idx:step:end]);
% ylabels = t.Properties.VariableNames([lr_idx, mt1_idx, rt1_idx, acc1_idx, first_ndi_idx:step:end]);
% 
% % Set up plot and measure-specific details.
% capsize = 0;
% marker = 'o';
% linewidth = 1.5;
% linestyle = 'none';
% markersize = 100;
% xtickvalues = [1 2 3 4];
% xlim_lo = min(xtickvalues)-0.5; xlim_hi = max(xtickvalues)+0.5;
% fontname = 'Arial';
% fontsize = 14;
% fontangle = 'italic';
% yticklength = 0;
% xticklength = 0.05;
% 
% % xaxis
% xax = get(gca, 'xaxis');
% % xax.Limits = [xlim_lo xlim_hi];
% xax.TickValues = 1:size(toplot, 2);
% xax.TickDirection = 'out';
% % xax.TickLength = [xticklength xticklength];
% % xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
% xax.FontName = fontname;
% xax.FontSize = fontsize;
% % xax.FontAngle = fontangle;
% xax.TickLabelRotation = 90;
% 
% % yaxis
% yax = get(gca,'yaxis');
% % yax.Limits = [ylimlo ylimhi];
% yax.TickValues = 1:size(toplot, 1);
% % yax.TickDirection = 'out';
% % yax.TickLength = [yticklength yticklength];
% yax.TickLabels = ylabels;
% yax.FontName = fontname;
% yax.FontSize = fontsize;
% yax.TickLabelRotation = 0;
% 
% title('NDI')


% mat = corr(m);
% imagesc(mat);
% colorbar;

figure(2);

% Fit a linear model: Does wm predict mt at day 1?
modelspec = 'mt_learningrate~leftparc_fa';
mdl = fitlm(t, modelspec);
mdl.Coefficients
anova(mdl, 'summary')

% Another way, easier to plot with: Extract data for this fit.
y = t.mt_learningrate;
x = t.leftparc_fa;

% Fit.
[f, f2] = fit(x, y, 'poly1'); %, 'Robust', 'Lar');

p=plot(f, x, y);
    hold on;

        % Customize scatter plot of data.
    p(1).Color = [0 128 128]/255; %teal
    p(1).MarkerSize = 40;
    % Customize line plot of data.
    p(2).Color = [0 128 128]/255; %teal
    p(2).LineWidth = linewidth;
    p(2).LineStyle = '--';
    
fontsize = 20;

xlim_lo = -2; xlim_hi = 2;
ylim_lo = -2; ylim_hi = 2;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [-2 -1 0 1 2];
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
yax.TickValues = [-2 -1 0 1 2];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.TickLabels = ylabels;
yax.FontName = fontname;
yax.FontSize = fontsize;
% yax.TickLabelRotation = 0;

% title('Left pArc')
a = gca;
a.YLabel.String = {'Learning Rate of Sensorimotor Skill'; '(z-scored)'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Fractional Anisotropy, Left pArc'; '(z-scored)'};
a.XLabel.FontSize = fontsize;

box off;
legend('box', 'off');
legend('location', 'southeast');
pbaspect([1 1 1])

n = size(x, 1);
text(-1.8, 1.8, ['adjr2 = ' num2str(f2.adjrsquare)])
text(-1.8, 1.6, ['rmse = ' num2str(f2.rmse)])
text(-1.8, 1.4, ['beta = ' num2str(f.p1)])
text(-1.8, 1.2, ['p = ' num2str(mdl.Coefficients.pValue(2))])
text(-1.8, 1.0, ['n = ' num2str(n)])

print(fullfile(rootDir, 'wml-wmpredictslearning-plots', ['plot_corr_mt_leftparc_n=' num2str(n)]), '-dpng')
print(fullfile(rootDir, 'wml-wmpredictslearning-plots', 'eps', ['plot_corr_mt_leftparc_n=' num2str(n)]), '-depsc')

hold off;

figure(3)

% Fit a linear model: Does wm predict mt at day 1?
modelspec = 'rt_learningrate~leftparc_fa';
mdl = fitlm(t, modelspec);
mdl.Coefficients
anova(mdl, 'summary')

% Another way, easier to plot with: Extract data for this fit.
y = t.rt_learningrate;
x = t.leftparc_fa;

% Fit.
[f, f2] = fit(x, y, 'poly1'); %, 'Robust', 'Lar');

p=plot(f, x, y);
    hold on;

        % Customize scatter plot of data.
    p(1).Color = [0 128 128]/255; %teal
    p(1).MarkerSize = 40;
    % Customize line plot of data.
    p(2).Color = [0 128 128]/255; %teal
    p(2).LineWidth = linewidth;
    p(2).LineStyle = '--';
    
fontsize = 20;

xlim_lo = -2; xlim_hi = 2;
ylim_lo = -2; ylim_hi = 2;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [-2 -1 0 1 2];
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
yax.TickValues = [-2 -1 0 1 2];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.TickLabels = ylabels;
yax.FontName = fontname;
yax.FontSize = fontsize;
% yax.TickLabelRotation = 0;

% title('Left pArc')
a = gca;
a.YLabel.String = {'Learning Rate of Visual Recognition'; '(z-scored)'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Fractional Anisotropy, Left pArc'; '(z-scored)'};
a.XLabel.FontSize = fontsize;

box off;
legend('box', 'off');
legend('location', 'southeast');
pbaspect([1 1 1])

n = size(x, 1);
text(-1.8, 1.8, ['adjr2 = ' num2str(f2.adjrsquare)])
text(-1.8, 1.6, ['rmse = ' num2str(f2.rmse)])
text(-1.8, 1.4, ['beta = ' num2str(f.p1)])
text(-1.8, 1.2, ['p = ' num2str(mdl.Coefficients.pValue(2))])
text(-1.8, 1.0, ['n = ' num2str(n)])

print(fullfile(rootDir, 'wml-wmpredictslearning-plots', ['plot_corr_rt_leftparc_n=' num2str(n)]), '-dpng')
print(fullfile(rootDir, 'wml-wmpredictslearning-plots', 'eps', ['plot_corr_rt_leftparc_n=' num2str(n)]), '-depsc')

figure(4);

% Fit a linear model: Does wm predict mt at day 1?
modelspec = 'mt_learningrate~leftslf1and2_fa';
mdl = fitlm(t, modelspec);
mdl.Coefficients
anova(mdl, 'summary')

% Another way, easier to plot with: Extract data for this fit.
y = t.mt_learningrate;
x = t.leftslf1and2_fa;

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

xlim_lo = -2; xlim_hi = 2;
ylim_lo = -2; ylim_hi = 2;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [-2 -1 0 1 2];
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
yax.TickValues = [-2 -1 0 1 2];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.TickLabels = ylabels;
yax.FontName = fontname;
yax.FontSize = fontsize;
% yax.TickLabelRotation = 0;

% title('Left SLF1/2')
a = gca;
a.YLabel.String = {'Learning Rate of Sensorimotor Skill'; '(z-scored)'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Fractional Anisotropy, Left SLF1/2'; '(z-scored)'};
a.XLabel.FontSize = fontsize;

box off;
legend('box', 'off');
legend('location', 'southeast');
pbaspect([1 1 1])

n = size(x, 1);
text(-1.8, 1.8, ['adjr2 = ' num2str(f2.adjrsquare)])
text(-1.8, 1.6, ['rmse = ' num2str(f2.rmse)])
text(-1.8, 1.4, ['beta = ' num2str(f.p1)])
text(-1.8, 1.2, ['p = ' num2str(mdl.Coefficients.pValue(2))])
text(-1.8, 1.0, ['n = ' num2str(n)])

print(fullfile(rootDir, 'wml-wmpredictslearning-plots', ['plot_corr_mt_leftslf1and2_n=' num2str(n)]), '-dpng')
print(fullfile(rootDir, 'wml-wmpredictslearning-plots', 'eps', ['plot_corr_mt_leftslf1and2_n=' num2str(n)]), '-depsc')

hold off;

figure(5);

% Fit a linear model: Does wm predict mt at day 1?
modelspec = 'rt_learningrate~leftslf1and2_fa';
mdl = fitlm(t, modelspec);
mdl.Coefficients
anova(mdl, 'summary')

% Another way, easier to plot with: Extract data for this fit.
y = t.rt_learningrate;
x = t.leftslf1and2_fa;

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

xlim_lo = -2; xlim_hi = 2;
ylim_lo = -2; ylim_hi = 2;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [-2 -1 0 1 2];
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
yax.TickValues = [-2 -1 0 1 2];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.TickLabels = ylabels;
yax.FontName = fontname;
yax.FontSize = fontsize;
% yax.TickLabelRotation = 0;

% title('Left SLF1/2')
a = gca;
a.YLabel.String = {'Learning Rate of Visual Recognition'; '(z-scored)'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Fractional Anisotropy, Left SLF1/2'; '(z-scored)'};
a.XLabel.FontSize = fontsize;

box off;
legend('box', 'off');
legend('location', 'southeast');
pbaspect([1 1 1])

n = size(x, 1);
text(-1.8, 1.8, ['adjr2 = ' num2str(f2.adjrsquare)])
text(-1.8, 1.6, ['rmse = ' num2str(f2.rmse)])
text(-1.8, 1.4, ['beta = ' num2str(f.p1)])
text(-1.8, 1.2, ['p = ' num2str(mdl.Coefficients.pValue(2))])
text(-1.8, 1.0, ['n = ' num2str(n)])

print(fullfile(rootDir, 'wml-wmpredictslearning-plots', ['plot_corr_rt_leftslf1and2_n=' num2str(n)]), '-dpng')
print(fullfile(rootDir, 'wml-wmpredictslearning-plots', 'eps', ['plot_corr_rt_leftslf1and2_n=' num2str(n)]), '-depsc')

hold off;

figure(6);

% Fit a linear model: Does wm predict mt at day 1?
modelspec = 'mt_learningrate~leftifof_fa';
mdl = fitlm(t, modelspec);
mdl.Coefficients
anova(mdl, 'summary')

% Another way, easier to plot with: Extract data for this fit.
y = t.mt_learningrate;
x = t.leftifof_fa;

% Fit.
[f, f2] = fit(x, y, 'poly1'); %, 'Robust', 'Lar');

p=plot(f, x, y);
    hold on;

        % Customize scatter plot of data.
    p(1).Color = [30 144 250]/255; %blue
    p(1).MarkerSize = 40;
    % Customize line plot of data.
    p(2).Color = [30 144 250]/255; %blue
    p(2).LineWidth = linewidth;
    p(2).LineStyle = '--';
    
fontsize = 20;

xlim_lo = -2; xlim_hi = 2;
ylim_lo = -2; ylim_hi = 2;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [-2 -1 0 1 2];
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
yax.TickValues = [-2 -1 0 1 2];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.TickLabels = ylabels;
yax.FontName = fontname;
yax.FontSize = fontsize;
% yax.TickLabelRotation = 0;

% title('Left IFOF')
a = gca;
a.YLabel.String = {'Learning Rate of Sensorimotor Skill'; '(z-scored)'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Fractional Anisotropy, Left IFOF'; '(z-scored)'};
a.XLabel.FontSize = fontsize;

box off;
legend('box', 'off');
legend('location', 'southeast');
pbaspect([1 1 1])

n = size(x, 1);
text(-1.8, 1.8, ['adjr2 = ' num2str(f2.adjrsquare)])
text(-1.8, 1.6, ['rmse = ' num2str(f2.rmse)])
text(-1.8, 1.4, ['beta = ' num2str(f.p1)])
text(-1.8, 1.2, ['p = ' num2str(mdl.Coefficients.pValue(2))])
text(-1.8, 1.0, ['n = ' num2str(n)])

print(fullfile(rootDir, 'wml-wmpredictslearning-plots', ['plot_corr_mt_leftifof_n=' num2str(n)]), '-dpng')
print(fullfile(rootDir, 'wml-wmpredictslearning-plots', 'eps', ['plot_corr_mt_leftifof_n=' num2str(n)]), '-depsc')

hold off;

figure(7);

% Fit a linear model: Does wm predict mt at day 1?
modelspec = 'rt_learningrate~leftifof_fa';
mdl = fitlm(t, modelspec);
mdl.Coefficients
anova(mdl, 'summary')

% Another way, easier to plot with: Extract data for this fit.
y = t.rt_learningrate;
x = t.leftifof_fa;

% Fit.
[f, f2] = fit(x, y, 'poly1'); %, 'Robust', 'Lar');

p=plot(f, x, y);
    hold on;

        % Customize scatter plot of data.
    p(1).Color = [30 144 250]/255; %blue
    p(1).MarkerSize = 40;
    % Customize line plot of data.
    p(2).Color = [30 144 250]/255; %blue
    p(2).LineWidth = linewidth;
    p(2).LineStyle = '--';
    
fontsize = 20;

xlim_lo = -2; xlim_hi = 2;
ylim_lo = -2; ylim_hi = 2;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [-2 -1 0 1 2];
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
yax.TickValues = [-2 -1 0 1 2];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.TickLabels = ylabels;
yax.FontName = fontname;
yax.FontSize = fontsize;
% yax.TickLabelRotation = 0;

% title('Left IFOF')
a = gca;
a.YLabel.String = {'Learning Rate of Visual Recognition'; '(z-scored)'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Fractional Anisotropy, Left IFOF'; '(z-scored)'};
a.XLabel.FontSize = fontsize;

box off;
legend('box', 'off');
legend('location', 'southeast');
pbaspect([1 1 1])

n = size(x, 1);
text(-1.8, 1.8, ['adjr2 = ' num2str(f2.adjrsquare)])
text(-1.8, 1.6, ['rmse = ' num2str(f2.rmse)])
text(-1.8, 1.4, ['beta = ' num2str(f.p1)])
text(-1.8, 1.2, ['p = ' num2str(mdl.Coefficients.pValue(2))])
text(-1.8, 1.0, ['n = ' num2str(n)])

print(fullfile(rootDir, 'wml-wmpredictslearning-plots', ['plot_corr_rt_leftifof_n=' num2str(n)]), '-dpng')
print(fullfile(rootDir, 'wml-wmpredictslearning-plots', 'eps', ['plot_corr_rt_leftifof_n=' num2str(n)]), '-depsc')

hold off;

figure(6);

% Fit a linear model: Does wm predict mt at day 1?
modelspec = 'mt_learningrate~leftilf_fa';
mdl = fitlm(t, modelspec);
mdl.Coefficients
anova(mdl, 'summary')

% Another way, easier to plot with: Extract data for this fit.
y = t.mt_learningrate;
x = t.leftilf_fa;

% Fit.
[f, f2] = fit(x, y, 'poly1'); %, 'Robust', 'Lar');

p=plot(f, x, y);
    hold on;

        % Customize scatter plot of data.
    p(1).Color = [30 144 250]/255; %blue
    p(1).MarkerSize = 40;
    % Customize line plot of data.
    p(2).Color = [30 144 250]/255; %blue
    p(2).LineWidth = linewidth;
    p(2).LineStyle = '--';
    
fontsize = 20;

xlim_lo = -2; xlim_hi = 2;
ylim_lo = -2; ylim_hi = 2;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [-2 -1 0 1 2];
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
yax.TickValues = [-2 -1 0 1 2];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.TickLabels = ylabels;
yax.FontName = fontname;
yax.FontSize = fontsize;
% yax.TickLabelRotation = 0;

% title('Left ILF')
a = gca;
a.YLabel.String = {'Learning Rate of Sensorimotor Skill'; '(z-scored)'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Fractional Anisotropy, Left ILF'; '(z-scored)'};
a.XLabel.FontSize = fontsize;

box off;
legend('box', 'off');
legend('location', 'southeast');
pbaspect([1 1 1])

n = size(x, 1);
text(-1.8, 1.8, ['adjr2 = ' num2str(f2.adjrsquare)])
text(-1.8, 1.6, ['rmse = ' num2str(f2.rmse)])
text(-1.8, 1.4, ['beta = ' num2str(f.p1)])
text(-1.8, 1.2, ['p = ' num2str(mdl.Coefficients.pValue(2))])
text(-1.8, 1.0, ['n = ' num2str(n)])

print(fullfile(rootDir, 'wml-wmpredictslearning-plots', ['plot_corr_mt_leftilf_n=' num2str(n)]), '-dpng')
print(fullfile(rootDir, 'wml-wmpredictslearning-plots', 'eps', ['plot_corr_mt_leftilf_n=' num2str(n)]), '-depsc')

hold off;

figure(7);

% Fit a linear model: Does wm predict mt at day 1?
modelspec = 'rt_learningrate~leftilf_fa';
mdl = fitlm(t, modelspec);
mdl.Coefficients
anova(mdl, 'summary')

% Another way, easier to plot with: Extract data for this fit.
y = t.rt_learningrate;
x = t.leftilf_fa;

% Fit.
[f, f2] = fit(x, y, 'poly1'); %, 'Robust', 'Lar');

p=plot(f, x, y);
    hold on;

        % Customize scatter plot of data.
    p(1).Color = [30 144 250]/255; %blue
    p(1).MarkerSize = 40;
    % Customize line plot of data.
    p(2).Color = [30 144 250]/255; %blue
    p(2).LineWidth = linewidth;
    p(2).LineStyle = '--';
    
fontsize = 20;

xlim_lo = -2; xlim_hi = 2;
ylim_lo = -2; ylim_hi = 2;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [-2 -1 0 1 2];
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
yax.TickValues = [-2 -1 0 1 2];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.TickLabels = ylabels;
yax.FontName = fontname;
yax.FontSize = fontsize;
% yax.TickLabelRotation = 0;

% title('Left ILF')
a = gca;
a.YLabel.String = {'Learning Rate of Visual Recognition'; '(z-scored)'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Fractional Anisotropy, Left ILF'; '(z-scored)'};
a.XLabel.FontSize = fontsize;

box off;
legend('box', 'off');
legend('location', 'southeast');
pbaspect([1 1 1])

n = size(x, 1);
text(-1.8, 1.8, ['adjr2 = ' num2str(f2.adjrsquare)])
text(-1.8, 1.6, ['rmse = ' num2str(f2.rmse)])
text(-1.8, 1.4, ['beta = ' num2str(f.p1)])
text(-1.8, 1.2, ['p = ' num2str(mdl.Coefficients.pValue(2))])
text(-1.8, 1.0, ['n = ' num2str(n)])

print(fullfile(rootDir, 'wml-wmpredictslearning-plots', ['plot_corr_rt_leftilf_n=' num2str(n)]), '-dpng')
print(fullfile(rootDir, 'wml-wmpredictslearning-plots', 'eps', ['plot_corr_rt_leftilf_n=' num2str(n)]), '-depsc')

hold off;


