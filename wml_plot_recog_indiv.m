clear all; close all; clc
%%%%%%%%%%% NOTE YET WORKING %%%%%%%%%%%%

% Set working directories.
rootDir = '/Volumes/Seagate/wml/';
remove = []; 
% 24 have only day 1 data
% 31 have only day 1 and day 2 data
% 32 had 0% accuracy by day 4 
% 50 missing day 3, experimenter error
% 66 in progress, as of 10/26/21

testortestgen = 'test'; %'test' or 'test_gen'
fittype = 'poly1'; % exp2, poly1

% Create date-specific file name that indicates how many subjects.
datestring = '20211026';
filename = sprintf('wml_beh_data_recog_%s_%s', testortestgen, datestring);

% Load data.
load(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-supportFiles', filename), 'data_recog', 'data_recog_rt_mean', 'data_recog_acc_mean');

% Remove subject outliers.
keep = find(~ismember(data_recog.subID, remove));
data_recog = data_recog(keep, :);

% Sort by subID, then by day, then by trial.
recog = sortrows(data_recog, [1 7 3]);

% Also, keep only RTs for accurate responses -- trying out true positives for now.
recog = recog(find(recog.truepositive == 1), :);

% Remove day 5 for now.
recog = recog(find(recog.day ~= 5), :);

% % Use day 1 only for now.
% recog = recog(find(recog.day == 1), :);

% Add overall trial counter.
recog.idx = transpose(1:size(recog, 1));

% Remove outliers based on subject/day means and timeouts -- right now assumed timeout if RT = 1. 
idx_remove = [];
subjectlist = unique(recog.subID);
subjectlist = subjectlist(~isnan(subjectlist));
for s = 1:length(subjectlist)
    
    % Find idx for this subject.
    sidx = find(recog.subID == subjectlist(s));
    
    if ~isempty(sidx)
        
        % Add overal trial counter for this subject.
        recog.trial_all(sidx(1):sidx(end)) = transpose(1:length(sidx));
        
        for d = 1:4
            
            didx = find(recog.subID == subjectlist(s) & recog.day == d);
            
            datanow = recog.RT(didx, :);
            m = nanmean(datanow); sd = nanstd(datanow);
            
            % Get lower and upper bound.
            idx_above = find(datanow > m+3*sd | datanow >= 1);
            idx_below = find(datanow < m-3*sd);
            % Remove outliers.
            if ~isempty(idx_above) | ~isempty(idx_below)
                idx = didx(cat(1, idx_above, idx_below), :);
            end
            
            % Remove
            if exist('idx')
                
                idx_remove = cat(1, idx_remove, idx);
                clear idx;
                
            end
            
        end
        
    end
    
end
clear datanow d didx;
% Remove the outliers.
recog(idx_remove, :) = [];

% Remove NaNs.
recog = recog(~isnan(recog.RT), :);

%% Plot scatter plot of individual data -- with all outliers removed.

% Plot trial x rt.
gscatter(recog.trial_all, recog.RT, recog.subID)

% Set up plot and measure-specific details.
capsize = 0;
marker = 'o';
linewidth = 1.5;
linestyle = 'none';
markersize = 100;
xtickvalues = 1:160;
xlim_lo = min(xtickvalues); xlim_hi = max(xtickvalues);
fontname = 'Arial';
fontsize = 15;
fontangle = 'italic';
yticklength = 0;
xticklength = 0.05;
coloralpha = .05;

ylimlo = .4; ylimhi = 1.0;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [40 80 120 160];
% xax.TickDirection = 'out';
% xax.TickLength = [xticklength xticklength];
xlabels = {'40', '80', '120', '160'};
% % xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylimlo ylimhi];
yax.TickValues = [ylimlo (ylimlo+ylimhi)/2 ylimhi];
yax.TickDirection = 'out';
yax.TickLength = [yticklength yticklength];
yax.TickLabels = {num2str(ylimlo, '%2.2f'), num2str((ylimlo+ylimhi)/2, '%2.2f'), num2str(ylimhi, '%2.2f')};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.FontAngle = fontangle;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off
%
legend('off')

a.YLabel.String = 'Draw Duration Per Symbol (seconds)';
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = 'Trial Number';
a.XLabel.FontSize = fontsize;
pbaspect([1 1 1])
n = length(subjectlist);
print(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-plots', ['plot_recog_scatter_n=' num2str(n)]), '-dpng')
print(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-plots', 'eps', ['plot_recog_scatter_n=' num2str(n)]), '-depsc')

hold off;

%% Fit a double exponential to each person's MT data to get learning rate.

co = colormap('lines');
steps = floor(size(co, 1)/length(subjectlist));
c = co(1:steps:end, :);

% One subject at a time: Estimate learning rate by fitting a double exponential function,
% Kahn et al., 2017, Cerebral Cortex.

% % Just for now, do this only for day 1.
% didx = find(mt.day == 1);
% mt = mt(didx, :);

for s = 1:length(subjectlist)
    
    % Find idx for this subject.
    sidx = find(recog.subID == subjectlist(s));
        
    % Extract data for this fit.
    y = recog.RT(sidx, :);
    x = recog.trial_all(sidx, :);
    d = recog.day(sidx, :);
    
%     % Regress out effect of day.
%     [f, f2] = fit(d, y, 'poly1'); %, 'Robust', 'Lar');
%     y_resid = y - f.p1*d;
    
    % Fit.
    [f, f2] = fit(x, y, fittype); %, 'Robust', 'Lar');
    
    % Collect the learning rate.
    if strcmp(fittype, 'exp2')
        
        kappa = f.d;
        lambda = f.b;
        learningrate(s) = kappa;
        
    elseif strcmp(fittype, 'poly1')
        
        learningrate(s) = f.p1;
        
    end
 
    figure(s)
    
    % Plot fit.
    p=plot(f, x, y);
    hold on;
    
    % Customize scatter plot of data.
    p(1).Color = c(s, :);
    
    % Customize line plot of data.
    p(2).Color = c(s, :);
    p(2).LineWidth = 1;
    
    title(['Participant ' num2str(subjectlist(s)) ', ' fittype ', adjr2 = ' num2str(f2.adjrsquare) ', rmse = ' num2str(f2.rmse)]);
    
    % xaxis
    xax = get(gca, 'xaxis');
    xax.Limits = [xlim_lo xlim_hi];
    xax.TickValues = [40 80 120 160];
    % xax.TickDirection = 'out';
    % xax.TickLength = [xticklength xticklength];
    xlabels = {'40', '80', '120', '160'};
    % % xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
    xax.TickLabels = xlabels;
    xax.FontName = fontname;
    xax.FontSize = fontsize;
    
    % yaxis
    yax = get(gca,'yaxis');
    yax.Limits = [ylimlo ylimhi];
    yax.TickValues = [ylimlo (ylimlo+ylimhi)/2 ylimhi];
    yax.TickDirection = 'out';
    yax.TickLength = [yticklength yticklength];
    yax.TickLabels = {num2str(ylimlo, '%2.2f'), num2str((ylimlo+ylimhi)/2, '%2.2f'), num2str(ylimhi, '%2.2f')};
    yax.FontName = fontname;
    yax.FontSize = fontsize;
    yax.FontAngle = fontangle;
    
    % general
    a = gca;
    %     a.TitleFontWeight = 'normal';
    box off
    %
    legend('off')
    
    a.YLabel.String = 'RT Per Symbol (seconds)';
    a.YLabel.FontSize = fontsize;
    a.YLabel.FontAngle = fontangle;
    a.XLabel.String = 'Trial Number';
    a.XLabel.FontSize = fontsize;
    pbaspect([1 1 1])
    
    print(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-plots', ['plot_rt_scatter_fit' fittype '_sub' num2str(subjectlist(s))]), '-dpng')
    print(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-plots', 'eps', ['plot_rt_scatter_fit' fittype '_sub' num2str(subjectlist(s))]), '-depsc')
    
    hold off;
    
end

% Remove subs from data_recog_rt_mean who were not in this rt analysis.
keep = find(~ismember(data_recog_rt_mean.Var1, remove));
data_recog_rt_mean = data_recog_rt_mean(keep, :);
n = size(data_recog_rt_mean, 1);

% Concatenate learning rate to data_write_mean.
data_recog_rt_mean = sortrows(data_recog_rt_mean, 1);
data_recog_rt_mean.learningrate = learningrate';

% Create date-specific file name that indicates how many subjects.
filename = sprintf('WML_beh_data_recog_rt_wlearningrate_%s_n=%d_%s', fittype, size(data_recog_rt_mean, 1), datestr(now,'yyyymmdd'));

% Save all variables.
save(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-supportFiles', filename), 'data_recog', 'data_recog_rt_mean');

% Save as a CSV files.
writetable(data_recog, fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-supportFiles', [filename '.csv']))
writetable(data_recog_rt_mean, fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-supportFiles', [filename '_mean.csv']))
