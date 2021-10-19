clear all; close all; clc

% Set working directories.
rootDir = '/Volumes/Seagate/wml/';

% Create date-specific file name that indicates how many subjects.
datestring = '20211019';
remove = []; %[22 27 35 41 60 62];
filename = sprintf('wml_beh_data_write_%s', datestring);

fittype = 'poly1'; % poly1, exp2

% Load data.
load(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-supportFiles', filename), 'data_write', 'data_write_mean');

% Remove subject outliers.
keep = find(~ismember(data_write.subID, remove));
data_write = data_write(keep, :);

% Sort by subID, then by day, then by trial.
mt = sortrows(data_write, [1 3 5 6]);

% Remove day 5 for now.
mt = mt(find(mt.day ~= 5), :);

% Add overall trial counter.
mt.idx = transpose(1:size(mt, 1));

% Remove outliers based on subject/day means and timeouts -- right now assumed timeout if drawduration = 4,
% which is trialduration. Future will incorporate data from visual inspection of the drawn images.
idx_remove = [];
subjectlist = unique(data_write.subID);
subjectlist = subjectlist(~isnan(subjectlist));
for s = 1:length(subjectlist)
    
    % Find idx for this subject.
    sidx = find(mt.subID == subjectlist(s));
    
    if ~isempty(sidx)
        
        % Add overal trial counter for this subject.
        mt.trial_all(sidx(1):sidx(end)) = transpose(1:length(sidx));
        
        for d = 1:4
            
            didx = find(mt.subID == subjectlist(s) & mt.day == d);
            
            datanow = mt.drawduration(didx, :);
            m = nanmean(datanow); sd = nanstd(datanow);
            
            % Get lower and upper bound.
            idx_above = find(datanow > m+3*sd | datanow >= 4);
            idx_below = find(datanow < m-3*sd);
            % Remove outliers in DI.
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
mt(idx_remove, :) = [];

% Remove NaNs.
mt = mt(~isnan(mt.drawduration), :);

%% Plot scatter plot of individual data -- with all outliers removed.

% Plot trial x drawduration.
gscatter(mt.trial_all, mt.drawduration, mt.subID)

% Set up plot and measure-specific details.
capsize = 0;
marker = 'o';
linewidth = 1.5;
linestyle = 'none';
markersize = 100;
xtickvalues = 1:1600;
xlim_lo = min(xtickvalues); xlim_hi = max(xtickvalues);
fontname = 'Arial';
fontsize = 15;
fontangle = 'italic';
yticklength = 0;
xticklength = 0.05;
coloralpha = .05;

ylimlo = 1.5; ylimhi = 4.0;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [200 400 600 800 1000 1200 1400 1600];
% xax.TickDirection = 'out';
% xax.TickLength = [xticklength xticklength];
xlabels = {'200', '400', '600', '800', '1000', '1200', '1400', '1600'};
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

print(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-plots', 'plot_write_scatter'), '-dpng')
print(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-plots', 'eps', 'plot_write_scatter'), '-depsc')

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
    sidx = find(mt.subID == subjectlist(s));
    
    % Extract data for this fit.
    y = mt.drawduration(sidx, :);
    x = mt.trial_all(sidx, :);
    
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
    xax.TickValues = [200 400 600 800 1000 1200 1400 1600];
    % xax.TickDirection = 'out';
    % xax.TickLength = [xticklength xticklength];
    xlabels = {'200', '400', '600', '800', '1000', '1200', '1400', '1600'};
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
    
    print(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-plots', ['plot_write_scatter_fit' fittype '_sub' num2str(subjectlist(s))]), '-dpng')
    print(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-plots', 'eps', ['plot_write_scatter_fit' fittype '_sub' num2str(subjectlist(s))]), '-depsc')
    
    hold off;
    
end

% Concatenate learning rate to data_write_mean.
data_write_mean = sortrows(data_write_mean, 1);
data_write_mean.learningrate = learningrate';

% Create date-specific file name that indicates how many subjects.
filename = sprintf('WML_beh_data_write_wlearningrate_%s_n=%d_%s', fittype, size(data_write_mean, 1), datestr(now,'yyyymmdd'));

% Save all variables.
save(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-supportFiles', filename), 'data_write', 'data_write_mean');

% Save as a CSV files.
writetable(data_write, fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-supportFiles', [filename '.csv']))
writetable(data_write_mean, fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-supportFiles', [filename '_mean.csv']))



