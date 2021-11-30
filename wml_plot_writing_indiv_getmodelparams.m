clear all; close all; clc

% Set working directories.
rootDir = '/Volumes/Seagate/wml/';

% Create date-specific file name that indicates how many subjects.
datestring = '20211119';
remove = [24 31 51 66]; %[22 27 35 41 60 62];
filename = sprintf('wml_beh_data_write_%s', datestring);

fittype = 'exp2'; % poly1, exp2

% Load data.
load(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-supportFiles', filename), 'data_write', 'data_write_mean');

% Remove subject outliers.
keep = find(~ismember(data_write.subID, remove));
data_write = data_write(keep, :);

% Sort by subID, then by day, then by block, then by trial.
mt = sortrows(data_write, [1 3 5 6]);

% Remove day 5 for now.
mt = mt(find(mt.day ~= 5), :);

% % Use day 1 only for now.
% mt = mt(find(mt.day == 1), :);

% Add overall trial counter.
mt.idx = transpose(1:size(mt, 1));

% Set up binning vector.
count = 1;
for b = 1:200
    binning(count:count+9) = b;
    count = count + 10; 
end

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
        
        % Add binning for this subject.
        mt.trial_bin(sidx(1):sidx(end)) = binning(1:length(sidx));
        
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
gscatter(mt.trial_bin, mt.drawduration, mt.subID)

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

ylimlo = 1.5; ylimhi = 4.0;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [20 40 60 80 100 120 140 160];
% xax.TickDirection = 'out';
% xax.TickLength = [xticklength xticklength];
xlabels = {'20', '40', '60', '80', '100', '120', '140', '160'};
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
a.XLabel.String = 'Bin Number';
a.XLabel.FontSize = fontsize;
pbaspect([1 1 1])
n = length(subjectlist);
print(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-plots', ['plot_write_binned_scatter_n=' num2str(n)]), '-dpng')
print(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-plots', 'eps', ['plot_write_binned_scatter_n=' num2str(n)]), '-depsc')

hold off;

%% Plot binned means -- with all outliers removed.

for b = 1:160
    
    binmean(b) = mean(mt.drawduration(find(mt.trial_bin == b)));
    
end

% Plot trial x drawduration.
s = scatter(1:160, binmean);
s.Marker = '.';
s.SizeData = 200;

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

ylimlo = 1.5; ylimhi = 4.0;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [20 40 60 80 100 120 140 160];
% xax.TickDirection = 'out';
% xax.TickLength = [xticklength xticklength];
xlabels = {'20', '40', '60', '80', '100', '120', '140', '160'};
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

a.YLabel.String = {'Draw Duration Per Symbol (seconds)'; 'Mean across subjects'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = 'Bin Number';
a.XLabel.FontSize = fontsize;
pbaspect([1 1 1])
n = length(subjectlist);
print(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-plots', ['plot_write_meanbinned_scatter_n=' num2str(n)]), '-dpng')
print(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-plots', 'eps', ['plot_write_meanbinned_scatter_n=' num2str(n)]), '-depsc')

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


% Find idx for subjects of interest.
% sidx = find(mt.subID == subjectlist);

% Extract data for this fit.
y = transpose(binmean);
x = transpose(1:160);

% Fit.
[f, f2] = fit(x, y, fittype); %, 'Robust', 'Lar');

% Collect the learning rate.
if strcmp(fittype, 'exp2')
    
    kappa = f.d;
    lambda = f.b;
    learningrate = kappa;
    
elseif strcmp(fittype, 'poly1')
    
    learningrate(s) = f.p1;
    
end

figure(3);

% Plot fit.
p=plot(f, x, y);
hold on;

% Customize scatter plot of data.
p(1).Color = c(1, :);

% Customize line plot of data.
p(2).Color = c(1, :);
p(2).LineWidth = 1;

title([fittype ', adjr2 = ' num2str(f2.adjrsquare) ', rmse = ' num2str(f2.rmse)]);

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [20 40 60 80 100 120 140 160];
% xax.TickDirection = 'out';
% xax.TickLength = [xticklength xticklength];
xlabels = {'20', '40', '60', '80', '100', '120', '140', '160'};
% % xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;

% yaxis
ylimlo = 2.5; ylimhi = 3.5;
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

a.YLabel.String = {'Draw Duration Per Symbol (seconds)'; 'Mean across subjects'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = 'Trial Number';
a.XLabel.FontSize = fontsize;
pbaspect([1 1 1])

print(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-plots', ['plot_write_meanbinnedwfit_scatter_n=' num2str(n)]), '-dpng')
print(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-plots', 'eps', ['plot_write_meanbinnedwfit_scatter_n=' num2str(n)]), '-depsc')

hold off;


