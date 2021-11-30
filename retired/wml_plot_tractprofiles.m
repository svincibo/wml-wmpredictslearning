% This script reads in FA, OD, and ICVF measures (from Brad Caron's
% TractProfiles App) for each of the tracts generated (from Dan Bullock's
% White Matter Segmentation App). It also reads in tract statistics (e.g.,
% number of streamlines for each tract (from Dan Bullock's Check Tract
% Quality App). It also reads in behavioral data collected as part of the
% LWX study.

clear all; close all; clc
format shortG

rootdir = '/Volumes/Seagate/wml/wml-wmpredictslearning';
blprojectid = 'proj-61609744fc8eb92cbcf82a5d';

measures = {'fa'}; %, 'md'};
streamline_min = 100;
alphastat = 0.01;
binsize = 5;

c_color = [204 0 204]/255; %pink [0.6350 0.0780 0.1840]; %red
a_color = [0 0 0]; %[75 75 75]/255; % gray[0 0 0]; %black
hold on;
linewidth = 1.5;
linestyle = 'none';
fontname = 'Arial';
fontsize = 16;
fontangle = 'italic';
xticklength = 0;
ylimlo = 0.35; ylimhi = 0.60;

save_figures = 'yes';

% Load wm data from lwx_qa_tractstats.m.
load(fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', 'wml_wmlpredictslearning_data_streamlinecount.mat'));

% % Load streamline count data from lwx_datacat.m.
% load(fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', 'wml_data_fa_singleshell.mat'));

% Should outliers be removed? If so, which subIDs?
remove_outliers = 'no';
if strcmp(remove_outliers, 'yes')

    % Identify outliers to be removed.
    outlier = [108 116 119 125 126 206 214 303 317 318];
    
else
    
    outlier = [];
    
end

% Should we include children only or all subjects?
include = 'all'; % options: childrenonly, all

fcount = 0;
for w = 1:length(measures)
    
    measure = measures{w};
    
    %% TRACTOGRAPHY.
    
    % Get contents of the directory where the tract measures for this subject are stored.
    grp_contents = dir(fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', 'tractprofiles', ['*' measure '*.csv']));
    
    % Load in each tract's tractography measures for this subject.
    sub_count = 0;
    for s = 1:size(grp_contents, 1)
        
        sub(s) = str2num(grp_contents(s).name(end-5:end-4));
        
        % Display current sub ID.
        disp(grp_contents(s).name)
        
%         % Update subject counter for when not all subjects are used/needed.
%         sub_count = sub_count + 1;
        
        for t = 1:size(grp_contents, 1)
            
            % Read in data for this subject and this tract.
            data_temp = readtable(fullfile(grp_contents(t).folder, grp_contents(t).name));
%             data_array = table2array(data_temp); 
            
            if t == 1
                
                wm = data_temp;
                
            else

                wm = [wm; data_temp];
                
            end
            
            if t == size(grp_contents, 1)
              
                wm_out = [s wm];
                
            end
                
            clear data_temp 
            
        end % sub_contents
        
    end % group_contents
    
%     % Find empty cells and fill with 'empty'.
%     t = find(cellfun(@isempty,tract));
%     tract(t) = {'empty'};
    

    
    % Get a list of unique sub IDs.
    subID = unique(sub);
    
    % Plot tract profiles for each tract.
    for k = 1:size(list_tract, 1)/2
        
        % Only plot for tracts of interest and non-empty tracts.
        if strcmp(list_tract{k},'leftSLF1And2') || strcmp(list_tract{k}, 'rightSLF1And2') ...
                || strcmp(list_tract{k},'leftSLF3') || strcmp(list_tract{k}, 'rightSLF3') ...
                || strcmp(list_tract{k}, 'leftILF') || strcmp(list_tract{k}, 'rightILF') ...
                || strcmp(list_tract{k}, 'leftIFOF') || strcmp(list_tract{k}, 'rightIFOF') ...
                || strcmp(list_tract{k}, 'leftAslant') || strcmp(list_tract{k}, 'rightAslant') ...
                || strcmp(list_tract{k}, 'leftTPC') || strcmp(list_tract{k}, 'rightTPC') ...
                || strcmp(list_tract{k}, 'leftpArc') || strcmp(list_tract{k}, 'rightpArc') ...
                || strcmp(list_tract{k}, 'leftMDLFspl') || strcmp(list_tract{k}, 'rightMDLFspl') ...
                || strcmp(list_tract{k}, 'leftVOF') || strcmp(list_tract{k}, 'rightVOF') ...
                || strcmp(list_tract{k}, 'leftMDLFang') || strcmp(list_tract{k}, 'rightMDLFang') ...
                && ~strcmp(list_tract{k}, 'empty')
            
            disp(list_tract{k})
            
            % Find entries that are for this tract, include both left and right.
            t_idx = endsWith(tract, list_tract{k}(5:end));
            
            % Open a new figure for this tract.
            fcount = fcount + 1;
            figure(fcount)
            hold on;
            
            c_count = 0; a_count = 0;
            for s = 1:length(subID)
                
%                 % Identify subjects who have less than 100 streamlines for this track.
%                 col_idx = find(strcmp(streamlinecounts.Properties.VariableNames, list_tract{k}));
%                 subID_stream = streamlinecounts.subID(find(table2array(streamlinecounts(:, col_idx)) < streamline_min));
                
                % Only include subjects who are not outliers and who have at least 100 streamlines for this tract.
                if subID(s)~=0 && ~ismember(subID(s), outlier) %&& ~ismember(subID(s), subID_stream)
                    
                    % Find entries that are for this subject.
                    s_idx = sub == subID(s);
                    
                    % Subset the thing so that we only plot for this tract and for this subject.
                    t_temp = wm(find(t_idx == 1 & s_idx == 1));
                    
                    if ~isempty(t_temp)
                        
%                         % Code the plot for subject and keep data for inspection (yc, oc, a).
%                         if subID(s) < 300
%                             
%                             c_count = c_count + 1;
%                             
%                             % Young child.
%                             plot(t_temp, 'LineWidth', 2, 'LineStyle', '-', 'Color', [c_color .2])
%                             
%                             % Collect.
%                             c(:, c_count) = t_temp;                           
%                             
%                         else
                            
                            a_count = a_count + 1;
                            
                            % Adult.
                            plot(t_temp, 'LineWidth', 2, 'LineStyle', '-', 'Color', [a_color .2])
                            
                            % Collect.
                            a(:, a_count) = t_temp;
                            
%                         end
                        
%                         hold on;
                        
                        clear t_temp;
                        
                    end % if subID{s} ~= 0
                    
%                     % Only include subjects who are not outliers and who
%                     % have at least 100 streamlines for this tract. Set
%                     % cells with less than 100 streamlines to NaN.
%                 elseif subID(s)~=0 && ~ismember(subID(s), outlier) && ismember(subID(s), subID_stream)
                    
                end %if exist
                
            end %sub
            
            % Plot means and 95% confidence intervals (calculated from
            % standard error: 1.96*SE). Use nanmean/nanstd because one oc subject is missing TPC.
%             plot(nanmean(c, 2), 'LineWidth', 4, 'LineStyle', '-', 'Color', c_color(1:3))
%             hi = nanmean(c, 2) + 1.96*nanstd(c, 0, 2)/sqrt(size(~isnan(c), 2)); lo = nanmean(c, 2) - 1.96*nanstd(c, 0, 2)/sqrt(size(~isnan(c), 2)); x = (1:size(nanmean(c, 2),1))';
%             hp1 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], c_color(1:3));
%             set(hp1, 'facecolor', c_color(1:3), 'edgecolor', 'none', 'facealpha', .2);
%             
            plot(nanmean(a, 2), 'LineWidth', 4, 'LineStyle', '-', 'Color', a_color(1:3))
            hi = nanmean(a, 2) + 1.96*nanstd(a, 0, 2)/sqrt(size(~isnan(a), 2)); lo = nanmean(a, 2) - 1.96*nanstd(a, 0, 2)/sqrt(size(~isnan(a), 2)); x = (1:size(nanmean(a, 2),1))';
            hp3 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], a_color(1:3));
            set(hp3, 'facecolor', a_color(1:3), 'edgecolor', 'none', 'facealpha', .2);
            
            % xaxis
            xax = get(gca, 'xaxis');
            xax.Limits = [0 160];
            xax.TickValues = [0 80 160];
            xax.TickLabels = {'20', '100', '180'};
            xax.TickDirection = 'out';
            xax.TickLength = [xticklength xticklength];
            xax.FontName = fontname;
            xax.FontSize = fontsize;
            xax.FontAngle = fontangle;
            
            % yaxis
            yax = get(gca,'yaxis');
            yax.Limits = [ylimlo ylimhi];
            yax.TickValues = [ylimlo (ylimlo+ylimhi)/2 ylimhi];
            yax.TickDirection = 'out';
            yax.TickLabels = {num2str(ylimlo, '%1.2f'), num2str((ylimlo+ylimhi)/2, '%1.2f'), num2str(ylimhi, '%1.2f')};
            yax.FontName = fontname;
            yax.FontSize = fontsize;
            
            % general
            g = gca;
            %     a.TitleFontWeight = 'normal';
            box off
            
            % legend({'Younger Children', 'Older Children', 'Adults'}, 'Location', 'southeast');
            % legend box off
            
            title(list_tract{k}(5:end))
            g.XLabel.String = 'Location along tract';
            g.XLabel.FontSize = fontsize;
            g.XLabel.FontAngle = fontangle;
            
            g.YLabel.String = 'Fractional Anisotropy (FA)';
            g.YLabel.FontSize = fontsize;
            pbaspect([1 1 1])
            
            print(fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', ['wml_plot_tractprofiles_' measure '_' list_tract{k}(5:end)]), '-dpng')
            print(fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', 'eps', ['wml_plot_tractprofiles_' measure '_' list_tract{k}(5:end)]), '-depsc')
            
            hold off;
            
%             % Open a new figure for the mean plot.
%             fcount = fcount + 1;
%             figure(fcount)
%             
%             gr = [1 2];
%             gscatter([0.5 1.5], [nanmean(c, 'all') nanmean(a, 'all')], gr, cat(1, c_color, a_color))
%             hold on;
%             plot([0.5 0.5], [nanmean(c, 2) + 1.96*nanstd(c, 0, 2)/sqrt(size(~isnan(c), 2)) nanmean(c, 2) - 1.96*nanstd(c, 0, 2)/sqrt(size(~isnan(c), 2))], 'LineStyle', '-', 'Color', c_color(1:3))
%             plot([2.5 2.5], [nanmean(a, 2) + 1.96*nanstd(a, 0, 2)/sqrt(size(~isnan(a), 2)) nanmean(a, 2) - 1.96*nanstd(a, 0, 2)/sqrt(size(~isnan(a), 2))], 'LineStyle', '-', 'Color', a_color(1:3))
%             
%             legend off
%             
%             % xaxis
%             xax = get(gca, 'xaxis');
%             xax.Limits = [0 3];
%             xax.TickValues = [0.5 1.5 2.5];
%             xax.TickLabels = {'Children', 'Adults'};
%             xax.TickDirection = 'out';
%             xax.TickLength = [xticklength xticklength];
%             xax.FontName = fontname;
%             xax.FontSize = fontsize;
%             xax.FontAngle = fontangle;
%             
%             % yaxis
%             yax = get(gca,'yaxis');
%             yax.Limits = [ylimlo ylimhi];
%             yax.TickValues = [ylimlo (ylimlo+ylimhi)/2 ylimhi];
%             yax.TickDirection = 'out';
%             yax.TickLabels = {num2str(ylimlo, '%1.2f'), num2str((ylimlo+ylimhi)/2, '%1.2f'), num2str(ylimhi, '%1.2f')};
%             yax.FontName = fontname;
%             yax.FontSize = fontsize;
%             
%             % general
%             g = gca;
%             %     a.TitleFontWeight = 'normal';
%             box off
%             
%             % legend({'Younger Children', 'Older Children', 'Adults'}, 'Location', 'southeast');
%             % legend box off
%             
%             title(list_tract{k})
%             
%             g.YLabel.String = 'Fractional Anisotropy (FA)';
%             g.YLabel.FontSize = fontsize;
%             pbaspect([1 1 1])
%             
%             print(fullfile(rootdir, 'plots-singleshell', ['plot_meancomparison_singleshell_' wm_measure '_' list_tract{k}]), '-dpng')
%             print(fullfile(rootdir, 'plots-singleshell', 'eps', ['plot_meancomparison_singleshell_' wm_measure '_' list_tract{k}]), '-depsc')
%             
%             hold off;
            
            clear a c 
            
        end % if toi
        
    end %tract
    
end % for w

% Get only non-empty cells.
figure
hold on;
idx = find(~isnan(z_keep));
z_new = z_keep(idx);
z_new_tractname = z_keep_tractname(idx);
tn = unique(z_new_tractname);

for r = 1:length(tn)
   
    idx_temp = find(strcmp(z_new_tractname, tn(r)));
    
    scatter(z_new(idx_temp), r*ones(size(z_new(idx_temp))), '*k', 'SizeData', 100);
    
    clear idx_temp
    
end

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [1.5 3];
xax.TickValues = [1.5 2.25 3];
xax.TickLabels = {'1.5', '2.25', '3.00'};
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;

% yaxis
yax = get(gca,'yaxis');
ylimlo = 0; ylimhi = 20;
yax.Limits = [ylimlo ylimhi];
yax.TickValues = 1:19;
yax.TickDirection = 'out';
yax.TickLabels = tn;
yax.FontName = fontname;
yax.FontSize = fontsize;

% general
g = gca;
%     a.TitleFontWeight = 'normal';
box off

% legend({'Younger Children', 'Older Children', 'Adults'}, 'Location', 'southeast');
% legend box off

g.XLabel.String = 'Effect Size (dprime)';
g.XLabel.FontSize = fontsize;

pbaspect([1 1 1])

print(fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', ['wml_plot_zhist_' measure]), '-dpng')
print(fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', 'eps', ['wml_plot_zhist_' measure]), '-depsc')

hold off;
