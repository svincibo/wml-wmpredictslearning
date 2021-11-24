% This script reads in SNR values and plots them according to session (pre-
% vs post-training) and group (expert=3, beginner=2, control=1).

clear all; close all; clc
format shortG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rootdir = '/Volumes/Seagate/wml/wml-wmpredictslearning';
blprojectid = 'proj-61609744fc8eb92cbcf82a5d';

% % Read in behavioral data.
% beh = readtable([rootdir 'supportFiles/LWX_all_groupings.csv'], 'TreatAsEmpty', {'.', 'na'});
% 
% % Get subIDs present in beh.
% beh_subjects = beh.SubjectID;
% 
% % Get groups from beh.
% beh_group = beh.group_age;

% Should outliers be removed? If so, which subIDs?
remove_outliers = 'no';
if strcmp(remove_outliers, 'yes')

    % Identify outliers to be removed.
    outlier = [108 116 125 126 203 206 212 214 315 316 318];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get contents of the directory where the tract measures for this subject are stored.
grp_contents = dir(fullfile(rootdir, blprojectid));

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

% Keep only names that are subject folders.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');

% Load in each tract's tractography measures for this subject.
sub_count = 0;
for s = 1:size(grp_contents, 1)
    
    % Only collect values for subjects that have both MRI and behaviora/demographic data.
%     if ~isempty(find((beh_subjects  == str2num(grp_contents(s).name(5:7)))))
        
        % Only collect values for subjects that have diffusion data (noted by whether or not an snr file was created).
        t = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, '/dt-raw.tag-snr*/product.json'));
        
        if numel(t) ~= 0
            %             exist(fullfile(grp_contents(s).folder, grp_contents(s).name, 'dt-raw.tag-snr-cc.id-5e8d9d8afa1b637b4f2f9350/product.json'))
            % Display current sub ID.
            disp(grp_contents(s).name)
            
            % Update subject counter for when not all subjects are used/needed.
            sub_count = sub_count + 1;
            
            % Get contents of the directory where the SNR values for this subject are stored.
            sub_contents_snr = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, '/dt-raw.tag-snr*/*product.json'));
            % Remove the '.' and '..' files.
            sub_contents_snr = sub_contents_snr(arrayfun(@(x) x.name(1), sub_contents_snr) ~= '.');
            
            % Get SNR for this subject.
            data_snr_temp = jsondecode(fileread([sub_contents_snr.folder filesep sub_contents_snr.name]));
            
            % Get SNR in b0 images.
            b0(sub_count) = str2num(data_snr_temp.SNRInB0_X_Y_Z{1});
            
            % Get mean SNR in X, Y, and Z directions.
            m(sub_count) = mean([str2num(data_snr_temp.SNRInB0_X_Y_Z{2}), str2num(data_snr_temp.SNRInB0_X_Y_Z{3}), str2num(data_snr_temp.SNRInB0_X_Y_Z{4})]);
            
            % Get subID.
            subID(sub_count) = str2num(grp_contents(s).name(5:7));
            
%             % Get training group.
%             group(sub_count) = beh_group(find((beh_subjects == str2num(grp_contents(s).name(5:7)))));
            
            clear data_snr_temp get_temp t
            
        end % if exist
        
%     end % end if isempty
    
end % end t

% Remove outliers.
if strcmp(remove_outliers, 'yes') && exist('outlier')
    
    % Get index for outliers to be removed.
    idx_outlier = ismember(subID, outlier);
    
    % Remove outliers.
    subID = subID(~idx_outlier);
    b0 = b0(~idx_outlier);
    m = m(~idx_outlier);
    group = group(~idx_outlier);
    
    % Set figure note.
    ttlstr = 'SNR outlier removed.';
    
else
    
    % Set figure note.
    ttlstr = 'SNR outlier retained.';
    
end

% Write out table for anova.
t_out = array2table(cat(2, subID', group', m', b0'), 'VariableNames', {'subID', 'group', 'm', 'b0'});
writetable(t_out, fullfile(rootdir, 'supportFiles', 'wml_data_snr.csv'));

% Group differences plot: b0
snr = b0;
figure(1)
hold on;
capsize = 0;
marker = 'o';
linewidth = 1.5;
linestyle = 'none';
markersize = 10;
fontname = 'Arial';
fontsize = 16;
fontangle = 'italic';
yticklength = 0;
xticklength = 0.05;
xtickvalues = [1 2 3];
alphablend = .8;
ylim_lo = 0;
ylim_hi = 50;

a_color = [0 0 0]; %[75 75 75]/255; % gray [.146 0 0]; % light black

b2 = bar(2, nanmean(snr), 'FaceColor', a_color, 'EdgeColor', a_color, 'FaceAlpha', alphablend);
plot([2 2], [nanmean(snr) - nanstd(snr) nanmean(snr) + nanstd(snr)], 'Color', a_color)

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0.5 2.5];
xax.TickValues = [1 2];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
% xlabels = {'Children', 'Adults'};
xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.0f'), '', num2str(ylim_hi, '%1.0f')};
yax.FontName = fontname;
yax.FontSize = fontsize;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

a.YLabel.String = 'SNR, b0 volumes';

a.YLabel.FontSize = fontsize;
pbaspect([1 1 1])

print(fullfile(rootdir, 'wml-wmpredictslearning-plots', ['wml_plot_barplot_snr_b0_outliers' remove_outliers]), '-dpng')
print(fullfile(rootdir, 'wml-wmpredictslearning-plots', 'eps', ['wml_plot_barplot_snr_b0_outliers' remove_outliers]), '-depsc')

hold off;
clear snr

% Group differences plot: weighted
snr = m;
figure(2)
hold on;

b2 = bar(2, nanmean(snr), 'FaceColor', a_color, 'EdgeColor', a_color, 'FaceAlpha', alphablend);
plot([2 2], [nanmean(snr) - nanstd(snr) nanmean(snr) + nanstd(snr)], 'Color', a_color)

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0.5 2.5];
xax.TickValues = [1 2];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
% xlabels = {'Children', 'Adults'};
xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.0f'), '', num2str(ylim_hi, '%1.0f')};
yax.FontName = fontname;
yax.FontSize = fontsize;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

a.YLabel.String = 'SNR, weighted volumes';

a.YLabel.FontSize = fontsize;
pbaspect([1 1 1])

print(fullfile(rootdir, 'wml-wmpredictslearning-plots', ['wml_plot_barplot_snr_weighted_outliers' remove_outliers]), '-dpng')
print(fullfile(rootdir, 'wml-wmpredictslearning-plots', 'eps', ['wml_plot_barplot_snr_outliers' remove_outliers]), '-depsc')

hold off;
