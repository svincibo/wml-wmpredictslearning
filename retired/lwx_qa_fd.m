clear all; close all; clc
format shortG

a_color = [0 0 0]; %[75 75 75]/255; % gray [.146 0 0]; % light black
c_color = [204 0 204]/255; %pink [178 34 34]/255; % firebrick red [0 .73 .73]; % turquoise

% Set working directories.
rootDir = '/Volumes/240/lwx/';

% Get bl project foldername.
blprojectid = 'proj-5e849e65952fef3dcd7a1700';

% Read in behavioral data.
beh = readtable([rootDir 'supportFiles/LWX_all_groupings.csv'], 'TreatAsEmpty', {'.', 'na'});

% Should outliers be removed? If so, which subIDs?
remove_outliers = 'yes';
if strcmp(remove_outliers, 'yes')
    
    % Identify outliers to be removed - conservative removal.
%         outlier = [108 126 214 318];
    % 108, snr is below 2 SD of group mean
    % 126, dwi image has major distortions, visual inspection
    % 214, major motion artifacts, visual inspection
    % 318, snr is below 2 SD of group mean and dwi image has major distortions, visual inspection
    
    % Identify outliers to be removed - liberal removal.
    outlier = [108 116 125 126 203 206 212 214 315 316 318];
    % 116, FD > 2
    % 119, FD > 2
    % 125, FD > 2
    % 206, FD > 2 
    % 303, SNR < 15 
    % 317, FD > 2
    
end

% Should we include children only or all subjects?
include = 'all'; % options: childrenonly, all

% Get contents of the directory where the tract measures for this subject are stored.
grp_contents = dir(fullfile(rootDir, blprojectid));

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

% Keep only names that are subject folders.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');

% For each subject.
sub_count = 0;
for s = 1:size(grp_contents, 1)
    
    % Only collect values for subjects that have both MRI and behaviora/demographic data.
    if ~isempty(find((beh.SubjectID == str2num(grp_contents(s).name(5:7)))))
        
        % Display current sub ID.
        disp(grp_contents(s).name)
        
        % Get contents of the directory where the SNR values for this subject are stored.
        sub_contents_motion = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, '/dt-neuro-dtiinit*/*ecXform.mat'));
        % Remove the '.' and '..' files.
        sub_contents_motion = sub_contents_motion(arrayfun(@(x) x.name(1), sub_contents_motion) ~= '.');
        
        if ~isempty(sub_contents_motion)
            
            % Update subject counter for when not all subjects are used/needed.
            sub_count = sub_count + 1;
            
            % Get contents of the directory where the SNR values for this subject are stored.
            sub_contents_motion = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, '/dt-neuro-dtiinit*/*ecXform.mat'));
            % Remove the '.' and '..' files.
            sub_contents_motion = sub_contents_motion(arrayfun(@(x) x.name(1), sub_contents_motion) ~= '.');
            
            % Get motion parameters for this subject.
            load(fullfile(sub_contents_motion.folder,sub_contents_motion.name));
            
            % Select only the translation/rotation parameters.
            mot = vertcat(xform(:).ecParams);
            mot = mot(:, 1:6); % xyz translation xyz rotation
            
            % Motion parameters represent the translation/rotation that must occur
            % to return the image at timepoint tn to the place that it was at timepoint
            % t0. Thus, to calculate instantaneouse parameters, we need a moving
            % difference.
            for m = 1:size(mot, 2)
                
                % Get moving difference for each parameter. Append row of zeros for t1, per convention (Power et al., 2014).
                % This step creates the fd timecourses for each motion parameter for each run in the motion struct and represents instantaneous movement.
                movingdifference(:, m) = [0 ; diff(mot(:, m), 1, 1)]';
                
            end
            
            % Get an overall fd for all 6 parameters for each run.
            % This step creates the fd summary statistic for all 6 parameters for each timepoint in a run for each run (e.g., scalar FD timecourse).
            motion(sub_count, :) = sum(abs(movingdifference), 2)';
            
            % Get subID.
            subID(sub_count) = str2num(grp_contents(s).name(5:7));
            
            % Get ylabel.
            lab{sub_count} = grp_contents(s).name;
            
            % Get session.
            age(sub_count) = beh.Age_months(find((beh.SubjectID == str2num(grp_contents(s).name(5:7)))));
            
            % Get training group.
            group(sub_count) = beh.group_age(find((beh.SubjectID == str2num(grp_contents(s).name(5:7)))));
            
        end % if ~isempty
        
    end % end if exist
    
end % end s


% Remove outliers.
if strcmp(remove_outliers, 'yes') && exist('outlier')
    
    % Get index for outliers to be removed.
    idx_outlier = ismember(subID, outlier);
    
    % Remove outliers.
    subID = subID(~idx_outlier);
    motion = motion(~idx_outlier, :);
    group = group(~idx_outlier);
    age = age(~idx_outlier);
    
    % Set figure note.
    ttlstr = 'Motion outlier removed.';
    
else
    
    % Set figure note.
    ttlstr = 'Motion outlier retained.';
    
end

meanmotion = mean(motion, 2);

% Write out table for anova.
t_out = array2table(cat(2, subID', age', group', meanmotion), 'VariableNames', {'subID', 'session', 'group', 'fd'});
writetable(t_out, fullfile(rootDir, 'supportFiles', 'spade_data_motion.csv'));

disp('Check for group differences in FD.')
[~, tableout, ~] = anova1(meanmotion, group', 'off');
disp(['F(' num2str(tableout{2, 3}) ', ' num2str(tableout{3, 3}) ') = ' num2str(tableout{2, 5}) ', p = ' num2str(tableout{2, 6}) '.'])

disp('Check for child vs. adult differences in FD.')
[h, p, ci stats] = ttest2(meanmotion(group ~= 3), meanmotion(group == 3));
disp(['t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', p = ' num2str(p) '.'])

% Visualize: group differences
figure(2)
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

b1 = bar(1, nanmean(meanmotion(group ~= 3)), 'FaceColor', c_color, 'EdgeColor', c_color, 'FaceAlpha', alphablend);
plot([1 1], [nanmean(meanmotion(group ~= 3)) - nanstd(meanmotion(group ~= 3)) nanmean(meanmotion(group ~= 3)) + nanstd(meanmotion(group ~= 3))], 'Color', c_color)
b2 = bar(2, nanmean(meanmotion(group == 3)), 'FaceColor', a_color, 'EdgeColor', a_color, 'FaceAlpha', alphablend);
plot([2 2], [nanmean(meanmotion(group == 3)) - nanstd(meanmotion(group == 3)) nanmean(meanmotion(group == 3)) + nanstd(meanmotion(group == 3))], 'Color', a_color)

% xlim_lo = min(age)-1;
% xlim_hi = max(age)+1;
ylim_lo = 0;
ylim_hi = 3; %max(meanmotion)+0.5;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0.5 2.5];
xax.TickValues = [1 2];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xlabels = {'Children', 'Adults'};
xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;

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

a.YLabel.String = 'Mean Framewise Displacement (FD)';
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
pbaspect([1 1 1])

print(fullfile(rootDir, 'plots-singleshell', 'plot_barplot_fd'), '-dpng')
print(fullfile(rootDir, 'plots-singleshell', 'eps', 'plot_barplot_fd'), '-depsc')

hold off;

% % Manually record outliers. Include observations with unusually high FD and any above 2mm.
% % (0 indicates no outliers)
% outliers.motion = 90;
% % outliers.motion = subID(meanmotion>2);
%
% save('devti_remove_motionoutliers.mat', 'outliers')
