% This script reads in SNR values and plots them according to session (pre-
% vs post-training) and group (expert=3, beginner=2, control=1).

clear all; close all; clc
format shortG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rootdir = '/Volumes/Seagate/wml/wml-wmpredictslearning';
blprojectid = 'proj-61609744fc8eb92cbcf82a5d';

% Read in behavioral data.
% beh = readtable(fullfile(rootdir, 'supportFiles', 'ping_combined_participantinformation.csv')); %, 'TreatAsEmpty', {'.', 'na'});

% % Get subIDs present in beh.
% beh_subjects = beh.subID;

% % Append age grouping variable that is double.
% idx_child = contains(beh.age_mri_tag, 'child_upto8.5');
% idx_preteen = contains(beh.age_mri_tag, 'child_8.6-13.5');
% idx_adolescent = contains(beh.age_mri_tag, 'child_13.6-18.5');
% idx_adult = contains(beh.age_mri_tag, 'adult');
% beh_out.gp_age = idx_child + 2*idx_preteen + 3*idx_adolescent + 4*idx_adult;
%     
% % Get groups from beh.
% beh_group = beh_out.gp_age;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get contents of the directory where the tract measures for this subject are stored.
grp_contents = dir(fullfile(rootdir, blprojectid));

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

% Keep only names that are subject folders.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');

% For each subject.
sub_count = 0;
for s = 1:size(grp_contents, 1)
    
    % Only collect values for subjects that have both MRI and behavioral/demographic data.
%     if sum(contains(beh_subjects, grp_contents(s).name(end-4:end))) > 0
        
        % Display current sub ID.
        disp(grp_contents(s).name)
        
        % Get contents of the directory where the SNR values for this subject are stored.
        sub_contents_snr = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, '/dt-raw.tag-snr*/output/snr.json'));
        % Remove the '.' and '..' files.
        sub_contents_snr = sub_contents_snr(arrayfun(@(x) x.name(1), sub_contents_snr) ~= '.');
        
        if ~isempty(sub_contents_snr)
            
            % Update subject counter for when not all subjects are used/needed.
            sub_count = sub_count + 1;
            
            % Get SNR for this subject.
            data_snr_temp = jsondecode(fileread(fullfile(sub_contents_snr.folder, sub_contents_snr.name)));
            
            % Get SNR in b0 images.
            b0(sub_count) = str2num(data_snr_temp.SNRInB0_X_Y_Z{1});
            
            % Get mean SNR in X, Y, and Z directions.
            m(sub_count) = mean([str2num(data_snr_temp.SNRInB0_X_Y_Z{2}), str2num(data_snr_temp.SNRInB0_X_Y_Z{3}), str2num(data_snr_temp.SNRInB0_X_Y_Z{4})]);
            
            % Get standard deviation of SNR in X, Y, and Z directions.
            sd(sub_count) = std([str2num(data_snr_temp.SNRInB0_X_Y_Z{2}), str2num(data_snr_temp.SNRInB0_X_Y_Z{3}), str2num(data_snr_temp.SNRInB0_X_Y_Z{4})]);
            
            % Get subID.
            subID{sub_count} = grp_contents(s).name(end-8:end-6);
            
            % Get training group.
%             group(sub_count) = beh_group(contains(beh_subjects, grp_contents(s).name(end-4:end)));
            
            clear data_snr_temp get_temp
            
%         end % if ~isempty
        
    end % end if exist
    
end % end s

% keep only group = 1 (children) and group = 4 (adults)

% Concatenate and sort according to group.
% idx = find(group == 1 | group == 4);
toplot = cat(2, m', sd', b0');
toplot = sortrows(toplot, [2 3 1], 'ascend');

% cheap trick
clear group m sd b0
m = toplot(:, 1)';
sd = toplot(:, 2)';
b0 = toplot(:, 3)';

% SNR plot
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
xticklength = 0;
alphablend = .8;

c_color = [0 0 0]; %[75 75 75]/255; % gray [.146 0 0]; % light black
a_color = [204 0 204]/255; %pink [178 34 34]/255; % firebrick red [0 .73 .73]; % turquoise

scatter(m, 1:length(m), 'MarkerEdgeColor', c_color, 'MarkerFaceColor', c_color, 'Marker', '.', 'SizeData', 1200)
hold on;
scatter(b0, 1:length(b0), 'MarkerEdgeColor', c_color, 'MarkerFaceColor', c_color, 'Marker', 'x', 'SizeData', 100)
plot([20 20], [0 length(m)], ':k')

for p = 1:length(m)
        
        plot([m(p) - abs(sd(p)) m(p) + abs(sd(p))], [p p], 'Color', c_color)
    
end

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0 150];
xax.TickValues = 0:25:150;
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xax.FontName = fontname;
xax.FontSize = fontsize;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [0 34];
yax.TickValues = 1:34;
yax.TickDirection = 'out';
yax.TickLabels = subID;
yax.FontName = fontname;
yax.FontSize = 8;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off
legend off

a.XLabel.String = 'Signal-to-Noise Ratio (SNR)';
a.XLabel.FontSize = fontsize;
a.XLabel.FontAngle = fontangle;

a.YLabel.String = 'Participant ID';
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
pbaspect([1 1 1])

print(fullfile(rootdir, 'wml-wmpredictslearning-plots', 'plot_snr'), '-dpng')

hold off;