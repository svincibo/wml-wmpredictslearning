clear all; close all; clc
format shortG

a_color = [0 0 0]; %[75 75 75]/255; % gray [.146 0 0]; % light black
c_color = [204 0 204]/255; %pink [178 34 34]/255; % firebrick red [0 .73 .73]; % turquoise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set working directories.
rootDir = '/Volumes/Seagate/ping';

% Get bl project foldername.
blprojectid = 'dtiinit-proj-60708cf9c7f80a684995e0b1';

% Read in behavioral data.
beh = readtable(fullfile(rootDir, 'supportFiles', 'ping_combined_participantinformation.csv')); %, 'TreatAsEmpty', {'.', 'na'});

% Get subIDs present in beh.
beh_subjects = beh.subID;

% Append age grouping variable that is double.
idx_child = contains(beh.age_mri_tag, 'child_upto8.5');
idx_preteen = contains(beh.age_mri_tag, 'child_8.6-13.5');
idx_adolescent = contains(beh.age_mri_tag, 'child_13.6-18.5');
idx_adult = contains(beh.age_mri_tag, 'adult');
beh_out.gp_age = idx_child + 2*idx_preteen + 3*idx_adolescent + 4*idx_adult;
    
% Get groups from beh.
beh_group = beh_out.gp_age;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get contents of the directory where the tract measures for this subject are stored.
grp_contents = dir(fullfile(rootDir, blprojectid));

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

% Keep only names that are subject folders.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');

% For each subject.
sub_count = 0;
            % Preallocate -- quick fix.
%             movingdifference = NaN(62, size(grp_contents, 1));
for s = 1:size(grp_contents, 1)
    
    % Only collect values for subjects that have both MRI and behaviora/demographic data.
    if sum(contains(beh_subjects, grp_contents(s).name(end-4:end)))
        
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
                movingdifference(1:length([0 ; diff(mot(:, m), 1, 1)]'), m) = [0 ; diff(mot(:, m), 1, 1)]';
                
            end
            
            % Get an overall fd for all 6 parameters for each run.
            % This step creates the fd summary statistic for all 6 parameters for each timepoint in a run for each run (e.g., scalar FD timecourse).
            motion(sub_count, :) = sum(abs(movingdifference), 2)';
            
            % Get subID.
            subID{sub_count} = grp_contents(s).name(end-4:end);
            
            % Get ylabel.
            lab{sub_count} = grp_contents(s).name(5:end);
            
            % Get session.
            age(sub_count) = beh.age_mri(contains(beh_subjects, grp_contents(s).name(end-4:end)));
            
            % Get training group.
            group(sub_count) = beh_group(contains(beh_subjects, grp_contents(s).name(end-4:end)));
            
        end % if ~isempty
        
    end % end if exist
    
end % end s

% keep only group = 1 (children) and group = 4 (adults)

% Concatenate and sort according to group.
idx = find(group == 1 | group == 4);

meanmotion = nanmean(motion, 2);
sd = std(motion, 0, 2);

% Concatenate and sort according to group.
toplot = cat(2, group(idx)', meanmotion(idx), sd(idx));
toplot = sortrows(toplot, [2 3 1], 'ascend');

% cheap trick
clear group m sd b0
group = toplot(:, 1)';
fd = toplot(:, 2)';
sd = toplot(:, 3)';

% FD plot
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

gscatter(fd, 1:length(fd), group, [c_color; a_color], '.', 30)

hold on;
for p = 1:length(fd)
    
    if group(p) == 1 
        
        plot([fd(p) - abs(sd(p)) fd(p) + abs(sd(p))], [p p], 'Color', c_color)
        
    elseif group(p) == 4 
        
        plot([fd(p) - abs(sd(p)) fd(p) + abs(sd(p))], [p p], 'Color', a_color)
        
    end
    
end

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0 6];
xax.TickValues = 0:1:6;
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [0 length(lab)+0.5];
yax.TickValues = 1:1:ceil(length(lab));
yax.TickDirection = 'out';
yax.TickLabels = subID;
yax.FontName = fontname;
yax.FontSize = 8;
plot([1 1], [0 length(subID)+0.5], ':k')

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

legend({'Children', 'Adults'}, 'Location', 'southeast');
legend box off

a.XLabel.String = 'Framewise Displacement (FD)';
a.XLabel.FontSize = fontsize;
a.XLabel.FontAngle = fontangle;

a.YLabel.String = 'Participant ID';
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
pbaspect([1 1 1])

print(fullfile(rootDir, 'plots', 'plot_fd'), '-dpng')

hold off;