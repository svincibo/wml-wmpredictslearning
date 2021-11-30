clear all; close all; clc
format shortG

a_color = [0 0 0]; %[75 75 75]/255; % gray [.146 0 0]; % light black
c_color = [204 0 204]/255; %pink [178 34 34]/255; % firebrick red [0 .73 .73]; % turquoise

blprojectid = 'proj-5e61139282b37f2cfe8fdb28';

% Set working directories.
rootDir = '/Volumes/240/lwx/';

% Get bl project foldername.
blprojectid = 'proj-5e849e65952fef3dcd7a1700';

% Read in behavioral data.
beh = readtable([rootDir 'supportFiles/LWX_all_groupings.csv'], 'TreatAsEmpty', {'.', 'na'});

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
            lab{sub_count} = grp_contents(s).name(5:end);
            
            % Get session.
            age(sub_count) = beh.Age_months(find((beh.SubjectID == str2num(grp_contents(s).name(5:7)))));
            
            % Get training group.
            group(sub_count) = beh.group_age(find((beh.SubjectID == str2num(grp_contents(s).name(5:7)))));
            
        end % if ~isempty
        
    end % end if exist
    
end % end s

% Fix group so that it is just children and adult.
group(find(group == 2)) = 1; % change older children to "children = 1""
group(find(group == 3)) = 2; % then change adults to "adults = 2"

meanmotion = mean(motion, 2);
sd = std(motion, 0, 2);

% Concatenate and sort according to group.
toplot = cat(2, subID', group', meanmotion, sd);
toplot = sortrows(toplot, [2 3 1], 'ascend');

% cheap trick
clear subID group m sd b0
subID = toplot(:, 1)';
group = toplot(:, 2)';
fd = toplot(:, 3)';
sd = toplot(:, 4)';

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
        
    elseif group(p) == 2 
        
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
plot([2 2], [0 length(subID)+0.5], ':k')

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

print(fullfile(rootDir, 'plots-singleshell', 'plot_fd'), '-dpng')

hold off;