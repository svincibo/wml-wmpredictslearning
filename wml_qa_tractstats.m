% This script reads in streamline count values (i.e., nfibers) and
% checks that the number of streamlines in each tract is correlated across
% subjects between sessions and checks that there are no significant
% differences in streamline count within tracts between groups at either session.
% A box plot is provided to help identify unusually low streamline counts within
% a particular tract and to compare sessions.

clear all; close all; clc
format shortG

rootdir = '/Volumes/Seagate/wml/wml-wmpredictslearning';
blprojectid = 'proj-61609744fc8eb92cbcf82a5d';

% Should outliers be removed? If so, which subIDs?
remove_outliers = 'no';
if strcmp(remove_outliers, 'yes')
     
    % Identify outliers to be removed - liberal removal.
    outlier = [108 116 125 126 203 206 212 214 315 316 318];
    
else
    
    outlier = [];
    
end

a_color = [0 0 0]; %[75 75 75]/255; % gray [.146 0 0]; % light black
c_color = [204 0 204]/255; %pink [178 34 34]/255; % firebrick red [0 .73 .73]; % turquoise

% % Read in behavioral data.
% beh_data_in_tbl = readtable([rootDir 'supportFiles/LWX_all_groupings.csv'], 'TreatAsEmpty', {'.', 'na'});

% Get contents of the directory where the tract measures for this subject are stored.
grp_contents = dir(fullfile(rootdir, blprojectid));

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

% Keep only names that are subject folders.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');

% Load in each tract's tractography measures for this subject.
sub_count = 0;
for t = 1:size(grp_contents, 1)
    
    % Only collect values for subjects that have both MRI and behaviora/demographic data.
%     if ~isempty(find((beh_data_in_tbl.SubjectID == str2num(grp_contents(t).name(5:7)))))
        
        % Get contents of the directory where the tract measures for this subject are stored.
        sub_contents_tractstats = dir(fullfile(grp_contents(t).folder, grp_contents(t).name, '/dt-neuro-tractmeasures.tag-cleaned*/*.csv'));
        
        % Remove the '.' and '..' files.
        sub_contents_tractstats = sub_contents_tractstats(arrayfun(@(x) x.name(1), sub_contents_tractstats) ~= '.');
        
        if ~isempty(sub_contents_tractstats)
            
            % Display current sub ID.
            disp(grp_contents(t).name)
            
            % Update subject counter for when not all subjects are used/needed.
            sub_count = sub_count + 1;
            
            % Read in data for this subject and this tract.
            data_tbl_in = readtable(fullfile(sub_contents_tractstats.folder, sub_contents_tractstats.name));
            
            % Convert into header for ease.
            data_all_in_header = data_tbl_in.structureID;
            
            % Get index matrices for hypothesis-driven grouping of WM tracts.
            for k = 1:length(data_all_in_header)
                
                % Indices of horizontal tracts.
                toi_idx(k) = strcmp(data_all_in_header{k}, 'leftSLF1And2') || strcmp(data_all_in_header{k}, 'rightSLF1And2') ...
                    || strcmp(data_all_in_header{k}, 'leftSLF3') || strcmp(data_all_in_header{k}, 'rightSLF3') ...
                    || strcmp(data_all_in_header{k}, 'leftILF') || strcmp(data_all_in_header{k}, 'rightILF') ...
                    || strcmp(data_all_in_header{k}, 'leftAslant') || strcmp(data_all_in_header{k}, 'rightAslant') ...
                    || strcmp(data_all_in_header{k}, 'leftTPC') || strcmp(data_all_in_header{k}, 'rightTPC') ...
                    || strcmp(data_all_in_header{k}, 'leftpArc') || strcmp(data_all_in_header{k}, 'rightpArc') ...
                    || strcmp(data_all_in_header{k}, 'leftMDLFspl') || strcmp(data_all_in_header{k}, 'rightMDLFspl') ...
                    || strcmp(data_all_in_header{k}, 'leftVOF') || strcmp(data_all_in_header{k}, 'rightVOF') ...
                    || strcmp(data_all_in_header{k}, 'leftMDLFang') || strcmp(data_all_in_header{k}, 'rightMDLFang') ...
                    || strcmp(data_all_in_header{k}, 'leftIFOF') || strcmp(data_all_in_header{k}, 'rightIFOF');
                
            end % end k
            
            % Get positions of tracts of interest in data_tbl.
            toi = find(toi_idx == 1);
            
            for t2 = 1:length(toi)
                
                % Read in mean stat.
                m(sub_count, t2) = data_tbl_in.StreamlineCount(toi(t2));
                
                % Grab tract name for grouping variable.
                tract{sub_count, t2} = char(data_tbl_in.structureID(toi(t2)));
                
                if m(sub_count, t2) <= 100
                    
%                     % Change entry to NaN so that it will not be included.
%                     m(t, t2) = NaN;
                      
                    % Alert user if tract has less than 300 streamlines.
                    disp(['Check data. Number of streamlines is ' num2str(m(sub_count, t2)) ' for ' tract{sub_count, t2} ' in subject ' grp_contents(t).name(5:7) '.'])
                    
                end % end if
                
            end % end t
            
            % Grab subID.
            sub(sub_count) = str2num(grp_contents(t).name(end-8:end-6));
            
%             % Get age group.
%             group(sub_count) = beh_data_in_tbl.group_age(find((beh_data_in_tbl.SubjectID == str2num(grp_contents(t).name(5:7)))));
%             
%             % Get age in months.
%             age(sub_count) = beh_data_in_tbl.group_age(find((beh_data_in_tbl.SubjectID == str2num(grp_contents(t).name(5:7)))));
%             
            clear data_temp toi toi_idx
            
        end % end if empty
        
%     end % end if empty
    
end % end t

% % Fix group so that it is just children and adult.
% group(find(group == 2)) = 1; % change older children to "children = 1""
% group(find(group == 3)) = 2; % then change adults to "adults = 2"

% Find empty cells.
idx = find(cellfun(@isempty,tract));

% Enter 'empty' in empty cells.
tract(idx) = {'empty'};

% Enter NaN for m in empty cells.
m(idx) = NaN;

% Get a list of unique tract names.
list_tract = tract(1, :);
list_tract = list_tract(~strcmp(list_tract, 'empty'));

% % % Determine which subIDs appear in both WM and BEH.
% sub_tract_beh = intersect(sub, beh_data_in_tbl.SubjectID);
% %
% % % Get indices of subjects who appear in both WM and BEH.
% sub_idx_wm = ismember(sub, sub_tract_beh);
% % sub_idx_beh = ismember(beh_data_in_tbl.No, sub_tract_beh);

% Select only subjects who appear in both WM and BEH.
% Concatenate into one data array and one header array.
% Remove redundant subID columns.
data_out = cat(2, sub',  m);
data_out_header = [{'subID'}, list_tract];

% Remove outliers.
if strcmp(remove_outliers, 'yes') && exist('outlier')
    
    % Get index for outliers to be removed.
    idx_outlier = ismember(data_out(:, find(strcmp(data_out_header, {'subID'}))), outlier);
    
    % Remove outliers.
    data_out = data_out(~idx_outlier, :);
    
end

data_tbl = array2table(data_out, 'VariableNames', data_out_header);

% Save all variables.
streamlinecounts = data_tbl;
save(fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', 'wml_wmlpredictslearning_data_streamlinecount.mat'), 'streamlinecounts')

% Write out table.
writetable(data_tbl, fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', 'wml_wmpredictslearning_data_streamlinecount.csv'));

% % Group differences test
% d = m(sub_idx_wm, :);
% for tn = 1:size(d, 2)
%     
% %     disp(list_tract{tn});
%     
%     disp(['Are there streamline count differences between age groups for the ' list_tract{tn} '?'])
%     tracttotest = d(:, tn);
%     [p, tableout, stats] = anova1(tracttotest, group, 'off');
%     disp(['F(' num2str(tableout{2, 3}) ', ' num2str(tableout{3, 3}) ') = ' num2str(tableout{2, 5}) ', p = ' num2str(tableout{2, 6}) '.'])
%     
% end % end tn

% Look at boxplot of streamline count.

% Convert into shorter format.
figure(1)
d = table2array(data_tbl);
boxplot(d(:, 2:end), 'colors', 'k', 'Labels', list_tract, 'PlotStyle', 'compact', 'Widths', 1, 'Orientation', 'horizontal')
hold on
xtickangle(90)
ylabel('Track Name')
xlabel('Streamline Count')
set(gca, 'XScale', 'log')
xlim([1 100000])
box off
print(fullfile(rootdir, 'wml-wmpredictslearning-plots', 'wml_boxplot_streamlinecount'), '-dpng')
print(fullfile(rootdir, 'wml-wmpredictslearning-plots', 'eps', 'wml_boxplot_streamlinecount'), '-depsc')