% This script reads in FA, MD, AD, and RD measures (from Brad Caron's
% TractProfiles App) for each of the tracts generated and concatenates
% those measurements with behavioral measurements.

clear all; close all; clc
format shortG

% remove_outliers = 'yes';

w_measures = {'fa', 'ad', 'md', 'rd', 'ndi', 'isovf', 'odi', 'map', 'T1', 'R1'};
rootdir = '/Volumes/Seagate/wml/wml-wmpredictslearning';
blprojectid = 'proj-61609744fc8eb92cbcf82a5d';

% Each session one at a time.
% for s = 1:3

%% TRACTOGRAPHY.

% Get contents of the directory where the tract measures for this subject are stored.
grp_contents = dir(fullfile(rootdir, blprojectid));

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

%     % Keep only names that are subject folders.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');

% Load in each tract's tractography measures for this subject.
for i = 1:size(grp_contents, 1)-12
    
    % Display current sub ID.
    disp(grp_contents(i).name)
    
    % Get contents of the directory where the tract measures for this subject are stored.
    sub_contents_tractprofiles = dir(fullfile(grp_contents(i).folder, grp_contents(i).name,  '/dt-neuro-tractmeasures.tag-profile*/*.csv'));
    
    % Read in data for this subject and this tract.
    data_temp = readtable(fullfile(sub_contents_tractprofiles.folder, sub_contents_tractprofiles.name));
    
    % Clip beginning and end to avoid partial volume effects.
    data_temp = data_temp(21:end-20, :);
    
    % Append session indicator: just say that all sessions = 1 because it doesn't matter for this analysis.
    data_temp.session = repmat(1, [size(data_temp, 1) 1]);
    
    if i == 1
        
        data = data_temp;
        
    else
        
        data = cat(1, data, data_temp);
        
    end
    
    clear data_temp;
    
end % end i

% end % end s

data_tbl = data;
%
% Save all variables.
save(fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', 'wml_data_mri_longform.mat'), 'data_tbl')
%
% Write out table.
writetable(data_tbl, fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', 'wml_data_mri_longform.csv'));

% Find empty cells.
t = find(cellfun(@isempty, data.structureID));

% Enter 'empty' in empty cells.
data.structureID(t) = {'empty'};

% Get a list of unique subID names.
list_subs = unique(data.subjectID);

% Get a list of unique tract names.
list_tract = unique(data.structureID);

% Get WM measurements for each tract (reorganizing so that each column is a tract).
for w = 1%:length(w_measures)
    
    wm_measure = w_measures{w};
    
    % Find column index for this wm_measure.
    idx = find(strcmp(data.Properties.VariableNames, wm_measure));
    data_here = data(:, idx);
    
    for ses = 1%:3
        
        for sub = 1:size(list_subs, 1)
            %
            %             % Get subID.
            %             subID = list_subs{sub};
            
            for t = 1:size(list_tract, 1)
                
                %                 % Get tractID.
                %                 tractID = list_tract{t};
                
                % Find the row indices for all values that are for this
                % session, this subject, and this tract.
                idx_here = find(data.subjectID == list_subs(sub) & ...
                    strcmp(data.structureID, list_tract{t}));
                
                % Select the wm_measurements for this tract from each subject.
                d = data_here(idx_here, 1);
                
                % Get the mean of the wm_measure for this tract.
                mean_measure = nanmean(table2array(d));
                
                % Append to the output table.
                data_out(sub, t) = mean_measure;
                
                clear temp
                
            end % end k
            
        end % end i
        
        data_out = array2table(data_out);
        data_out = cat(2, array2table(list_subs), data_out);
        data_out.Properties.VariableNames = vertcat('subID', list_tract);
        
        %         % Remove outliers.
        %         if strcmp(remove_outliers, 'yes') && exist('outlier')
        %
        %             % Get index for outliers to be removed.
        %             idx_outlier = ismember(data_all(:, find(strcmp(data_all_header, {'subID'}))), outlier);
        %
        %             % Remove outliers.
        %             data_all = data_all(~idx_outlier, :);
        %
        %         end
        %
        data_tbl = data_out;
        %
        % Save all variables.
        save(fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', ['wml_data_mri_' wm_measure '_shortform.mat']), 'data_tbl')
        %
        % Write out table.
        writetable(data_out, fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', ['wml_data_mri_' wm_measure '_shortform.csv']));
        
        %         % Reset for next loop.
        %         clearvars -except w rootdir beh_data_in_tbl beh_data_in_header beh_data_in blprojectid remove_outliers w_measures outlier
        %
        clear data_out;
        
    end % end ses
    
end % end w