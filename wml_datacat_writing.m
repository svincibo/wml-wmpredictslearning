clear all; close all; clc

% Set working directories.
rootDir = '/Volumes/Seagate/wml/';

% Get contents of the directory where the measures are stored.
grp_contents = dir(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-data-beh', '*train_*.txt'));

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

% % Remove the temp not in use files.
% grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= 'n');

% scount = 0;
for s = 1:size(grp_contents, 1)
    
%     % Get filenames of letter/shape data for this subject.
%     sub_contents = dir(fullfile(grp_contents(s).folder, grp_contents(s).name));
%     
%     % Remove the '.' and '..' files.
%     sub_contents = sub_contents(arrayfun(@(x) x.name(1), sub_contents) ~= '.');
    
    if ~isempty(grp_contents(s))
        
%         scount = scount + 1;
        
%         % Grab subID.
%         sub(scount) = str2num(grp_contents(s).name(4:5));
        
        % Display current sub ID.
        disp(grp_contents(s).name)
        
%         % Read in writing data for this subject for each file (i.e., day).
%         for day = 1:size(sub_contents, 1)
            
            data_temp = readtable(fullfile(grp_contents(s).folder, grp_contents(s).name));
            
            % If data_recog exists, append; if not, create.
            if s == 1
                
                % Create data_out array.
                data_write = data_temp;
                
            else
                
                % Concatenate this array with the previous subject's array.
                data_write = cat(1, data_write, data_temp);
                
            end
            
            clear data_temp
            
%         end
        
    end
    
end

% Add header because readtable isn't recognizing the header in the txt files for some reason (Oct 2020).
data_write.Properties.VariableNames = {'subID', 'group', 'day', 'symbolname', 'block', 'trial', 'drawduration', 'trialduration'};

% Find DI observations.
% idx_DI = (data_write.group == 1);
data_write_DI = data_write;%(idx_DI, :);
% Get lower and upper bound.
idx_above = find(data_write_DI.drawduration > (nanmean(data_write_DI.drawduration)+3*nanstd(data_write_DI.drawduration)));
idx_below = find(data_write_DI.drawduration < (nanmean(data_write_DI.drawduration)-3*nanstd(data_write_DI.drawduration)));
% Remove outliers in DI.
data_write_DI(cat(1, idx_above, idx_below), :) = [];

% Recombine.
clear data_write
data_write = cat(1, data_write_DI); %, data_write_DnI);
clear data_write_DI %data_write_DnI

%% Output mean data

% Get individual subject means for each day.
subjectlist = unique(data_write.subID);
for sub = 1:length(subjectlist)
    
    for day = 1:length(unique(data_write.day))
    
        clear idx;
        idx = find(data_write.subID == subjectlist(sub) & data_write.day == day);
        
        if isempty(idx)
            
            drawduration(sub, day) = NaN;
            
        else
            
            % Mean
            drawduration(sub, day) = nanmean(data_write.drawduration(idx));
%             
%             % Estimate learning rate by fitting a double exponential function
%             x = data_write.
% %             learningrate(sub, day) = fit(x, y, 'exp2', 'Robust', 'Lar');
%             
        end
    
    end
    
end
A = cat(2, subjectlist, drawduration);

% Remove NaN columns and rows.
A = A(:,any(~isnan(A)));  % for columns
data_write_mean = array2table(A(any(~isnan(A),2),:));   %for rows
clear A;

% Create date-specific file name that indicates how many subjects.
filename = sprintf('wml_beh_data_write_%s', datestr(now,'yyyymmdd'));

% Save all variables.
save(fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-supportFiles', filename), 'data_write', 'data_write_mean');

% Save as a CSV files.
writetable(data_write, fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-supportFiles', [filename '.csv']))
writetable(data_write_mean, fullfile(rootDir, 'wml-wmpredictslearning', 'wml-wmpredictslearning-supportFiles', [filename '_mean.csv']))


