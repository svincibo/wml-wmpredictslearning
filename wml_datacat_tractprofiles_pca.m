% Gathers all data of interest for ping study and produces .mat and .csv
% files that contain the data of interest.

% tractsofinterest = {'leftpArc', 'leftMDLFspl', 'leftMDLFang', 'leftTPC', 'leftIFOF', 'leftILF', 'leftSLF1And2', 'leftSLF3', 'leftVOF', 'leftAslant', ...
%     'rightpArc', 'rightMDLFspl', 'rightMDLFang', 'rightTPC', 'rightIFOF', 'rightILF', 'rightSLF1And2', 'rightSLF3', 'rightVOF', 'rightAslant'};

clear all; close all; clc;

rootdir = '/Volumes/Seagate/wml/wml-wmpredictslearning';
blprojectid = 'proj-61609744fc8eb92cbcf82a5d';

% % Read in demographics data.
% tdem = readtable(fullfile(rootdir, 'supportFiles/ping_combined_participantinformation.csv'));

% Find all subject folders in the blproject directory.
subfolders = dir(fullfile(rootdir, blprojectid));

% Remove the '.' and '..' files.
subfolders = subfolders(arrayfun(@(x) x.name(1), subfolders) ~= '.');

% Keep only names that are subject folders.
subfolders = subfolders(arrayfun(@(x) x.name(1), subfolders) == 's');

% Initiate data table.
d = table();

for s = 1:length(subfolders)
    
    % Grab demographics data.
    subID = str2num(subfolders(s).name(end-8:end-6));
    
    % Read in tractprofiles data: fa.
    proffolders = dir(fullfile(subfolders(s).folder, subfolders(s).name, '*tractmeasures*'));
    ttck = readtable(fullfile(proffolders.folder, proffolders.name, 'tractmeasures.csv'));
    
    % Sort by structureID and then by nodeID.
    ttck = sortrows(ttck, [2 3]);
    
    % Get tractnames.
    tractnames = unique(ttck.structureID);
    
    for t = 1:length(tractnames)
        
        idx = find(contains(ttck.structureID, tractnames{t}));
        
        profile = ttck.fa(idx(21:end-20));
        
        out(:, t) = profile;
        
        clear idx; clear profile;
        
    end % tractnames
    
    sub_profiles_out = array2table(out, 'VariableNames', tractnames);
    
    writetable(sub_profiles_out, fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', 'tractprofiles', ['WML_mri_data_tractprofiles_sub' num2str(subID) '.csv']))
    save(fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', 'tractprofiles', ['WML_mri_data_tractprofiles_sub' num2str(subID) '.mat']), 'sub_profiles_out')
    
    clear out;
    
end

