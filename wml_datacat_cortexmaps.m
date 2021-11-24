% Gathers all data of interest for ping study and produces .mat and .csv
% files that contain the data of interest.

% tractsofinterest = {'leftpArc', 'leftMDLFspl', 'leftMDLFang', 'leftTPC', 'leftIFOF', 'leftILF', 'leftSLF1And2', 'leftSLF3', 'leftVOF', 'leftAslant', ...
%     'rightpArc', 'rightMDLFspl', 'rightMDLFang', 'rightTPC', 'rightIFOF', 'rightILF', 'rightSLF1And2', 'rightSLF3', 'rightVOF', 'rightAslant'};

clear all; close all; clc;

rootdir = '/Volumes/Seagate/wml/wml-wmpredictslearning';
blprojectid = 'proj-61609744fc8eb92cbcf82a5d';
mp2rage = 'no';

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

if strcmp(mp2rage, 'no')
    for f = 1:length(subfolders)
        
        % Grab demographics data.
        d.subID(f) = str2num(subfolders(f).name(end-8:end-6));
        %     d.subID{f} = tdem.subID(find(ismember(tdem.subID, subfolders(f).name(end-4:end))));
        %     d.age(f) = tdem.age_mri(find(ismember(tdem.subID, subfolders(f).name(end-4:end))))/12;
        %     temp = tdem.sex(find(ismember(tdem.subID, subfolders(f).name(end-4:end))));
        %     if strcmp(temp, 'F')
        %         d.sex(f) = 1;
        %     elseif strcmp(temp, 'M')
        %         d.sex(f) = 2;
        %     end
        clear temp
        
        % Read in cortex mapping data: dwi.
        cmfolders_dwi = dir(fullfile(subfolders(f).folder, subfolders(f).name, '*parc-stats.tag-cortex_mapping*'));
        tcm = readtable(fullfile(cmfolders_dwi.folder, cmfolders_dwi.name, '/parc-stats/tracts_MEAN.csv'));
        
        %fa
        d.leftparc_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftpArc_gaussian_3mm_LPI_FiberEndpoint')));
        d.leftparc_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftpArc_gaussian_3mm_RAS_FiberEndpoint')));
        d.rightparc_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightpArc_gaussian_3mm_LPI_FiberEndpoint')));
        d.rightparc_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightpArc_gaussian_3mm_RAS_FiberEndpoint')));
        
        d.lefttpc_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftTPC_gaussian_3mm_LPI_FiberEndpoint')));
        d.lefttpc_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftTPC_gaussian_3mm_RAS_FiberEndpoint')));
        d.righttpc_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightTPC_gaussian_3mm_LPI_FiberEndpoint')));
        d.righttpc_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightTPC_gaussian_3mm_RAS_FiberEndpoint')));
        
        d.leftmdlfspl_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftMDLFspl_gaussian_3mm_LPI_FiberEndpoint')));
        d.leftmdlfspl_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftMDLFspl_gaussian_3mm_RAS_FiberEndpoint')));
        d.rightmdlfspl_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightMDLFspl_gaussian_3mm_LPI_FiberEndpoint')));
        d.rightmdlfspl_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightMDLFspl_gaussian_3mm_RAS_FiberEndpoint')));
        
        d.leftmdlfang_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftMDLFang_gaussian_3mm_LPI_FiberEndpoint')));
        d.leftmdlfang_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftMDLFang_gaussian_3mm_RAS_FiberEndpoint')));
        d.rightmdlfang_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightMDLFang_gaussian_3mm_LPI_FiberEndpoint')));
        d.rightmdlfang_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightMDLFang_gaussian_3mm_RAS_FiberEndpoint')));
        
        d.leftifof_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftIFOF_gaussian_3mm_LPI_FiberEndpoint')));
        d.leftifof_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftIFOF_gaussian_3mm_RAS_FiberEndpoint')));
        d.rightifof_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightIFOF_gaussian_3mm_LPI_FiberEndpoint')));
        d.rightifof_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightIFOF_gaussian_3mm_RAS_FiberEndpoint')));
        
        d.leftilf_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftILF_gaussian_3mm_LPI_FiberEndpoint')));
        d.leftilf_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftILF_gaussian_3mm_RAS_FiberEndpoint')));
        d.rightilf_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightILF_gaussian_3mm_LPI_FiberEndpoint')));
        d.rightilf_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightILF_gaussian_3mm_RAS_FiberEndpoint')));
        
        d.leftslf1and2_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftSLF1And2_gaussian_3mm_LPI_FiberEndpoint')));
        d.leftslf1and2_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftSLF1And2_gaussian_3mm_RAS_FiberEndpoint')));
        d.rightslf1and2_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightSLF1And2_gaussian_3mm_LPI_FiberEndpoint')));
        d.rightslf1and2_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightSLF1And2_gaussian_3mm_RAS_FiberEndpoint')));
        
        d.leftslf3_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftSLF3_gaussian_3mm_LPI_FiberEndpoint')));
        d.leftslf3_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftSLF3_gaussian_3mm_RAS_FiberEndpoint')));
        d.rightslf3_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightSLF3_gaussian_3mm_LPI_FiberEndpoint')));
        d.rightslf3_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightSLF3_gaussian_3mm_RAS_FiberEndpoint')));
        
        d.leftfat_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftAslant_gaussian_3mm_LPI_FiberEndpoint')));
        d.leftfat_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftAslant_gaussian_3mm_RAS_FiberEndpoint')));
        d.rightfat_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightAslant_gaussian_3mm_LPI_FiberEndpoint')));
        d.rightfat_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightAslant_gaussian_3mm_RAS_FiberEndpoint')));
        
        d.leftvof_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftVOF_gaussian_3mm_LPI_FiberEndpoint')));
        d.leftvof_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftVOF_gaussian_3mm_RAS_FiberEndpoint')));
        d.rightvof_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightVOF_gaussian_3mm_LPI_FiberEndpoint')));
        d.rightvof_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightVOF_gaussian_3mm_RAS_FiberEndpoint')));
        
        d.leftuncinate_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftUncinate_gaussian_3mm_LPI_FiberEndpoint')));
        d.leftuncinate_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftUncinate_gaussian_3mm_RAS_FiberEndpoint')));
        d.rightuncinate_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightUncinate_gaussian_3mm_LPI_FiberEndpoint')));
        d.rightuncinate_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightUncinate_gaussian_3mm_RAS_FiberEndpoint')));
        
    end
    
    data_cortexmaps_mean = d;
    
    writetable(d, fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', 'WML_mri_data_cortexmaps_nomp2rage_mean.csv'))
    save(fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', 'WML_mri_data_cortexmaps_nomp2rage.mat'), 'data_cortexmaps_mean')
end