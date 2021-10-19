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
    
    % Read in tractprofiles data: fa and md and t1/t2 ratio
    proffolders = dir(fullfile(subfolders(f).folder, subfolders(f).name, '*tractmeasures*'));
    ttck = readtable(fullfile(proffolders.folder, proffolders.name, 'tractmeasures.csv'));
    
    idx = find(contains(ttck.structureID, 'leftpArc'));
    d.leftparc_fa(f) = nanmean(ttck.fa(idx(20:end-20))); d.leftparc_md(f) = nanmean(ttck.md(idx(20:end-20)));
    d.leftparc_t1t2(f) = nanmean(ttck.map(idx(20:end-20))); 
%     d.leftparc_T1(f) = nanmean(ttck.T1(idx(20:end-20))); d.leftparc_R1(f) = nanmean(ttck.R1(idx(20:end-20)));
    d.leftparc_ndi(f) = nanmean(ttck.ndi(idx(20:end-20))); d.leftparc_odi(f) = nanmean(ttck.odi(idx(20:end-20))); d.leftparc_isovf(f) = nanmean(ttck.isovf(idx(20:end-20)));
    idx = find(contains(ttck.structureID, 'rightpArc'));
    d.rightparc_fa(f) = nanmean(ttck.fa(idx(20:end-20))); d.rightparc_md(f) = nanmean(ttck.md(idx(20:end-20)));
    d.rightparc_t1t2(f) = nanmean(ttck.map(idx(20:end-20)));
%     d.rightparc_T1(f) = nanmean(ttck.T1(idx(20:end-20))); d.rightparc_R1(f) = nanmean(ttck.R1(idx(20:end-20)));
    d.rightparc_ndi(f) = nanmean(ttck.ndi(idx(20:end-20))); d.rightparc_odi(f) = nanmean(ttck.odi(idx(20:end-20))); d.rightparc_isovf(f) = nanmean(ttck.isovf(idx(20:end-20)));
    idx = find(contains(ttck.structureID, 'leftTPC'));
    d.lefttpc_fa(f) = nanmean(ttck.fa(idx(20:end-20))); d.lefttpc_md(f) = nanmean(ttck.md(idx(20:end-20)));
    d.lefttpc_t1t2(f) = nanmean(ttck.map(idx(20:end-20)));
%         d.lefttpc_T1(f) = nanmean(ttck.T1(idx(20:end-20))); d.lefttpc_R1(f) = nanmean(ttck.R1(idx(20:end-20)));
    d.lefttpc_ndi(f) = nanmean(ttck.ndi(idx(20:end-20))); d.lefttpc_odi(f) = nanmean(ttck.odi(idx(20:end-20))); d.lefttpc_isovf(f) = nanmean(ttck.isovf(idx(20:end-20)));
    
    idx = find(contains(ttck.structureID, 'rightTPC'));
    d.righttpc_fa(f) = nanmean(ttck.fa(idx(20:end-20))); d.righttpc_md(f) = nanmean(ttck.md(idx(20:end-20)));
    d.righttpc_t1t2(f) = nanmean(ttck.map(idx(20:end-20)));
%             d.righttpc_T1(f) = nanmean(ttck.T1(idx(20:end-20))); d.righttpc_R1(f) = nanmean(ttck.R1(idx(20:end-20)));
    d.righttpc_ndi(f) = nanmean(ttck.ndi(idx(20:end-20))); d.righttpc_odi(f) = nanmean(ttck.odi(idx(20:end-20))); d.righttpc_isovf(f) = nanmean(ttck.isovf(idx(20:end-20)));
    
    idx = find(contains(ttck.structureID, 'leftMDLFspl'));
    d.leftmdlfspl_fa(f) = nanmean(ttck.fa(idx(20:end-20))); d.leftmdlfspl_md(f) = nanmean(ttck.md(idx(20:end-20)));
    d.leftmdlfspl_t1t2(f) = nanmean(ttck.map(idx(20:end-20)));
%             d.leftmdlfspl_T1(f) = nanmean(ttck.T1(idx(20:end-20))); d.leftmdlfspl_R1(f) = nanmean(ttck.R1(idx(20:end-20)));
    d.leftmdlfspl_ndi(f) = nanmean(ttck.ndi(idx(20:end-20))); d.leftmdlfspl_odi(f) = nanmean(ttck.odi(idx(20:end-20))); d.leftmdlfspl_isovf(f) = nanmean(ttck.isovf(idx(20:end-20)));
    
    idx = find(contains(ttck.structureID, 'rightMDLFspl'));
    d.rightmdlfspl_fa(f) = nanmean(ttck.fa(idx(20:end-20))); d.rightmdlfspl_md(f) = nanmean(ttck.md(idx(20:end-20)));
    d.rightmdlfspl_t1t2(f) = nanmean(ttck.map(idx(20:end-20)));
%             d.rightmdlfspl_T1(f) = nanmean(ttck.T1(idx(20:end-20))); d.rightmdlfspl_R1(f) = nanmean(ttck.R1(idx(20:end-20)));
    d.rightmdlfspl_ndi(f) = nanmean(ttck.ndi(idx(20:end-20))); d.rightmdlfspl_odi(f) = nanmean(ttck.odi(idx(20:end-20))); d.rightmdlfspl_isovf(f) = nanmean(ttck.isovf(idx(20:end-20)));
    
    idx = find(contains(ttck.structureID, 'leftMDLFang'));
    d.leftmdlfang_fa(f) = nanmean(ttck.fa(idx(20:end-20))); d.leftmdlfang_md(f) = nanmean(ttck.md(idx(20:end-20)));
    d.leftmdlfang_t1t2(f) = nanmean(ttck.map(idx(20:end-20)));
%             d.leftmdlfang_T1(f) = nanmean(ttck.T1(idx(20:end-20))); d.leftmdlfang_R1(f) = nanmean(ttck.R1(idx(20:end-20)));
    d.leftmdlfang_ndi(f) = nanmean(ttck.ndi(idx(20:end-20))); d.leftmdlfang_odi(f) = nanmean(ttck.odi(idx(20:end-20))); d.leftmdlfang_isovf(f) = nanmean(ttck.isovf(idx(20:end-20)));
    
    idx = find(contains(ttck.structureID, 'rightMDLFang'));
    d.rightmdlfang_fa(f) = nanmean(ttck.fa(idx(20:end-20))); d.rightmdlfang_md(f) = nanmean(ttck.md(idx(20:end-20)));
    d.rightmdlfang_t1t2(f) = nanmean(ttck.map(idx(20:end-20)));
%             d.rightmdlfang_T1(f) = nanmean(ttck.T1(idx(20:end-20))); d.rightmdlfang_R1(f) = nanmean(ttck.R1(idx(20:end-20)));
    d.rightmdlfang_ndi(f) = nanmean(ttck.ndi(idx(20:end-20))); d.rightmdlfang_odi(f) = nanmean(ttck.odi(idx(20:end-20))); d.rightmdlfang_isovf(f) = nanmean(ttck.isovf(idx(20:end-20)));
    
    
    idx = find(contains(ttck.structureID, 'leftIFOF'));
    d.leftifof_fa(f) = nanmean(ttck.fa(idx(20:end-20))); d.leftifof_md(f) = nanmean(ttck.md(idx(20:end-20)));
    d.leftifof_t1t2(f) = nanmean(ttck.map(idx(20:end-20)));
%             d.leftifof_T1(f) = nanmean(ttck.T1(idx(20:end-20))); d.leftifof_R1(f) = nanmean(ttck.R1(idx(20:end-20)));
    d.leftifof_ndi(f) = nanmean(ttck.ndi(idx(20:end-20))); d.leftifof_odi(f) = nanmean(ttck.odi(idx(20:end-20))); d.leftifof_isovf(f) = nanmean(ttck.isovf(idx(20:end-20)));
    
    idx = find(contains(ttck.structureID, 'rightIFOF'));
    d.rightifof_fa(f) = nanmean(ttck.fa(idx(20:end-20))); d.rightifof_md(f) = nanmean(ttck.md(idx(20:end-20)));
    d.rightifof_t1t2(f) = nanmean(ttck.map(idx(20:end-20)));
%             d.rightifof_T1(f) = nanmean(ttck.T1(idx(20:end-20))); d.rightifof_R1(f) = nanmean(ttck.R1(idx(20:end-20)));
    d.rightifof_ndi(f) = nanmean(ttck.ndi(idx(20:end-20))); d.rightifof_odi(f) = nanmean(ttck.odi(idx(20:end-20))); d.rightifof_isovf(f) = nanmean(ttck.isovf(idx(20:end-20)));
    
    idx = find(contains(ttck.structureID, 'leftILF'));
    d.leftilf_fa(f) = nanmean(ttck.fa(idx(20:end-20))); d.leftilf_md(f) = nanmean(ttck.md(idx(20:end-20)));
    d.leftilf_t1t2(f) = nanmean(ttck.map(idx(20:end-20)));
%             d.leftilf_T1(f) = nanmean(ttck.T1(idx(20:end-20))); d.leftilf_R1(f) = nanmean(ttck.R1(idx(20:end-20)));
    d.leftilf_ndi(f) = nanmean(ttck.ndi(idx(20:end-20))); d.leftilf_odi(f) = nanmean(ttck.odi(idx(20:end-20))); d.leftilf_isovf(f) = nanmean(ttck.isovf(idx(20:end-20)));
    
    idx = find(contains(ttck.structureID, 'rightILF'));
    d.rightilf_fa(f) = nanmean(ttck.fa(idx(20:end-20))); d.rightilf_md(f) = nanmean(ttck.md(idx(20:end-20)));
    d.rightilf_t1t2(f) = nanmean(ttck.map(idx(20:end-20)));
%             d.rightilf_T1(f) = nanmean(ttck.T1(idx(20:end-20))); d.rightilf_R1(f) = nanmean(ttck.R1(idx(20:end-20)));
    d.rightilf_ndi(f) = nanmean(ttck.ndi(idx(20:end-20))); d.rightilf_odi(f) = nanmean(ttck.odi(idx(20:end-20))); d.rightilf_isovf(f) = nanmean(ttck.isovf(idx(20:end-20)));
    
    
    idx = find(contains(ttck.structureID, 'leftSLF1And2'));
    d.leftslf1and2_fa(f) = nanmean(ttck.fa(idx(20:end-20))); d.leftslf1and2_md(f) = nanmean(ttck.md(idx(20:end-20)));
    d.leftslf1and2_t1t2(f) = nanmean(ttck.map(idx(20:end-20)));
%             d.leftslf1and2_T1(f) = nanmean(ttck.T1(idx(20:end-20))); d.leftslf1and2_R1(f) = nanmean(ttck.R1(idx(20:end-20)));
    d.leftslf1and2_ndi(f) = nanmean(ttck.ndi(idx(20:end-20))); d.leftslf1and2_odi(f) = nanmean(ttck.odi(idx(20:end-20))); d.leftslf1and2_isovf(f) = nanmean(ttck.isovf(idx(20:end-20)));
    
    idx = find(contains(ttck.structureID, 'rightSLF1And2'));
    d.rightslf1and2_fa(f) = nanmean(ttck.fa(idx(20:end-20))); d.rightslf1and2_md(f) = nanmean(ttck.md(idx(20:end-20)));
    d.rightslf1and2_t1t2(f) = nanmean(ttck.map(idx(20:end-20)));
%             d.rightslf1and2_T1(f) = nanmean(ttck.T1(idx(20:end-20))); d.rightslf1and2_R1(f) = nanmean(ttck.R1(idx(20:end-20)));
    d.rightslf1and2_ndi(f) = nanmean(ttck.ndi(idx(20:end-20))); d.rightslf1and2_odi(f) = nanmean(ttck.odi(idx(20:end-20))); d.rightslf1and2_isovf(f) = nanmean(ttck.isovf(idx(20:end-20)));
    
    idx = find(contains(ttck.structureID, 'leftSLF3'));
    d.leftslf3_fa(f) = nanmean(ttck.fa(idx(20:end-20))); d.leftslf3_md(f) = nanmean(ttck.md(idx(20:end-20)));
    d.leftslf3_t1t2(f) = nanmean(ttck.map(idx(20:end-20)));
%             d.leftslf3_T1(f) = nanmean(ttck.T1(idx(20:end-20))); d.leftslf3_R1(f) = nanmean(ttck.R1(idx(20:end-20)));
    d.leftslf3_ndi(f) = nanmean(ttck.ndi(idx(20:end-20))); d.leftslf3_odi(f) = nanmean(ttck.odi(idx(20:end-20))); d.leftslf3_isovf(f) = nanmean(ttck.isovf(idx(20:end-20)));
    
    idx = find(contains(ttck.structureID, 'rightSLF3'));
    d.rightslf3_fa(f) = nanmean(ttck.fa(idx(20:end-20))); d.rightslf3_md(f) = nanmean(ttck.md(idx(20:end-20)));
    d.rightslf3_t1t2(f) = nanmean(ttck.map(idx(20:end-20)));
%             d.rightslf3_T1(f) = nanmean(ttck.T1(idx(20:end-20))); d.rightslf3_R1(f) = nanmean(ttck.R1(idx(20:end-20)));
    d.rightslf3_ndi(f) = nanmean(ttck.ndi(idx(20:end-20))); d.rightslf3_odi(f) = nanmean(ttck.odi(idx(20:end-20))); d.rightslf3_isovf(f) = nanmean(ttck.isovf(idx(20:end-20)));
    
    
    idx = find(contains(ttck.structureID, 'leftAslant'));
    d.leftfat_fa(f) = nanmean(ttck.fa(idx(20:end-20))); d.leftfat_md(f) = nanmean(ttck.md(idx(20:end-20)));
    d.leftfat_t1t2(f) = nanmean(ttck.map(idx(20:end-20)));
%             d.leftfat_T1(f) = nanmean(ttck.T1(idx(20:end-20))); d.leftfat_R1(f) = nanmean(ttck.R1(idx(20:end-20)));
    d.leftfat_ndi(f) = nanmean(ttck.ndi(idx(20:end-20))); d.leftfat_odi(f) = nanmean(ttck.odi(idx(20:end-20))); d.leftfat_isovf(f) = nanmean(ttck.isovf(idx(20:end-20)));
    
    idx = find(contains(ttck.structureID, 'rightAslant'));
    d.rightfat_fa(f) = nanmean(ttck.fa(idx(20:end-20))); d.rightfat_md(f) = nanmean(ttck.md(idx(20:end-20)));
    d.rightfat_t1t2(f) = nanmean(ttck.map(idx(20:end-20)));
%             d.rightfat_T1(f) = nanmean(ttck.T1(idx(20:end-20))); d.rightfat_R1(f) = nanmean(ttck.R1(idx(20:end-20)));
    d.rightfat_ndi(f) = nanmean(ttck.ndi(idx(20:end-20))); d.rightfat_odi(f) = nanmean(ttck.odi(idx(20:end-20))); d.rightfat_isovf(f) = nanmean(ttck.isovf(idx(20:end-20)));
    
    idx = find(contains(ttck.structureID, 'leftVOF'));
    d.leftvof_fa(f) = nanmean(ttck.fa(idx(20:end-20))); d.leftvof_md(f) = nanmean(ttck.md(idx(20:end-20)));
    d.leftvof_t1t2(f) = nanmean(ttck.map(idx(20:end-20)));
%             d.leftvof_T1(f) = nanmean(ttck.T1(idx(20:end-20))); d.leftvof_R1(f) = nanmean(ttck.R1(idx(20:end-20)));
    d.leftvof_ndi(f) = nanmean(ttck.ndi(idx(20:end-20))); d.leftvof_odi(f) = nanmean(ttck.odi(idx(20:end-20))); d.leftvof_isovf(f) = nanmean(ttck.isovf(idx(20:end-20)));
    
    idx = find(contains(ttck.structureID, 'rightVOF'));
    d.rightvof_fa(f) = nanmean(ttck.fa(idx(20:end-20))); d.rightvof_md(f) = nanmean(ttck.md(idx(20:end-20)));
    d.rightvof_t1t2(f) = nanmean(ttck.map(idx(20:end-20)));
%             d.rightvof_T1(f) = nanmean(ttck.T1(idx(20:end-20))); d.rightvof_R1(f) = nanmean(ttck.R1(idx(20:end-20)));
    d.rightvof_ndi(f) = nanmean(ttck.ndi(idx(20:end-20))); d.rightvof_odi(f) = nanmean(ttck.odi(idx(20:end-20))); d.rightvof_isovf(f) = nanmean(ttck.isovf(idx(20:end-20)));
    
    
    idx = find(contains(ttck.structureID, 'leftUncinate'));
    d.leftuncinate_fa(f) = nanmean(ttck.fa(idx(20:end-20))); d.leftuncinate_md(f) = nanmean(ttck.md(idx(20:end-20)));
    d.leftuncinate_t1t2(f) = nanmean(ttck.map(idx(20:end-20)));
%             d.leftuncinate_T1(f) = nanmean(ttck.T1(idx(20:end-20))); d.leftuncinate_R1(f) = nanmean(ttck.R1(idx(20:end-20)));
    d.leftuncinate_ndi(f) = nanmean(ttck.ndi(idx(20:end-20))); d.leftuncinate_odi(f) = nanmean(ttck.odi(idx(20:end-20))); d.leftuncinate_isovf(f) = nanmean(ttck.isovf(idx(20:end-20)));
    
    idx = find(contains(ttck.structureID, 'rightUncinate'));
    d.rightuncinate_fa(f) = nanmean(ttck.fa(idx(20:end-20))); d.rightuncinate_md(f) = nanmean(ttck.md(idx(20:end-20)));
    d.rightuncinate_t1t2(f) = nanmean(ttck.map(idx(20:end-20)));
%             d.rightuncinate_T1(f) = nanmean(ttck.T1(idx(20:end-20))); d.rightuncinate_R1(f) = nanmean(ttck.R1(idx(20:end-20)));
    d.rightuncinate_ndi(f) = nanmean(ttck.ndi(idx(20:end-20))); d.rightuncinate_odi(f) = nanmean(ttck.odi(idx(20:end-20))); d.rightuncinate_isovf(f) = nanmean(ttck.isovf(idx(20:end-20)));
    
    
%     % Read in cortex mapping data: dwi.
%     cmfolders_dwi = dir(fullfile(subfolders(f).folder, subfolders(f).name, '*parc-stats.tag-cortex_mapping*'));
%     tcm = readtable(fullfile(cmfolders_dwi.folder, cmfolders_dwi.name, '/tracts_MEAN.csv'));
    
%     %fa
%     d.leftparc_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftpArc_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftparc_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftpArc_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightparc_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightpArc_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightparc_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightpArc_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.lefttpc_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftTPC_gaussian_5mm_LPI_FiberEndpoint')));
%     d.lefttpc_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftTPC_gaussian_5mm_RAS_FiberEndpoint')));
%     d.righttpc_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightTPC_gaussian_5mm_LPI_FiberEndpoint')));
%     d.righttpc_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightTPC_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.leftmdlfspl_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftMDLFspl_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftmdlfspl_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftMDLFspl_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightmdlfspl_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightMDLFspl_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightmdlfspl_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightMDLFspl_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.leftmdlfang_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftMDLFang_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftmdlfang_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftMDLFang_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightmdlfang_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightMDLFang_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightmdlfang_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightMDLFang_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.leftifof_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftIFOF_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftifof_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftIFOF_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightifof_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightIFOF_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightifof_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightIFOF_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.leftilf_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftILF_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftilf_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftILF_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightilf_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightILF_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightilf_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightILF_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.leftslf1and2_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftSLF1And2_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftslf1and2_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftSLF1And2_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightslf1and2_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightSLF1And2_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightslf1and2_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightSLF1And2_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.leftslf3_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftSLF3_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftslf3_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftSLF3_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightslf3_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightSLF3_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightslf3_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightSLF3_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.leftfat_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftAslant_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftfat_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftAslant_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightfat_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightAslant_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightfat_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightAslant_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.leftvof_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftVOF_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftvof_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftVOF_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightvof_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightVOF_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightvof_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightVOF_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.leftuncinate_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftUncinate_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftuncinate_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'lh.leftUncinate_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightuncinate_ventral_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightUncinate_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightuncinate_dorsal_fa(f) = tcm.fa(find(ismember(tcm.structureID, 'rh.rightUncinate_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     % Read in cortex mapping data: t1t2.
%     cmfolders_t1t2 = dir(fullfile(subfolders(f).folder, subfolders(f).name, '*parc-stats.tag-myelin_mapping*'));
%     tcm = readtable(fullfile(cmfolders_t1t2.folder, cmfolders_t1t2.name, '/tracts_MEAN.csv'));
%     
%     %myelinmap
%     d.leftparc_ventral_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'lh.leftpArc_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftparc_dorsal_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'lh.leftpArc_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightparc_ventral_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'rh.rightpArc_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightparc_dorsal_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'rh.rightpArc_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.lefttpc_ventral_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'lh.leftTPC_gaussian_5mm_LPI_FiberEndpoint')));
%     d.lefttpc_dorsal_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'lh.leftTPC_gaussian_5mm_RAS_FiberEndpoint')));
%     d.righttpc_ventral_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'rh.rightTPC_gaussian_5mm_LPI_FiberEndpoint')));
%     d.righttpc_dorsal_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'rh.rightTPC_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.leftmdlfspl_ventral_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'lh.leftMDLFspl_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftmdlfspl_dorsal_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'lh.leftMDLFspl_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightmdlfspl_ventral_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'rh.rightMDLFspl_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightmdlfspl_dorsal_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'rh.rightMDLFspl_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.leftifof_ventral_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'lh.leftIFOF_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftifof_dorsal_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'lh.leftIFOF_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightifof_ventral_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'rh.rightIFOF_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightifof_dorsal_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'rh.rightIFOF_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.leftilf_ventral_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'lh.leftILF_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftilf_dorsal_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'lh.leftILF_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightilf_ventral_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'rh.rightILF_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightilf_dorsal_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'rh.rightILF_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.leftslf1and2_ventral_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'lh.leftSLF1And2_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftslf1and2_dorsal_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'lh.leftSLF1And2_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightslf1and2_ventral_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'rh.rightSLF1And2_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightslf1and2_dorsal_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'rh.rightSLF1And2_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.leftslf3_ventral_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'lh.leftSLF3_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftslf3_dorsal_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'lh.leftSLF3_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightslf3_ventral_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'rh.rightSLF3_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightslf3_dorsal_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'rh.rightSLF3_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.leftfat_ventral_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'lh.leftAslant_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftfat_dorsal_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'lh.leftAslant_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightfat_ventral_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'rh.rightAslant_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightfat_dorsal_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'rh.rightAslant_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.leftvof_ventral_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'lh.leftVOF_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftvof_dorsal_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'lh.leftVOF_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightvof_ventral_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'rh.rightVOF_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightvof_dorsal_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'rh.rightVOF_gaussian_5mm_RAS_FiberEndpoint')));
%            
%     d.leftuncinate_ventral_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'lh.leftUncinate_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftuncinate_dorsal_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'lh.leftUncinate_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightuncinate_ventral_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'rh.rightUncinate_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightuncinate_dorsal_t1t2(f) = tcm.myelinmap(find(ismember(tcm.structureID, 'rh.rightUncinate_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     %thickness
%     d.leftparc_ventral_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'lh.leftpArc_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftparc_dorsal_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'lh.leftpArc_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightparc_ventral_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'rh.rightpArc_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightparc_dorsal_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'rh.rightpArc_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.lefttpc_ventral_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'lh.leftTPC_gaussian_5mm_LPI_FiberEndpoint')));
%     d.lefttpc_dorsal_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'lh.leftTPC_gaussian_5mm_RAS_FiberEndpoint')));
%     d.righttpc_ventral_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'rh.rightTPC_gaussian_5mm_LPI_FiberEndpoint')));
%     d.righttpc_dorsal_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'rh.rightTPC_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.leftmdlfspl_ventral_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'lh.leftMDLFspl_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftmdlfspl_dorsal_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'lh.leftMDLFspl_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightmdlfspl_ventral_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'rh.rightMDLFspl_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightmdlfspl_dorsal_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'rh.rightMDLFspl_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.leftmdlfang_ventral_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'lh.leftMDLFang_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftmdlfang_dorsal_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'lh.leftMDLFang_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightmdlfang_ventral_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'rh.rightMDLFang_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightmdlfang_dorsal_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'rh.rightMDLFang_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.leftifof_ventral_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'lh.leftIFOF_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftifof_dorsal_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'lh.leftIFOF_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightifof_ventral_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'rh.rightIFOF_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightifof_dorsal_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'rh.rightIFOF_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.leftilf_ventral_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'lh.leftILF_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftilf_dorsal_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'lh.leftILF_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightilf_ventral_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'rh.rightILF_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightilf_dorsal_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'rh.rightILF_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.leftslf1and2_ventral_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'lh.leftSLF1And2_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftslf1and2_dorsal_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'lh.leftSLF1And2_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightslf1and2_ventral_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'rh.rightSLF1And2_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightslf1and2_dorsal_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'rh.rightSLF1And2_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.leftslf3_ventral_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'lh.leftSLF3_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftslf3_dorsal_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'lh.leftSLF3_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightslf3_ventral_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'rh.rightSLF3_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightslf3_dorsal_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'rh.rightSLF3_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.leftfat_ventral_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'lh.leftAslant_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftfat_dorsal_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'lh.leftAslant_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightfat_ventral_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'rh.rightAslant_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightfat_dorsal_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'rh.rightAslant_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.leftvof_ventral_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'lh.leftVOF_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftvof_dorsal_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'lh.leftVOF_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightvof_ventral_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'rh.rightVOF_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightvof_dorsal_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'rh.rightVOF_gaussian_5mm_RAS_FiberEndpoint')));
%     
%     d.leftuncinate_ventral_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'lh.leftUncinate_gaussian_5mm_LPI_FiberEndpoint')));
%     d.leftuncinate_dorsal_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'lh.leftUncinate_gaussian_5mm_RAS_FiberEndpoint')));
%     d.rightuncinate_ventral_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'rh.rightUncinate_gaussian_5mm_LPI_FiberEndpoint')));
%     d.rightuncinate_dorsal_thickness(f) = tcm.thickness(find(ismember(tcm.structureID, 'rh.rightUncinate_gaussian_5mm_RAS_FiberEndpoint')));
%     
    
end

data_tractprofiles_mean = d;

writetable(d, fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', 'WML_mri_data_tractprofiles_mean.csv'))
save(fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', 'WML_mri_data_tractprofiles.mat'), 'data_tractprofiles_mean')
