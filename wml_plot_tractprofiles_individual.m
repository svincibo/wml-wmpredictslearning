% This script reads in FA and MD measures (from Brad Caron's
% TractProfiles App) for each of the tracts generated (from Dan Bullock's
% White Matter Segmentation App).
% It also reads in behavioral data (e.g., experimental group) collected as part of the
% spade study.
%it runs twice, one for producing the plots for mena Fa and difference (S2-S1) and one for mean MD and
%Md difference (S2-S1)

clear all; close all; clc
format shortG

w_measures = {'fa', 'md', 'ad', 'rd', 'odi', 'ndi', 'isovf', 'R1', 'T1'};
rootdir = '/Volumes/Seagate/wml/wml-wmpredictslearning';
blprojectid = 'proj-61609744fc8eb92cbcf82a5d';

for w = 1:length(w_measures)
    
    measure = w_measures{w};
    
    di_color = [0.6350 0.0780 0.1840]; %red
    
    hold on;
    linewidth = 1.5;
    linestyle = 'none';
    fontname = 'Arial';
    fontsize = 20;
    fontangle = 'italic';
    xticklength = 0;
    
    save_figures = 'yes';
    
    % Should outliers be removed? If so, which subIDs?
    remove_outliers = 'yes';
    if strcmp(remove_outliers, 'yes')
        
        % Identify outliers to be removed - e.g., outlier = [108 126 212 214 318];
        outlier = [];
        
    else
        
        outlier = [];
        
    end
    
    % % Read in behavioral data.
    % beh_tbl = readtable(fullfile(rootdir, 'supportFiles', 'wml_demographics.csv'), 'TreatAsEmpty', {'.', 'na'});
    
    % Read in mri data.
    load(fullfile(rootdir, 'wml-wmpredictslearning-supportFiles', 'wml_data_mri_longform.mat'))
    mri = data_tbl; clear data_tbl;
    
    % Set index for clipping beginning and end of tractprofile.
    idx = 21:length(unique(mri.nodeID))-20;
    mri.clipped = ismember(mri.nodeID, idx); clear idx;
    
    % Select clipped tractprofiles for each session.
    mri_1 = mri(find(mri.session == 1 & mri.clipped == 1), :);
    % mri_2 = mri(find(mri.session == 2 & mri.clipped == 1), :);
    % mri_3 = mri(find(mri.session == 3 & mri.clipped == 1), :);
    
    % Select y-axis limits.
    if strcmp(measure, 'fa')
        ylim_lo = 0.30; ylim_hi = 1;
        ylabel = 'Fractional Anisotropy (FA)';
    elseif strcmp(measure, 'md')
        ylim_lo = 0.5; ylim_hi = 1.5;
        ylabel = 'Mean Diffusivity (MD)';
    elseif strcmp(measure, 'ad')
        ylim_lo = 1.2; ylim_hi = 2.4;
        ylabel = 'Axial Diffusivity (AD)';
    elseif strcmp(measure, 'rd')
        ylim_lo = 0.3; ylim_hi = 0.8;
        ylabel = 'Radial Diffusivity (RD)';
    elseif strcmp(measure, 'odi')
        ylim_lo = 0.10; ylim_hi = 0.30;
        ylabel = 'Orientation Dispersion Index (ODI)';
    elseif strcmp(measure, 'ndi')
        ylim_lo = 0.5; ylim_hi = 0.8;
        ylabel = 'Neurite Density Index (NDI)';
    elseif strcmp(measure, 'isovf')
        ylim_lo = 0; ylim_hi = 0.30;
        ylabel = 'Isometric Volume Fraction (ISOVF)';
    elseif strcmp(measure, 'R1')
        ylim_lo = 0.5; ylim_hi = 1.5;
        ylabel = 'Quantitative MRI (R1)';
    elseif strcmp(measure, 'T1')
        ylim_lo = 0.5; ylim_hi = 1.5;
        ylabel = 'Quantitative MRI (T1)';
    end
    
    %% TRACTOGRAPHY.
    
    list_tracts = unique(mri.structureID);
    list_subs = unique(mri.subjectID);
    
    % Make plot for each tract.
    for t = 1:length(list_tracts)
        
        figure(t);
        hold on;
        
        tractname = list_tracts{t};
        
        % Plot individual data.
        c = colormap(parula); c = c(1:size(c, 1)/length(list_subs):size(c, 1), :);
        for s = 1:length(list_subs)
            
            % Select data for this round: session 1.
            t_idx = find(contains(mri_1.structureID, list_tracts{t}));
            s_idx = find(mri_1.subjectID == list_subs(s));
            m_idx = find(ismember(mri_1.Properties.VariableNames, measure));
            session1 = table2array(mri_1(intersect(t_idx, s_idx), m_idx));
            clear t_idx s_idx m_idx;
            
            %         % Select data for this round: session 2.
            %         t_idx = find(contains(mri_2.structureID, list_tracts{t}));
            %         s_idx = find(mri_2.subjectID == list_subs(s));
            %         m_idx = find(contains(mri_2.Properties.VariableNames, measure));
            %         session2 = table2array(mri_2(intersect(t_idx, s_idx), m_idx));
            %         clear t_idx s_idx m_idx;
            %
            %         % Select data for this round: session 3.
            %         t_idx = find(contains(mri_3.structureID, list_tracts{t}));
            %         s_idx = find(mri_3.subjectID == list_subs(s));
            %         m_idx = find(contains(mri_3.Properties.VariableNames, measure));
            %         session3 = table2array(mri_3(intersect(t_idx, s_idx), m_idx));
            %         clear t_idx s_idx m_idx;
            
            % Session 1
            if ~isempty(session1)
                p1 = plot(1:size(session1), session1);
                if ~isempty(p1)
                    p1.Color = c(s, :); p1.LineStyle = '-';
                end
                clear p1;
            end
            
            %         % Session 2
            %         p2 = plot(1:size(session2), session2);
            %         if ~isempty(p2)
            %             p2.Color = c(s, :); p2.LineStyle = '-.';
            %         end
            %
            %         % Session 3
            %         p3 = plot(1:size(session3), session3);
            %         if ~isempty(p3)
            %             p3.Color = c(s, :); p3.LineStyle = ':';
            %         end
            
        end % end subject
        %
        %     % Plot group mean and standard deviation.
        %     submean = nanmean([session1 session2 session3], 2);
        %     p = plot(1:size(session1), submean);
        %     p.Color = 'k'; p.LineStyle = '-'; p.LineWidth = 3;
        %     hi = nanmean(submean, 2) + 1.96*nanstd(submean, 0, 2)/sqrt(size(~isnan(submean), 2)); lo = nanmean(submean, 2) - 1.96*nanstd(submean, 0, 2)/sqrt(size(~isnan(submean), 2)); x = (1:size(nanmean(submean, 2),1))';
        %     hp1 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], 'k');
        %     set(hp1, 'facecolor', 'k', 'edgecolor', 'none', 'facealpha', .2);
        %     clear submean;
        
        
        g = gca;
        g.YLabel.String = ylabel;
        g.YLabel.FontSize = fontsize;
        
        % xaxis
        xax = get(gca, 'xaxis');
        xax.Limits = [0 160];
        xax.TickValues = [0 80 160];
        xax.TickLabels = {'20', '100', '180'};
        xax.TickDirection = 'out';
        xax.TickLength = [xticklength xticklength];
        xax.FontName = fontname;
        xax.FontSize = fontsize;
        % xax.FontAngle = fontangle;
        
        % yaxis
        yax = get(gca,'yaxis');
        yax.Limits = [ylim_lo ylim_hi];
        yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
        yax.TickDirection = 'out';
        yax.TickLabels = {num2str(ylim_lo, '%2.2f'), num2str((ylim_lo+ylim_hi)/2, '%2.2f'), num2str(ylim_hi, '%2.2f')};
        yax.FontName = fontname;
        yax.FontSize = fontsize;
        
        % general
        g = gca;
        box off
        g.XLabel.String = 'Location along tract';
        g.XLabel.FontSize = fontsize;
        g.XLabel.FontAngle = fontangle;
        pbaspect([1 1 1])
        %             legend ('s1', '', 's2', '',  'Location', 'southoutside')
        
        title(list_tracts{t})
        
        print(fullfile(rootdir, 'wml-wmpredictslearning-plots', ['plot_tractprofiles_indiv_' measure '_' list_tracts{t}]), '-dpng')
        print(fullfile(rootdir, 'wml-wmpredictslearning-plots', 'eps', ['plot_tractprofiles_indiv_' measure '_' list_tracts{t}]), '-depsc')
        
        hold off;
        
    end % end t
    
    close all;
    
end

