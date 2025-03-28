function runGroupAnalysisWithInterpolatedData(experiment_path, session_id)
    % RUNGROUPANALYSISWITHINTERPOLATEDDATA - Performs group-level analysis with interpolated normalized ISPC
    %
    % This function creates group averages for ACTIVE and SHAM conditions
    % using the interpolated normalized ISPC values
    %
    % Inputs:
    %   experiment_path - Base path for the experiment data
    %   session_id - Session ID (e.g., 'N1')
    
    % Set default session if not provided
    if nargin < 2
        session_id = 'N1';
    end
    
    fprintf('Running group analysis with interpolated normalized data for session %s...\n', session_id);
    
    % Get base paths
    code_path = '/Users/idohaber/Git-Projects/entrainment';
    src_path = fullfile(code_path, 'src');
    assets_path = fullfile(code_path, 'src/assets');
    
    % Add EEGLAB to the path
    eeglab_path = '/Users/idohaber/Documents/MATLAB/eeglab2024.0/';
    addpath(eeglab_path);
    eeglab nogui;
    
    % Add utility paths
    utilities_path = fullfile(code_path, 'src/utils');
    addpath(utilities_path);
    addpath(src_path);
    
    % Load excluded channels
    exclude_json = fullfile(assets_path, 'exclude.json');
    fprintf('Loading excluded channels from: %s\n', exclude_json);
    
    fid = fopen(exclude_json, 'r');
    if fid == -1
        warning('Could not open exclude.json file. No channels will be excluded.');
        excluded_channels = {};
    else
        try
            excluded_channels = jsondecode(fileread(exclude_json));
            fprintf('Loaded %d excluded channels\n', length(excluded_channels));
        catch err
            fclose(fid);
            warning('Error reading exclude.json file: %s', err.message);
            excluded_channels = {};
        end
        fclose(fid);
    end
    
    % Define electrode file path
    electrode_file = fullfile(code_path, 'src/utils/egi256_GSN_HydroCel.sfp');
    fprintf('Using electrode locations from: %s\n', electrode_file);
    
    % Define group results directory
    group_results_dir = fullfile(experiment_path, 'group_results', session_id, 'normalized_interp');
    if ~exist(group_results_dir, 'dir')
        [success, msg] = mkdir(group_results_dir);
        if ~success
            error('Failed to create group results directory %s: %s', group_results_dir, msg);
        end
        fprintf('Created group results directory: %s\n', group_results_dir);
    end
    
    % Load subject condition assignments
    subject_condition_file = fullfile(assets_path, 'subject_condition.json');
    fprintf('Loading subject condition data from: %s\n', subject_condition_file);
    
    fid = fopen(subject_condition_file, 'r');
    if fid == -1
        error('Could not open subject condition file: %s', subject_condition_file);
    end
    
    try
        condition_data = jsondecode(fileread(subject_condition_file));
        subjects = condition_data.subjects;
    catch err
        fclose(fid);
        error('Error reading subject condition file: %s', err.message);
    end
    fclose(fid);
    
    % Group subjects by condition
    active_subjects = {};
    sham_subjects = {};
    
    for i = 1:length(subjects)
        subject = subjects(i);
        if strcmpi(subject.condition, 'ACTIVE')
            active_subjects{end+1} = subject.id;
        elseif strcmpi(subject.condition, 'SHAM')
            sham_subjects{end+1} = subject.id;
        end
    end
    
    fprintf('Found %d ACTIVE subjects and %d SHAM subjects\n', ...
        length(active_subjects), length(sham_subjects));
    
    % Create group averages for each condition
    createInterpolatedGroupAverage('ACTIVE', active_subjects, session_id, experiment_path, group_results_dir, electrode_file, excluded_channels);
    createInterpolatedGroupAverage('SHAM', sham_subjects, session_id, experiment_path, group_results_dir, electrode_file, excluded_channels);
    
    fprintf('Group analysis with interpolated normalized data completed successfully.\n');
end

function createInterpolatedGroupAverage(condition, subjects, session_id, experiment_path, group_results_dir, electrode_file, excluded_channels)
    % Create average across subjects within a condition using interpolated normalized data
    fprintf('Creating group average for %s condition (%d subjects)...\n', ...
        condition, length(subjects));
    
    % Initialize data structures
    all_norm_topos = struct();
    all_pct_changes = struct();
    all_global_ispc = [];
    all_global_std = [];
    count = 0;
    
    % Load data from each subject
    for i = 1:length(subjects)
        subject_id = subjects{i};
        subject_base_dir = fullfile(experiment_path, subject_id, session_id);
        subject_results_dir = fullfile(subject_base_dir, 'output', ...
            sprintf('entrainment_%s_%s', subject_id, session_id));
        subject_data_file = fullfile(subject_results_dir, sprintf('subject_%s_normalized_interp.mat', subject_id));
        
        % Load the global ISPC values for this subject
        global_file = fullfile(subject_results_dir, 'global_ispc.mat');
        global_ispc_value = NaN;
        global_std_value = NaN;
        
        if exist(global_file, 'file')
            try
                global_data = load(global_file);
                if isfield(global_data, 'globalISPC')
                    global_ispc_value = global_data.globalISPC;
                    fprintf('Loaded globalISPC: %.4f for subject %s\n', global_ispc_value, subject_id);
                end
                if isfield(global_data, 'globalStd')
                    global_std_value = global_data.globalStd;
                    fprintf('Loaded globalStd: %.4f for subject %s\n', global_std_value, subject_id);
                end
            catch err
                warning('Error loading global ISPC data for subject %s: %s', subject_id, err.message);
            end
        else
            fprintf('Global ISPC file not found for subject %s: %s\n', subject_id, global_file);
        end
        
        if exist(subject_data_file, 'file')
            fprintf('Loading interpolated normalized data for subject %s...\n', subject_id);
            
            try
                % Load subject's interpolated normalized data
                data = load(subject_data_file);
                
                % First subject initializes the data structures
                if count == 0
                    % Get segment types from the first subject
                    segment_types = fieldnames(data.interp_norm_topos);
                    change_types = fieldnames(data.interp_pct_topos);
                    
                    % Initialize data arrays
                    for s = 1:length(segment_types)
                        segment = segment_types{s};
                        all_norm_topos.(segment) = [];
                    end
                    
                    for d = 1:length(change_types)
                        change = change_types{d};
                        all_pct_changes.(change) = [];
                    end
                end
                
                % Store globalISPC and globalStd values for this subject
                all_global_ispc = [all_global_ispc; global_ispc_value];
                all_global_std = [all_global_std; global_std_value];
                
                % Add this subject's data to the accumulation
                for s = 1:length(segment_types)
                    segment = segment_types{s};
                    if isfield(data.interp_norm_topos, segment)
                        all_norm_topos.(segment) = [all_norm_topos.(segment), data.interp_norm_topos.(segment)];
                    end
                end
                
                for d = 1:length(change_types)
                    change = change_types{d};
                    if isfield(data.interp_pct_topos, change)
                        all_pct_changes.(change) = [all_pct_changes.(change), data.interp_pct_topos.(change)];
                    end
                end
                
                count = count + 1;
                
            catch err
                warning('Error loading data for subject %s: %s', subject_id, err.message);
            end
        else
            warning('Interpolated normalized data file not found for subject %s: %s', subject_id, subject_data_file);
        end
    end
    
    % If no data loaded, exit
    if count == 0
        warning('No valid interpolated normalized data found for any subjects in %s condition. Skipping group average.', condition);
        return;
    end
    
    fprintf('Successfully loaded interpolated normalized data from %d/%d subjects for %s condition\n', ...
        count, length(subjects), condition);
    
    % Calculate average of global ISPC and globalStd values
    mean_global_ispc = mean(all_global_ispc, 'omitnan');
    mean_global_std = mean(all_global_std, 'omitnan');
    
    fprintf('Mean global ISPC across subjects: %.4f, Mean global StdDev: %.4f\n', mean_global_ispc, mean_global_std);
    
    % Calculate group averages
    group_norm_topos = struct();
    group_pct_changes = struct();
    
    segment_types = fieldnames(all_norm_topos);
    for s = 1:length(segment_types)
        segment = segment_types{s};
        if ~isempty(all_norm_topos.(segment))
            group_norm_topos.(segment) = mean(all_norm_topos.(segment), 2, 'omitnan');
        end
    end
    
    change_types = fieldnames(all_pct_changes);
    for d = 1:length(change_types)
        change = change_types{d};
        if ~isempty(all_pct_changes.(change))
            group_pct_changes.(change) = mean(all_pct_changes.(change), 2, 'omitnan');
        end
    end
    
    % Calculate Z-score versions using the mean global standard deviation
    group_z_score_topos = struct();
    for s = 1:length(segment_types)
        segment = segment_types{s};
        if isfield(group_norm_topos, segment)
            % Use the actual formula for Z-scores: (x - mean) / std
            % For normalized data, mean is 1.0 by definition
            group_z_score_topos.(segment) = (group_norm_topos.(segment) - 1.0) / mean_global_std;
        end
    end
    
    % Calculate Z-score differences
    z_diff_topos = struct();
    
    if isfield(group_z_score_topos, 'pre_stim')
        if isfield(group_z_score_topos, 'early_stim')
            z_diff_topos.early_vs_pre = group_z_score_topos.early_stim - group_z_score_topos.pre_stim;
        end
        
        if isfield(group_z_score_topos, 'late_stim')
            z_diff_topos.late_vs_pre = group_z_score_topos.late_stim - group_z_score_topos.pre_stim;
        end
        
        if isfield(group_z_score_topos, 'post_stim')
            z_diff_topos.post_vs_pre = group_z_score_topos.post_stim - group_z_score_topos.pre_stim;
        end
    end
    
    if isfield(group_z_score_topos, 'early_stim') && isfield(group_z_score_topos, 'late_stim')
        z_diff_topos.late_vs_early = group_z_score_topos.late_stim - group_z_score_topos.early_stim;
    end
    
    % Save group average data
    group_data_file = fullfile(group_results_dir, sprintf('group_%s_normalized_average.mat', lower(condition)));
    save(group_data_file, 'group_norm_topos', 'group_pct_changes', 'group_z_score_topos', 'z_diff_topos', 'subjects', 'mean_global_ispc', 'mean_global_std', 'count');
    fprintf('Saved interpolated group average data to: %s\n', group_data_file);
    
    % Load electrode locations for visualization
    try
        % Read electrode locations
        chanlocs = readlocs(electrode_file);
        
        % Create EEG structure with channels
        EEG = struct('chanlocs', chanlocs, 'nbchan', length(chanlocs));
        
        % Create group-level topoplots
        createInterpolatedGroupTopoplots(EEG, group_norm_topos, group_pct_changes, group_z_score_topos, z_diff_topos, condition, group_results_dir, excluded_channels, count, mean_global_std);
    catch err
        warning('Failed to create topoplots: %s', err.message);
    end
end

function createInterpolatedGroupTopoplots(EEG, group_norm_topos, group_pct_changes, group_z_score_topos, z_diff_topos, condition, group_results_dir, excluded_channels, subject_count, mean_global_std)
    % Create topoplots for interpolated normalized group average data
    fprintf('Creating topoplots for %s group average...\n', condition);
    
    % Define standard segment types and their display names
    standard_segments = {'pre_stim', 'early_stim', 'late_stim', 'post_stim'};
    display_names = {'Pre-Stim', 'Early-Stim', 'Late-Stim', 'Post-Stim'};
    
    % Define color limits for normalized topoplots
    norm_limits = [0.5 1.5];
    pct_change_limits = [-10 10];
    
    % Create a list of channels to exclude
    excluded_indices = [];
    for i = 1:length(excluded_channels)
        for c = 1:length(EEG.chanlocs)
            if strcmpi(EEG.chanlocs(c).labels, excluded_channels{i})
                excluded_indices = [excluded_indices, c];
                break;
            end
        end
    end
    
    % Create a list of channels to include (all minus excluded)
    included_indices = setdiff(1:length(EEG.chanlocs), excluded_indices);
    
    fprintf('Found %d excluded channels, using %d channels for topoplots\n', ...
        length(excluded_indices), length(included_indices));
    
    %% PART 1: Create normalized topoplots
    % Create comprehensive figure with segments and percent changes
    figure('Name', sprintf('%s Interpolated Group Average', condition), 'Position', [50, 50, 1200, 800], 'visible', 'off');
    
    % Calculate total number of plots
    num_segments = sum(isfield(group_norm_topos, standard_segments));
    num_changes = length(fieldnames(group_pct_changes));
    total_plots = num_segments + num_changes;
    
    % Calculate a suitable grid layout
    grid_cols = ceil(sqrt(total_plots));
    grid_rows = ceil(total_plots / grid_cols);
    
    % Keep track of current subplot index
    plot_idx = 1;
    
    % First plot the normalized segments
    for i = 1:length(standard_segments)
        segment_type = standard_segments{i};
        if isfield(group_norm_topos, segment_type)
            subplot(grid_rows, grid_cols, plot_idx);
            
            % Prepare data for topoplot
            topo_data = group_norm_topos.(segment_type);
            
            try
                % Create topoplot with only valid channels
                topoplot(topo_data, EEG.chanlocs, 'maplimits', norm_limits, ...
                    'electrodes', 'on', 'plotchans', included_indices);
                
                title([display_names{i} ' (Normalized)']);
                colorbar;
            catch err
                fprintf('Error in topoplot for %s: %s\n', segment_type, err.message);
                % Create fallback visualization
                bar(topo_data(~isnan(topo_data)));
                title([display_names{i} ' (Normalized) - Simple View']);
            end
            
            plot_idx = plot_idx + 1;
        end
    end
    
    % Then plot the percent changes
    change_display_names = {
        'early_vs_pre', 'Early-Stim vs Pre-Stim';
        'late_vs_pre', 'Late-Stim vs Pre-Stim';
        'post_vs_pre', 'Post-Stim vs Pre-Stim';
        'late_vs_early', 'Late-Stim vs Early-Stim'
    };
    
    change_fields = fieldnames(group_pct_changes);
    for i = 1:length(change_fields)
        change_type = change_fields{i};
        
        subplot(grid_rows, grid_cols, plot_idx);
        
        % Find display name for this change
        display_idx = find(strcmp(change_display_names(:,1), change_type));
        if ~isempty(display_idx)
            display_name = change_display_names{display_idx, 2};
        else
            display_name = strrep(change_type, '_', ' ');
        end
        
        % Prepare data for topoplot
        topo_data = group_pct_changes.(change_type);
        
        try
            % Create topoplot with only valid channels
            topoplot(topo_data, EEG.chanlocs, 'maplimits', pct_change_limits, ...
                'electrodes', 'on', 'plotchans', included_indices);
            
            title(['% Change: ' display_name]);
            colorbar;
        catch err
            fprintf('Error in topoplot for change %s: %s\n', change_type, err.message);
            % Create fallback visualization
            bar(topo_data(~isnan(topo_data)));
            title(['% Change: ' display_name ' - Simple View']);
        end
        
        plot_idx = plot_idx + 1;
    end
    
    % Update the title to show correct subject count
    sgtitle(sprintf('%s Interpolated Group Average (n=%d)', condition, subject_count), 'FontSize', 16);
    
    % Save the comprehensive figure
    saveas(gcf, fullfile(group_results_dir, sprintf('group_%s_interpolated_topoplots.png', lower(condition))));
    fprintf('Saved interpolated group topoplots to: %s\n', fullfile(group_results_dir, sprintf('group_%s_interpolated_topoplots.png', lower(condition))));
    
    % Close the figure
    close(gcf);
    
    %% PART 2: Create Z-scored topoplots
    
    % Create a new figure for Z-score topoplots
    figure('Name', sprintf('%s Z-Scored Group Average', condition), 'Position', [50, 50, 1200, 800], 'visible', 'off');
    
    % Calculate total number of Z-score plots
    num_z_segments = sum(isfield(group_z_score_topos, standard_segments));
    num_z_diffs = length(fieldnames(z_diff_topos));
    total_z_plots = num_z_segments + num_z_diffs;
    
    % Calculate grid layout for Z-score plots
    z_grid_cols = ceil(sqrt(total_z_plots));
    z_grid_rows = ceil(total_z_plots / z_grid_cols);
    
    z_plot_idx = 1;
    
    % Define z-score color limits - now using +/- 2 standard deviations
    z_limits = [-2 2];
    
    % Plot Z-scored segments
    for i = 1:length(standard_segments)
        segment_type = standard_segments{i};
        if isfield(group_z_score_topos, segment_type)
            subplot(z_grid_rows, z_grid_cols, z_plot_idx);
            
            try
                % Create topoplot with only valid channels
                topoplot(group_z_score_topos.(segment_type), EEG.chanlocs, 'maplimits', z_limits, ...
                    'electrodes', 'on', 'plotchans', included_indices);
                
                title([display_names{i} ' (Z-Score)']);
                colorbar;
            catch err
                fprintf('Error in Z-score topoplot for %s: %s\n', segment_type, err.message);
                % Create fallback visualization
                bar(group_z_score_topos.(segment_type)(~isnan(group_z_score_topos.(segment_type))));
                title([display_names{i} ' (Z-Score) - Simple View']);
            end
            
            z_plot_idx = z_plot_idx + 1;
        end
    end
    
    % Plot Z-score differences
    z_diff_fields = fieldnames(z_diff_topos);
    for i = 1:length(z_diff_fields)
        z_diff_type = z_diff_fields{i};
        
        subplot(z_grid_rows, z_grid_cols, z_plot_idx);
        
        % Find display name for this difference
        display_idx = find(strcmp(change_display_names(:,1), z_diff_type));
        if ~isempty(display_idx)
            display_name = change_display_names{display_idx, 2};
        else
            display_name = strrep(z_diff_type, '_', ' ');
        end
        
        try
            % Create topoplot with only valid channels
            topoplot(z_diff_topos.(z_diff_type), EEG.chanlocs, 'maplimits', z_limits, ...
                'electrodes', 'on', 'plotchans', included_indices);
            
            title(['Z-Score Diff: ' display_name]);
            colorbar;
        catch err
            fprintf('Error in Z-score diff topoplot for %s: %s\n', z_diff_type, err.message);
            % Create fallback visualization
            bar(z_diff_topos.(z_diff_type)(~isnan(z_diff_topos.(z_diff_type))));
            title(['Z-Score Diff: ' display_name ' - Simple View']);
        end
        
        z_plot_idx = z_plot_idx + 1;
    end
    
    % Update the title to show correct subject count and mean global StdDev
    sgtitle(sprintf('%s Z-Scored Group Average (n=%d, Mean Global StdDev: %.4f)', ...
        condition, subject_count, mean_global_std), 'FontSize', 16);
    
    % Save the Z-score figure
    saveas(gcf, fullfile(group_results_dir, sprintf('group_%s_zscored_topoplots.png', lower(condition))));
    fprintf('Saved Z-scored group topoplots to: %s\n', fullfile(group_results_dir, sprintf('group_%s_zscored_topoplots.png', lower(condition))));
    
    % Close the figure
    close(gcf);
end
