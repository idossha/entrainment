function compareActiveAndShamTopoplots(experiment_path, session_id)
    % COMPAREACTIVESHAMTOPOPLOTS - Compares ACTIVE and SHAM group topoplots
    %
    % This function calculates and visualizes the difference between ACTIVE and
    % SHAM group topoplots for both normalized and Z-scored data
    %
    % Inputs:
    %   experiment_path - Base path for the experiment data
    %   session_id - Session ID (e.g., 'N1')
    
    % Set default session if not provided
    if nargin < 2
        session_id = 'N1';
    end
    
    fprintf('Comparing ACTIVE vs SHAM group topoplots for session %s...\n', session_id);
    
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
    
    % Define group results directory
    group_results_dir = fullfile(experiment_path, 'group_results', session_id, 'normalized_interp');
    if ~exist(group_results_dir, 'dir')
        error('Group results directory not found: %s', group_results_dir);
    end
    
    % Define electrode file path
    electrode_file = fullfile(code_path, 'src/utils/egi256_GSN_HydroCel.sfp');
    fprintf('Using electrode locations from: %s\n', electrode_file);
    
    % Define comparison results directory
    comparison_dir = fullfile(experiment_path, 'group_results', session_id, 'active_vs_sham');
    if ~exist(comparison_dir, 'dir')
        [success, msg] = mkdir(comparison_dir);
        if ~success
            error('Failed to create comparison directory %s: %s', comparison_dir, msg);
        end
        fprintf('Created comparison directory: %s\n', comparison_dir);
    end
    
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
    
    % Load ACTIVE group data
    active_file = fullfile(group_results_dir, 'group_active_normalized_average.mat');
    if ~exist(active_file, 'file')
        error('ACTIVE group data file not found: %s', active_file);
    end
    
    fprintf('Loading ACTIVE group data...\n');
    active_data = load(active_file);
    
    % Load SHAM group data
    sham_file = fullfile(group_results_dir, 'group_sham_normalized_average.mat');
    if ~exist(sham_file, 'file')
        error('SHAM group data file not found: %s', sham_file);
    end
    
    fprintf('Loading SHAM group data...\n');
    sham_data = load(sham_file);
    
    % Calculate the subject counts for each group
    active_count = active_data.count;
    sham_count = sham_data.count;
    
    % Calculate the difference between normalized topoplots (ACTIVE - SHAM)
    norm_diff_topos = struct();
    norm_segments = fieldnames(active_data.group_norm_topos);
    
    for i = 1:length(norm_segments)
        segment = norm_segments{i};
        if isfield(sham_data.group_norm_topos, segment)
            norm_diff_topos.(segment) = active_data.group_norm_topos.(segment) - sham_data.group_norm_topos.(segment);
        end
    end
    
    % Calculate the difference between percent change topoplots (ACTIVE - SHAM)
    pct_change_diff = struct();
    pct_change_fields = fieldnames(active_data.group_pct_changes);
    
    for i = 1:length(pct_change_fields)
        change = pct_change_fields{i};
        if isfield(sham_data.group_pct_changes, change)
            pct_change_diff.(change) = active_data.group_pct_changes.(change) - sham_data.group_pct_changes.(change);
        end
    end
    
    % Calculate the difference between Z-scored topoplots (ACTIVE - SHAM)
    z_score_diff_topos = struct();
    z_segments = fieldnames(active_data.group_z_score_topos);
    
    for i = 1:length(z_segments)
        segment = z_segments{i};
        if isfield(sham_data.group_z_score_topos, segment)
            z_score_diff_topos.(segment) = active_data.group_z_score_topos.(segment) - sham_data.group_z_score_topos.(segment);
        end
    end
    
    % Calculate the difference between Z-score differences (ACTIVE - SHAM)
    z_diff_diff_topos = struct();
    z_diff_fields = fieldnames(active_data.z_diff_topos);
    
    for i = 1:length(z_diff_fields)
        diff_type = z_diff_fields{i};
        if isfield(sham_data.z_diff_topos, diff_type)
            z_diff_diff_topos.(diff_type) = active_data.z_diff_topos.(diff_type) - sham_data.z_diff_topos.(diff_type);
        end
    end
    
    % Calculate the mean global standard deviation (pooled from both groups)
    active_global_std = active_data.mean_global_std;
    sham_global_std = sham_data.mean_global_std;
    pooled_global_std = sqrt((active_global_std^2 * active_count + sham_global_std^2 * sham_count) / (active_count + sham_count));
    
    fprintf('ACTIVE Global StdDev: %.4f (n=%d)\n', active_global_std, active_count);
    fprintf('SHAM Global StdDev: %.4f (n=%d)\n', sham_global_std, sham_count);
    fprintf('Pooled Global StdDev: %.4f (n=%d)\n', pooled_global_std, active_count + sham_count);
    
    % Save the difference data
    diff_file = fullfile(comparison_dir, 'active_vs_sham_diffs.mat');
    save(diff_file, 'norm_diff_topos', 'pct_change_diff', 'z_score_diff_topos', 'z_diff_diff_topos', ...
        'active_count', 'sham_count', 'active_global_std', 'sham_global_std', 'pooled_global_std');
    fprintf('Saved difference data to: %s\n', diff_file);
    
    % Create visualizations
    try
        % Read electrode locations
        chanlocs = readlocs(electrode_file);
        
        % Create EEG structure with channels
        EEG = struct('chanlocs', chanlocs, 'nbchan', length(chanlocs));
        
        % Create comparison topoplots
        createComparisonTopoplots(EEG, norm_diff_topos, pct_change_diff, z_score_diff_topos, z_diff_diff_topos, ...
            comparison_dir, excluded_channels, active_count, sham_count, pooled_global_std);
    catch err
        warning('Failed to create comparison topoplots: %s', err.message);
    end
    
    fprintf('Comparison of ACTIVE vs SHAM topoplots completed.\n');
end

function createComparisonTopoplots(EEG, norm_diff_topos, pct_change_diff, z_score_diff_topos, z_diff_diff_topos, ...
    comparison_dir, excluded_channels, active_count, sham_count, pooled_global_std)
    % Create topoplots to visualize the differences between ACTIVE and SHAM groups
    
    % Define standard segment types and their display names
    standard_segments = {'pre_stim', 'early_stim', 'late_stim', 'post_stim'};
    display_names = {'Pre-Stim', 'Early-Stim', 'Late-Stim', 'Post-Stim'};
    
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
    
    %% PART 1: Create normalized difference topoplots
    
    % Define color limits for different plot types
    norm_limits = [-0.5 0.5];     % Normalized differences
    pct_change_limits = [-20 20]; % Percent change differences
    z_score_limits = [-2 2];      % Z-score differences
    
    % Create figure for normalized differences (ACTIVE - SHAM)
    figure('Name', 'ACTIVE vs SHAM Normalized Differences', 'Position', [50, 50, 1200, 600], 'visible', 'off');
    
    % Calculate total number of plots
    num_segments = sum(isfield(norm_diff_topos, standard_segments));
    
    % Calculate a suitable grid layout
    grid_cols = ceil(sqrt(num_segments));
    grid_rows = ceil(num_segments / grid_cols);
    
    % Plot normalized segment differences
    for i = 1:length(standard_segments)
        segment_type = standard_segments{i};
        if isfield(norm_diff_topos, segment_type)
            subplot(grid_rows, grid_cols, i);
            
            try
                % Create topoplot with only valid channels
                topoplot(norm_diff_topos.(segment_type), EEG.chanlocs, 'maplimits', norm_limits, ...
                    'electrodes', 'on', 'plotchans', included_indices);
                
                title([display_names{i} ' (ACTIVE - SHAM)']);
                colorbar;
            catch err
                fprintf('Error in difference topoplot for %s: %s\n', segment_type, err.message);
                % Create fallback visualization
                bar(norm_diff_topos.(segment_type)(~isnan(norm_diff_topos.(segment_type))));
                title([display_names{i} ' (ACTIVE - SHAM) - Simple View']);
            end
        end
    end
    
    % Update the title to show subject counts
    sgtitle(sprintf('Normalized Difference: ACTIVE (n=%d) - SHAM (n=%d)', active_count, sham_count), 'FontSize', 16);
    
    % Save the figure
    saveas(gcf, fullfile(comparison_dir, 'active_vs_sham_normalized_diff.png'));
    fprintf('Saved normalized difference topoplots to: %s\n', fullfile(comparison_dir, 'active_vs_sham_normalized_diff.png'));
    
    % Close the figure
    close(gcf);
    
    %% PART 2: Create percent change difference topoplots
    
    % Create figure for percent change differences (ACTIVE - SHAM)
    figure('Name', 'ACTIVE vs SHAM Percent Change Differences', 'Position', [50, 50, 1200, 600], 'visible', 'off');
    
    % Calculate total number of plots
    num_changes = length(fieldnames(pct_change_diff));
    
    % Calculate a suitable grid layout
    grid_cols = ceil(sqrt(num_changes));
    grid_rows = ceil(num_changes / grid_cols);
    
    % Define change display names for better labeling
    change_display_names = {
        'early_vs_pre', 'Early-Stim vs Pre-Stim';
        'late_vs_pre', 'Late-Stim vs Pre-Stim';
        'post_vs_pre', 'Post-Stim vs Pre-Stim';
        'late_vs_early', 'Late-Stim vs Early-Stim'
    };
    
    % Plot percent change differences
    change_fields = fieldnames(pct_change_diff);
    for i = 1:length(change_fields)
        change_type = change_fields{i};
        subplot(grid_rows, grid_cols, i);
        
        % Find display name for this change
        display_idx = find(strcmp(change_display_names(:,1), change_type));
        if ~isempty(display_idx)
            display_name = change_display_names{display_idx, 2};
        else
            display_name = strrep(change_type, '_', ' ');
        end
        
        try
            % Create topoplot with only valid channels
            topoplot(pct_change_diff.(change_type), EEG.chanlocs, 'maplimits', pct_change_limits, ...
                'electrodes', 'on', 'plotchans', included_indices);
            
            title(['% Change Diff: ' display_name]);
            colorbar;
        catch err
            fprintf('Error in percent change diff topoplot for %s: %s\n', change_type, err.message);
            % Create fallback visualization
            bar(pct_change_diff.(change_type)(~isnan(pct_change_diff.(change_type))));
            title(['% Change Diff: ' display_name ' - Simple View']);
        end
    end
    
    % Update the title to show subject counts
    sgtitle(sprintf('Percent Change Difference: ACTIVE (n=%d) - SHAM (n=%d)', active_count, sham_count), 'FontSize', 16);
    
    % Save the figure
    saveas(gcf, fullfile(comparison_dir, 'active_vs_sham_pct_change_diff.png'));
    fprintf('Saved percent change difference topoplots to: %s\n', fullfile(comparison_dir, 'active_vs_sham_pct_change_diff.png'));
    
    % Close the figure
    close(gcf);
    
    %% PART 3: Create Z-scored difference topoplots
    
    % Create figure for Z-scored differences (ACTIVE - SHAM)
    figure('Name', 'ACTIVE vs SHAM Z-Scored Differences', 'Position', [50, 50, 1200, 600], 'visible', 'off');
    
    % Calculate total number of plots
    num_z_segments = sum(isfield(z_score_diff_topos, standard_segments));
    
    % Calculate a suitable grid layout
    z_grid_cols = ceil(sqrt(num_z_segments));
    z_grid_rows = ceil(num_z_segments / z_grid_cols);
    
    % Plot Z-scored segment differences
    for i = 1:length(standard_segments)
        segment_type = standard_segments{i};
        if isfield(z_score_diff_topos, segment_type)
            subplot(z_grid_rows, z_grid_cols, i);
            
            try
                % Create topoplot with only valid channels
                topoplot(z_score_diff_topos.(segment_type), EEG.chanlocs, 'maplimits', z_score_limits, ...
                    'electrodes', 'on', 'plotchans', included_indices);
                
                title([display_names{i} ' (ACTIVE - SHAM)']);
                colorbar;
            catch err
                fprintf('Error in Z-score diff topoplot for %s: %s\n', segment_type, err.message);
                % Create fallback visualization
                bar(z_score_diff_topos.(segment_type)(~isnan(z_score_diff_topos.(segment_type))));
                title([display_names{i} ' (ACTIVE - SHAM) - Simple View']);
            end
        end
    end
    
    % Update the title to show subject counts and pooled standard deviation
    sgtitle(sprintf('Z-Scored Difference: ACTIVE (n=%d) - SHAM (n=%d), Pooled StdDev: %.4f', ...
        active_count, sham_count, pooled_global_std), 'FontSize', 16);
    
    % Save the figure
    saveas(gcf, fullfile(comparison_dir, 'active_vs_sham_zscored_diff.png'));
    fprintf('Saved Z-scored difference topoplots to: %s\n', fullfile(comparison_dir, 'active_vs_sham_zscored_diff.png'));
    
    % Close the figure
    close(gcf);
    
    %% PART 4: Create Z-score difference of differences topoplots
    
    % Create figure for Z-score differences of differences (ACTIVE - SHAM)
    figure('Name', 'ACTIVE vs SHAM Z-Score Difference of Differences', 'Position', [50, 50, 1200, 600], 'visible', 'off');
    
    % Calculate total number of plots
    num_z_diffs = length(fieldnames(z_diff_diff_topos));
    
    % Calculate a suitable grid layout
    z_diff_grid_cols = ceil(sqrt(num_z_diffs));
    z_diff_grid_rows = ceil(num_z_diffs / z_diff_grid_cols);
    
    % Plot Z-score difference of differences
    z_diff_fields = fieldnames(z_diff_diff_topos);
    for i = 1:length(z_diff_fields)
        z_diff_type = z_diff_fields{i};
        subplot(z_diff_grid_rows, z_diff_grid_cols, i);
        
        % Find display name for this difference
        display_idx = find(strcmp(change_display_names(:,1), z_diff_type));
        if ~isempty(display_idx)
            display_name = change_display_names{display_idx, 2};
        else
            display_name = strrep(z_diff_type, '_', ' ');
        end
        
        try
            % Create topoplot with only valid channels
            topoplot(z_diff_diff_topos.(z_diff_type), EEG.chanlocs, 'maplimits', z_score_limits, ...
                'electrodes', 'on', 'plotchans', included_indices);
            
            title(['Z-Score Diff: ' display_name]);
            colorbar;
        catch err
            fprintf('Error in Z-score diff of diffs topoplot for %s: %s\n', z_diff_type, err.message);
            % Create fallback visualization
            bar(z_diff_diff_topos.(z_diff_type)(~isnan(z_diff_diff_topos.(z_diff_type))));
            title(['Z-Score Diff: ' display_name ' - Simple View']);
        end
    end
    
    % Update the title to show subject counts
    sgtitle(sprintf('Z-Score Difference of Differences: ACTIVE (n=%d) - SHAM (n=%d)', ...
        active_count, sham_count), 'FontSize', 16);
    
    % Save the figure
    saveas(gcf, fullfile(comparison_dir, 'active_vs_sham_zscore_diff_of_diffs.png'));
    fprintf('Saved Z-score difference of differences topoplots to: %s\n', ...
        fullfile(comparison_dir, 'active_vs_sham_zscore_diff_of_diffs.png'));
    
    % Close the figure
    close(gcf);
end
