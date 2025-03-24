function runProtocolGroupAnalysis(experiment_path, session_id, max_stim)
    % RUNPROTOCOLGROUPANALYSIS - Creates group averages for specific protocols
    %
    % This function creates group averages for ACTIVE and SHAM conditions
    % for each stimulation protocol (Stim_1, Stim_2, etc.) separately
    %
    % Inputs:
    %   experiment_path - Base path for the experiment data
    %   session_id - Session ID (e.g., 'N1')
    %   max_stim - Maximum number of stimulation protocols to process (default: 10)
    
    % Set default session if not provided
    if nargin < 2
        session_id = 'N1';
    end
    
    % Set default max stimulation number
    if nargin < 3
        max_stim = 10;
    end
    
    % Global variables for consistent color scaling
    global GLOBAL_SEGMENT_MIN GLOBAL_SEGMENT_MAX GLOBAL_DIFF_ABS_MAX;
    
    fprintf('Running protocol-specific group analysis for session %s (up to Stim_%d)...\n', ...
        session_id, max_stim);
    
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
    addpath(fullfile(src_path, 'eeglab'));
    
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
    
    % Define configuration for access to electrode file
    config = struct();
    config.assets_path = assets_path;
    config.utilities_path = utilities_path;
    
    % Get electrode file with prioritization
    config.electrode_file = fullfile(config.assets_path, 'GSN-HydroCel-256.sfp');
    if ~exist(config.electrode_file, 'file')
        config.electrode_file = fullfile(config.utilities_path, 'egi256_GSN_HydroCel.sfp');
        if ~exist(config.electrode_file, 'file')
            error('No electrode position file found. Cannot continue without electrode locations.');
        end
    end
    electrode_file = config.electrode_file;
    fprintf('Using electrode locations from: %s\n', electrode_file);
    
    % Define base group results directory
    base_group_dir = fullfile(experiment_path, 'group_results', session_id);
    if ~exist(base_group_dir, 'dir')
        [success, msg] = mkdir(base_group_dir);
        if ~success
            error('Failed to create group results directory %s: %s', base_group_dir, msg);
        end
        fprintf('Created group results directory: %s\n', base_group_dir);
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
    
    % Initialize color scale variables
    GLOBAL_SEGMENT_MIN = [];
    GLOBAL_SEGMENT_MAX = [];
    GLOBAL_DIFF_ABS_MAX = [];
    
    % Process each stimulation protocol
    for stim_num = 1:max_stim
        stim_id = sprintf('stim_%d', stim_num);
        fprintf('\n==== Processing %s ====\n', stim_id);
        
        % Create directories for this stimulation
        stim_group_dir = fullfile(base_group_dir, stim_id);
        if ~exist(stim_group_dir, 'dir')
            [success, msg] = mkdir(stim_group_dir);
            if ~success
                warning('Failed to create directory for %s: %s. Skipping.', stim_id, msg);
                continue;
            end
            fprintf('Created directory: %s\n', stim_group_dir);
        end
        
        % Initialize data structures for both conditions
        active_data = struct();
        sham_data = struct();
        
        % Load all data first to find global min/max for this protocol
        [active_data.group_topos, active_data.group_diffs, active_data.count] = loadProtocolGroupData('ACTIVE', active_subjects, session_id, experiment_path, electrode_file, excluded_channels, stim_num, stim_id);
        [sham_data.group_topos, sham_data.group_diffs, sham_data.count] = loadProtocolGroupData('SHAM', sham_subjects, session_id, experiment_path, electrode_file, excluded_channels, stim_num, stim_id);
        
        % Calculate protocol-specific min/max for consistent color scales between ACTIVE/SHAM
        segment_values = [];
        standard_segments = {'pre_stim', 'early_stim', 'late_stim', 'post_stim'};
        for i = 1:length(standard_segments)
            segment_type = standard_segments{i};
            if isfield(active_data.group_topos, segment_type)
                segment_values = [segment_values; active_data.group_topos.(segment_type)];
            end
            if isfield(sham_data.group_topos, segment_type)
                segment_values = [segment_values; sham_data.group_topos.(segment_type)];
            end
        end
        
        % Find min/max for segments
        protocol_segment_min = min(segment_values(~isnan(segment_values)));
        protocol_segment_max = max(segment_values(~isnan(segment_values)));
        
        % Find min/max for difference maps
        diff_values = [];
        diff_fields = fieldnames(active_data.group_diffs);
        for i = 1:length(diff_fields)
            if isfield(active_data.group_diffs, diff_fields{i})
                diff_values = [diff_values; active_data.group_diffs.(diff_fields{i})];
            end
            if isfield(sham_data.group_diffs, diff_fields{i})
                diff_values = [diff_values; sham_data.group_diffs.(diff_fields{i})];
            end
        end
        
        protocol_diff_abs_max = max(abs(diff_values(~isnan(diff_values))));
        
        fprintf('Protocol %s color scale: segments [%.4f, %.4f], diffs [%.4f, %.4f]\n', ...
            stim_id, protocol_segment_min, protocol_segment_max, -protocol_diff_abs_max, protocol_diff_abs_max);
            
        % Set color scales for this protocol only
        GLOBAL_SEGMENT_MIN = protocol_segment_min;
        GLOBAL_SEGMENT_MAX = protocol_segment_max;
        GLOBAL_DIFF_ABS_MAX = protocol_diff_abs_max;
        
        % Now create visualizations using consistent color scaling
        if active_data.count > 0
            createProtocolVisualizations('ACTIVE', active_data, stim_group_dir, electrode_file, excluded_channels, stim_id);
        end
        
        if sham_data.count > 0
            createProtocolVisualizations('SHAM', sham_data, stim_group_dir, electrode_file, excluded_channels, stim_id);
        end
    end
    
    fprintf('\nProtocol-specific group analysis completed successfully.\n');
end

function [group_topos, group_diffs, count] = loadProtocolGroupData(condition, subjects, session_id, experiment_path, electrode_file, excluded_channels, stim_num, stim_id)
    % Load and process protocol-specific data for a condition
    fprintf('Loading data for %s condition, %s (%d subjects)...\n', condition, stim_id, length(subjects));
    
    % Load standard electrode locations for interpolation
    chanlocs = readlocs(electrode_file);
    template_EEG = struct();
    template_EEG.chanlocs = chanlocs;
    template_EEG.nbchan = length(chanlocs);
    
    % Initialize data structures
    all_topos = struct();
    all_diffs = struct();
    count = 0;
    
    % Load data from each subject
    for i = 1:length(subjects)
        subject_id = subjects{i};
        subject_dir = fullfile(experiment_path, subject_id, session_id, 'output', ...
            sprintf('entrainment_%s_%s', subject_id, session_id), stim_id);
        subject_data_file = fullfile(subject_dir, 'results.mat');
        
        if exist(subject_data_file, 'file')
            fprintf('Loading data for subject %s, %s...\n', subject_id, stim_id);
            
            try
                % Load subject's stimulation data
                data = load(subject_data_file);
                
                % Check if the expected fields exist
                if ~isfield(data, 'ispc_results') || ~isfield(data, 'stim_segments')
                    warning('Subject %s: %s data file does not contain expected fields. Skipping.', ...
                        subject_id, stim_id);
                    continue;
                end
                
                % Get ISPC results and segments
                ispc_results = data.ispc_results;
                stim_segments = data.stim_segments;
                
                % Standard segment names (normalize names by removing any numeric suffixes)
                standard_segments = {'pre_stim', 'early_stim', 'late_stim', 'post_stim'};
                
                % Create 'final_topos' structure similar to the subject average files
                subject_topos = struct();
                segment_names = fieldnames(stim_segments);
                
                for s = 1:length(standard_segments)
                    segment_type = standard_segments{s};
                    
                    % Find the matching segment (may have a suffix)
                    matching_segment = '';
                    for seg_idx = 1:length(segment_names)
                        if contains(segment_names{seg_idx}, segment_type)
                            matching_segment = segment_names{seg_idx};
                            break;
                        end
                    end
                    
                    if ~isempty(matching_segment)
                        % Extract segment range
                        segment_range = stim_segments.(matching_segment);
                        
                        % Calculate mean ISPC over time for this segment
                        subject_topos.(segment_type) = mean(ispc_results(:, segment_range(1):segment_range(2)), 2);
                    end
                end
                
                % Calculate differences for this subject
                subject_diffs = struct();
                
                if isfield(subject_topos, 'pre_stim')
                    if isfield(subject_topos, 'early_stim')
                        subject_diffs.early_minus_pre = subject_topos.early_stim - subject_topos.pre_stim;
                    end
                    
                    if isfield(subject_topos, 'late_stim')
                        subject_diffs.late_minus_pre = subject_topos.late_stim - subject_topos.pre_stim;
                    end
                    
                    if isfield(subject_topos, 'post_stim')
                        subject_diffs.post_minus_pre = subject_topos.post_stim - subject_topos.pre_stim;
                    end
                end
                
                if isfield(subject_topos, 'early_stim') && isfield(subject_topos, 'late_stim')
                    subject_diffs.late_minus_early = subject_topos.late_stim - subject_topos.early_stim;
                end
                
                % Interpolate missing channels
                interp_topos = interpolateData(subject_topos, template_EEG, excluded_channels);
                interp_diffs = interpolateData(subject_diffs, template_EEG, excluded_channels);
                
                % First subject initializes the data structures
                if count == 0
                    % Get segment types 
                    segment_types = fieldnames(interp_topos);
                    diff_types = fieldnames(interp_diffs);
                    
                    % Initialize data arrays
                    for s = 1:length(segment_types)
                        segment = segment_types{s};
                        all_topos.(segment) = [];
                    end
                    
                    for d = 1:length(diff_types)
                        diff = diff_types{d};
                        all_diffs.(diff) = [];
                    end
                end
                
                % Add this subject's data to the accumulation
                for s = 1:length(segment_types)
                    segment = segment_types{s};
                    if isfield(interp_topos, segment)
                        all_topos.(segment) = [all_topos.(segment), interp_topos.(segment)];
                    end
                end
                
                for d = 1:length(diff_types)
                    diff = diff_types{d};
                    if isfield(interp_diffs, diff)
                        all_diffs.(diff) = [all_diffs.(diff), interp_diffs.(diff)];
                    end
                end
                
                count = count + 1;
                
            catch err
                warning('Error processing subject %s, %s: %s', subject_id, stim_id, err.message);
            end
        else
            % For diagnostics, let's check if the parent directory exists
            subject_results_dir = fullfile(experiment_path, subject_id, session_id, 'output', ...
                sprintf('entrainment_%s_%s', subject_id, session_id));
            
            if exist(subject_results_dir, 'dir')
                % Check if any stim directories exist
                stim_dirs = dir(fullfile(subject_results_dir, 'stim_*'));
                if ~isempty(stim_dirs)
                    fprintf('Subject %s has stimulation data but not for %s\n', subject_id, stim_id);
                else
                    fprintf('Subject %s does not have any stimulation data\n', subject_id);
                end
            else
                fprintf('Subject %s does not have a results directory\n', subject_id);
            end
        end
    end
    
    % If no data loaded, exit
    if count == 0
        warning('No valid data found for any subjects in %s condition for %s. Skipping group average.', ...
            condition, stim_id);
        group_topos = struct();
        group_diffs = struct();
        return;
    end
    
    fprintf('Successfully loaded data from %d/%d subjects for %s condition, %s\n', ...
        count, length(subjects), condition, stim_id);
    
    % Calculate group averages
    group_topos = struct();
    group_diffs = struct();
    
    segment_types = fieldnames(all_topos);
    for s = 1:length(segment_types)
        segment = segment_types{s};
        if ~isempty(all_topos.(segment))
            group_topos.(segment) = mean(all_topos.(segment), 2);
        end
    end
    
    diff_types = fieldnames(all_diffs);
    for d = 1:length(diff_types)
        diff = diff_types{d};
        if ~isempty(all_diffs.(diff))
            group_diffs.(diff) = mean(all_diffs.(diff), 2);
        end
    end
end

function interp_data = interpolateData(data_struct, template_EEG, excluded_channels)
    % Interpolate data structure to match template
    interp_data = struct();
    fields = fieldnames(data_struct);
    
    for f = 1:length(fields)
        field_name = fields{f};
        original_data = data_struct.(field_name);
        n_channels = length(original_data);
        
        % Get indices of excluded channels in the template
        excluded_indices = [];
        for e = 1:length(excluded_channels)
            for c = 1:length(template_EEG.chanlocs)
                if strcmpi(template_EEG.chanlocs(c).labels, excluded_channels{e})
                    excluded_indices = [excluded_indices, c];
                    break;
                end
            end
        end
        
        % Create a full data array with NaNs
        full_data = nan(length(template_EEG.chanlocs), 1);
        
        % Fill in the data we have (assuming it starts from channel 1)
        if n_channels <= length(template_EEG.chanlocs)
            full_data(1:n_channels) = original_data;
        else
            % If we have more channels than template, truncate
            full_data = original_data(1:length(template_EEG.chanlocs));
        end
        
        % Interpolate NaN values using nearest neighbors
        valid_indices = find(~isnan(full_data));
        nan_indices = find(isnan(full_data));
        
        % Remove excluded channels from interpolation targets
        nan_indices = setdiff(nan_indices, excluded_indices);
        
        if ~isempty(valid_indices) && ~isempty(nan_indices)
            % Get coordinates for valid and NaN points
            valid_coords = zeros(length(valid_indices), 3);
            nan_coords = zeros(length(nan_indices), 3);
            
            for i = 1:length(valid_indices)
                idx = valid_indices(i);
                valid_coords(i, 1) = template_EEG.chanlocs(idx).X;
                valid_coords(i, 2) = template_EEG.chanlocs(idx).Y;
                valid_coords(i, 3) = template_EEG.chanlocs(idx).Z;
            end
            
            for i = 1:length(nan_indices)
                idx = nan_indices(i);
                nan_coords(i, 1) = template_EEG.chanlocs(idx).X;
                nan_coords(i, 2) = template_EEG.chanlocs(idx).Y;
                nan_coords(i, 3) = template_EEG.chanlocs(idx).Z;
            end
            
            % Perform interpolation (simplified inverse distance weighting)
            for i = 1:length(nan_indices)
                idx = nan_indices(i);
                
                % Calculate distances to all valid channels
                distances = zeros(length(valid_indices), 1);
                for j = 1:length(valid_indices)
                    valid_idx = valid_indices(j);
                    
                    % Calculate Euclidean distance
                    distances(j) = sqrt(sum((nan_coords(i, :) - valid_coords(j, :)).^2));
                end
                
                % Find the 5 nearest channels (or fewer if less available)
                [sorted_dist, nearest_idx] = sort(distances);
                nearest_idx = nearest_idx(1:min(5, length(nearest_idx)));
                sorted_dist = sorted_dist(1:min(5, length(sorted_dist)));
                
                % Inverse distance weighting
                weights = 1 ./ sorted_dist;
                weights = weights / sum(weights);
                
                % Calculate weighted average
                full_data(idx) = sum(weights .* full_data(valid_indices(nearest_idx)));
            end
        end
        
        % Make sure excluded channels are set to NaN
        full_data(excluded_indices) = NaN;
        
        % Store interpolated data
        interp_data.(field_name) = full_data;
    end
end

function createProtocolVisualizations(condition, data, group_results_dir, electrode_file, excluded_channels, stim_id)
    % Creates visualizations for protocol-specific data with consistent color scales
    fprintf('Creating visualizations for %s condition, %s...\n', condition, stim_id);
    
    global GLOBAL_SEGMENT_MIN GLOBAL_SEGMENT_MAX GLOBAL_DIFF_ABS_MAX;
    
    % Save group average data
    group_data_file = fullfile(group_results_dir, sprintf('group_%s_average.mat', lower(condition)));
    group_topos = data.group_topos;
    group_diffs = data.group_diffs;
    count = data.count;
    save(group_data_file, 'group_topos', 'group_diffs', 'count');
    fprintf('Saved group average data to: %s\n', group_data_file);
    
    % Create topoplots
    try
        EEG = struct();
        EEG.chanlocs = readlocs(electrode_file);
        EEG.nbchan = length(EEG.chanlocs);
        createProtocolTopoplots(EEG, group_topos, group_diffs, condition, group_results_dir, excluded_channels, count, stim_id);
    catch err
        warning('Failed to create topoplots: %s', err.message);
    end
end

function createProtocolTopoplots(EEG, group_topos, group_diffs, condition, group_results_dir, excluded_channels, subject_count, stim_id)
    % Create topoplots for protocol-specific group average data
    fprintf('Creating topoplots for %s condition group average (%s)...\n', condition, stim_id);
    
    global GLOBAL_SEGMENT_MIN GLOBAL_SEGMENT_MAX GLOBAL_DIFF_ABS_MAX;
    
    % Define standard segment types and their display names
    standard_segments = {'pre_stim', 'early_stim', 'late_stim', 'post_stim'};
    display_names = {'Pre-Stim', 'Early-Stim', 'Late-Stim', 'Post-Stim'};
    
    % Use global color limits for consistent scaling across conditions
    segment_min = GLOBAL_SEGMENT_MIN;
    segment_max = GLOBAL_SEGMENT_MAX;
    diff_limit = [-GLOBAL_DIFF_ABS_MAX GLOBAL_DIFF_ABS_MAX]; % Symmetric colormap for differences
    
    % Create comprehensive figure with segments and differences
    figure('Name', sprintf('%s Group Average - %s', condition, stim_id), 'Position', [50, 50, 1200, 800], 'visible', 'off');
    
    % Calculate total number of plots
    num_segments = sum(isfield(group_topos, standard_segments));
    num_diffs = length(fieldnames(group_diffs));
    total_plots = num_segments + num_diffs;
    
    % Calculate a suitable grid layout
    grid_cols = ceil(sqrt(total_plots));
    grid_rows = ceil(total_plots / grid_cols);
    
    % Keep track of current subplot index
    plot_idx = 1;
    
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
    
    % First plot the segments
    for i = 1:length(standard_segments)
        segment_type = standard_segments{i};
        if isfield(group_topos, segment_type)
            subplot(grid_rows, grid_cols, plot_idx);
            
            % Prepare data for topoplot
            topo_data = group_topos.(segment_type);
            
            try
                % Try to use standard EEGLAB topoplot first
                topoplot(topo_data, EEG.chanlocs, 'maplimits', [segment_min segment_max], ...
                    'electrodes', 'labels', 'efontsize', 6, 'plotchans', included_indices);
                
                title(display_names{i});
                colorbar;
            catch err
                fprintf('Error in topoplot for %s: %s\n', segment_type, err.message);
                
                % Fallback to simple visualization - just plot the values
                bar(topo_data(~isnan(topo_data)));
                title([display_names{i} ' (Simple View)']);
                ylabel('Value');
                xlabel('Channel Index');
            end
            
            plot_idx = plot_idx + 1;
        end
    end
    
    % Then plot the differences
    diff_fields = fieldnames(group_diffs);
    diff_display_names = {
        'early_minus_pre', 'Early-Stim - Pre-Stim';
        'late_minus_pre', 'Late-Stim - Pre-Stim';
        'post_minus_pre', 'Post-Stim - Pre-Stim';
        'late_minus_early', 'Late-Stim - Early-Stim'
    };
    
    for i = 1:length(diff_fields)
        diff_type = diff_fields{i};
        
        subplot(grid_rows, grid_cols, plot_idx);
        
        % Find display name for this difference
        display_idx = find(strcmp(diff_display_names(:,1), diff_type));
        if ~isempty(display_idx)
            display_name = diff_display_names{display_idx, 2};
        else
            display_name = strrep(diff_type, '_', ' ');
        end
        
        % Prepare data for topoplot
        topo_data = group_diffs.(diff_type);
        
        try
            % Try to use standard EEGLAB topoplot first
            topoplot(topo_data, EEG.chanlocs, 'maplimits', diff_limit, ...
                'electrodes', 'labels', 'efontsize', 6, 'plotchans', included_indices);
            
            title(display_name);
            colorbar;
        catch err
            fprintf('Error in topoplot for diff %s: %s\n', diff_type, err.message);
            
            % Fallback to simple visualization - just plot the values
            bar(topo_data(~isnan(topo_data)));
            title([display_name ' (Simple View)']);
            ylabel('Value');
            xlabel('Channel Index');
        end
        
        plot_idx = plot_idx + 1;
    end
    
    % Update the title to show correct subject count and protocol
    sgtitle(sprintf('%s Group Average - %s (n=%d)', condition, stim_id, subject_count), 'FontSize', 16);
    
    % Save the comprehensive figure
    saveas(gcf, fullfile(group_results_dir, sprintf('group_%s_topoplots.png', lower(condition))));
    fprintf('Saved group topoplots to: %s\n', fullfile(group_results_dir, sprintf('group_%s_topoplots.png', lower(condition))));
    
    % Create a list of visible channels and save it
    visible_channels = {EEG.chanlocs(included_indices).labels};
    channel_list_file = fullfile(group_results_dir, sprintf('group_%s_visible_channels.txt', lower(condition)));
    fid = fopen(channel_list_file, 'w');
    if fid ~= -1
        fprintf(fid, 'Visible channels for %s group analysis - %s (%d channels):\n\n', condition, stim_id, length(visible_channels));
        for i = 1:length(visible_channels)
            fprintf(fid, '%s', visible_channels{i});
            if mod(i, 10) == 0
                fprintf(fid, '\n');
            else
                fprintf(fid, ', ');
            end
        end
        fclose(fid);
        fprintf('Saved visible channel list to: %s\n', channel_list_file);
    end
    
    % Close the figure
    close(gcf);
end