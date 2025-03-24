function runGroupAnalysisWithInterp(experiment_path, session_id)
    % RUNGROUPANALYSISWITHINTERP - Performs group-level analysis across subjects
    % using interpolated data
    %
    % This function creates group averages for ACTIVE and SHAM conditions
    % using the interpolated subject data files
    %
    % Inputs:
    %   experiment_path - Base path for the experiment data
    %   session_id - Session ID (e.g., 'N1')
    
    % Set default session if not provided
    if nargin < 2
        session_id = 'N1';
    end
    
    % Global variables to store color limits for consistent scaling
    global GLOBAL_SEGMENT_MIN GLOBAL_SEGMENT_MAX GLOBAL_DIFF_ABS_MAX;
    
    fprintf('Running group analysis with interpolated data for session %s...\n', session_id);
    
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
    
    % Define group results directory
    group_results_dir = fullfile(experiment_path, 'group_results', session_id);
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
    
    % Initialize data structures to store all group data
    active_data = struct();
    sham_data = struct();
    
    % First pass: load and process data for both conditions to find global min/max
    [active_data.group_topos, active_data.group_diffs, active_data.count] = loadGroupData('ACTIVE', active_subjects, session_id, experiment_path, electrode_file, excluded_channels);
    [sham_data.group_topos, sham_data.group_diffs, sham_data.count] = loadGroupData('SHAM', sham_subjects, session_id, experiment_path, electrode_file, excluded_channels);
    
    % Calculate min/max for consistent color scales across ACTIVE/SHAM
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
    GLOBAL_SEGMENT_MIN = min(segment_values(~isnan(segment_values)));
    GLOBAL_SEGMENT_MAX = max(segment_values(~isnan(segment_values)));
    
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
    
    GLOBAL_DIFF_ABS_MAX = max(abs(diff_values(~isnan(diff_values))));
    
    fprintf('Global analysis color scale: segments [%.4f, %.4f], diffs [%.4f, %.4f]\n', ...
        GLOBAL_SEGMENT_MIN, GLOBAL_SEGMENT_MAX, -GLOBAL_DIFF_ABS_MAX, GLOBAL_DIFF_ABS_MAX);
    
    % Create group averages for each condition using consistent color scales
    createGroupVisualizations('ACTIVE', active_data, group_results_dir, electrode_file, excluded_channels);
    createGroupVisualizations('SHAM', sham_data, group_results_dir, electrode_file, excluded_channels);
    
    fprintf('Group analysis completed successfully.\n');
end

function [group_topos, group_diffs, count] = loadGroupData(condition, subjects, session_id, experiment_path, electrode_file, excluded_channels)
    % Loads and processes data for a condition without creating visualizations
    fprintf('Loading data for %s condition (%d subjects)...\n', condition, length(subjects));
    
    % Initialize data structures
    all_topos = struct();
    all_diffs = struct();
    count = 0;
    
    % Load data from each subject
    for i = 1:length(subjects)
        subject_id = subjects{i};
        subject_results_dir = fullfile(experiment_path, subject_id, session_id, 'output', ...
            sprintf('entrainment_%s_%s', subject_id, session_id));
        subject_data_file = fullfile(subject_results_dir, sprintf('subject_%s_interp.mat', subject_id));
        
        if exist(subject_data_file, 'file')
            fprintf('Loading interpolated data for subject %s...\n', subject_id);
            
            try
                % Load subject's interpolated data
                data = load(subject_data_file);
                
                % First time: initialize data structures
                if count == 0
                    % Get segment types from the first subject
                    segment_types = fieldnames(data.interp_final_topos);
                    diff_types = fieldnames(data.interp_diff_topos);
                    
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
                    if isfield(data.interp_final_topos, segment)
                        all_topos.(segment) = [all_topos.(segment), data.interp_final_topos.(segment)];
                    end
                end
                
                for d = 1:length(diff_types)
                    diff = diff_types{d};
                    if isfield(data.interp_diff_topos, diff)
                        all_diffs.(diff) = [all_diffs.(diff), data.interp_diff_topos.(diff)];
                    end
                end
                
                count = count + 1;
                
            catch err
                warning('Error loading data for subject %s: %s', subject_id, err.message);
            end
        else
            warning('Interpolated data file not found for subject %s: %s', subject_id, subject_data_file);
        end
    end
    
    % If no data loaded, exit
    if count == 0
        warning('No valid interpolated data found for any subjects in %s condition. Skipping group average.', condition);
        group_topos = struct();
        group_diffs = struct();
        return;
    end
    
    fprintf('Successfully loaded data from %d/%d subjects for %s condition\n', ...
        count, length(subjects), condition);
    
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

function createGroupVisualizations(condition, data, group_results_dir, electrode_file, excluded_channels)
    % Creates visualizations using the processed data with consistent color scales
    fprintf('Creating visualizations for %s condition...\n', condition);
    
    global GLOBAL_SEGMENT_MIN GLOBAL_SEGMENT_MAX GLOBAL_DIFF_ABS_MAX;
    
    % Save group average data
    group_data_file = fullfile(group_results_dir, sprintf('group_%s_average_interp.mat', lower(condition)));
    group_topos = data.group_topos;
    group_diffs = data.group_diffs;
    count = data.count;
    save(group_data_file, 'group_topos', 'group_diffs', 'count');
    fprintf('Saved group average data to: %s\n', group_data_file);
    
    % Try to load electrode locations
    try
        EEG = struct();
        EEG.chanlocs = readlocs(electrode_file);
        EEG.nbchan = length(EEG.chanlocs);
        createGroupTopoplots(EEG, group_topos, group_diffs, condition, group_results_dir, excluded_channels, count);
    catch err
        warning('Failed to create topoplots: %s', err.message);
    end
end

function createGroupTopoplots(EEG, group_topos, group_diffs, condition, group_results_dir, excluded_channels, subject_count)
    % Create topoplots for group average data using interpolated data
    fprintf('Creating topoplots for %s group average from interpolated data...\n', condition);
    
    global GLOBAL_SEGMENT_MIN GLOBAL_SEGMENT_MAX GLOBAL_DIFF_ABS_MAX;
    
    % Define standard segment types and their display names
    standard_segments = {'pre_stim', 'early_stim', 'late_stim', 'post_stim'};
    display_names = {'Pre-Stim', 'Early-Stim', 'Late-Stim', 'Post-Stim'};
    
    % Use global color limits for consistent scaling across conditions
    segment_min = GLOBAL_SEGMENT_MIN;
    segment_max = GLOBAL_SEGMENT_MAX;
    diff_limit = [-GLOBAL_DIFF_ABS_MAX GLOBAL_DIFF_ABS_MAX]; % Symmetric colormap for differences
    
    % Create comprehensive figure with segments and differences
    figure('Name', sprintf('%s Group Average', condition), 'Position', [50, 50, 1200, 800], 'visible', 'off');
    
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
                % 'labels', 'efontsize', 6,
                topoplot(topo_data, EEG.chanlocs, 'maplimits', [segment_min segment_max], ...
                    'electrodes', 'on', 'plotchans', included_indices);
                
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
                'electrodes', 'on', 'plotchans', included_indices);
            
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
    
    % Update the title to show correct subject count
    sgtitle(sprintf('%s Group Average (n=%d)', condition, subject_count), 'FontSize', 16);
    
    % Save the comprehensive figure
    saveas(gcf, fullfile(group_results_dir, sprintf('group_%s_topoplots_interp.png', lower(condition))));
    fprintf('Saved group topoplots to: %s\n', fullfile(group_results_dir, sprintf('group_%s_topoplots_interp.png', lower(condition))));
    
    % Create a list of visible channels and save it
    visible_channels = {EEG.chanlocs(included_indices).labels};
    channel_list_file = fullfile(group_results_dir, sprintf('group_%s_visible_channels.txt', lower(condition)));
    fid = fopen(channel_list_file, 'w');
    if fid ~= -1
        fprintf(fid, 'Visible channels for %s group analysis (%d channels):\n\n', condition, length(visible_channels));
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
