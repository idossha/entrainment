function runGroupAnalysis(experiment_path, session_id)
    % RUNGROUPANALYSIS Performs group-level analysis across subjects
    %
    % This function creates group averages for ACTIVE and SHAM conditions
    %
    % Inputs:
    %   experiment_path - Base path for the experiment data
    %   session_id - Session ID (e.g., 'N1')
    
    % Set default session if not provided
    if nargin < 2
        session_id = 'N1';
    end
    
    fprintf('Running group analysis for session %s...\n', session_id);
    
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
    
    % Create group averages for each condition
    createGroupAverage('ACTIVE', active_subjects, session_id, experiment_path, group_results_dir, electrode_file, excluded_channels);
    createGroupAverage('SHAM', sham_subjects, session_id, experiment_path, group_results_dir, electrode_file, excluded_channels);
    
    fprintf('Group analysis completed successfully.\n');
end

function createGroupAverage(condition, subjects, session_id, experiment_path, group_results_dir, electrode_file, excluded_channels)
    % Create average across subjects within a condition
    fprintf('Creating group average for %s condition (%d subjects)...\n', ...
        condition, length(subjects));
    
    % Initialize data structures
    all_topos = struct();
    all_diffs = struct();
    count = 0;
    
    % Load data from each subject
    for i = 1:length(subjects)
        subject_id = subjects{i};
        subject_results_dir = fullfile(experiment_path, subject_id, session_id, 'output', ...
            sprintf('entrainment_%s_%s', subject_id, session_id));
        subject_data_file = fullfile(subject_results_dir, sprintf('subject_%s_average.mat', subject_id));
        
        if exist(subject_data_file, 'file')
            fprintf('Loading data for subject %s...\n', subject_id);
            
            try
                % Load subject's average data
                data = load(subject_data_file);
                
                % First subject initializes the data structures
                if count == 0
                    % Get segment types from the first subject
                    segment_types = fieldnames(data.final_topos);
                    diff_types = fieldnames(data.diff_topos);
                    
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
                    if isfield(data.final_topos, segment)
                        all_topos.(segment) = [all_topos.(segment), data.final_topos.(segment)];
                    end
                end
                
                for d = 1:length(diff_types)
                    diff = diff_types{d};
                    if isfield(data.diff_topos, diff)
                        all_diffs.(diff) = [all_diffs.(diff), data.diff_topos.(diff)];
                    end
                end
                
                count = count + 1;
                
            catch err
                warning('Error loading data for subject %s: %s', subject_id, err.message);
            end
        else
            warning('Data file not found for subject %s: %s', subject_id, subject_data_file);
        end
    end
    
    % If no data loaded, exit
    if count == 0
        warning('No valid data found for any subjects in %s condition. Skipping group average.', condition);
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
    
    % Save group average data
    group_data_file = fullfile(group_results_dir, sprintf('group_%s_average.mat', lower(condition)));
    save(group_data_file, 'group_topos', 'group_diffs', 'subjects');
    fprintf('Saved group average data to: %s\n', group_data_file);
    
    % Load an EEG file from the first subject to get channel info
    if ~isempty(subjects)
        subject_id = subjects{1};
        eeg_file = fullfile(experiment_path, subject_id, session_id, sprintf('Strength_%s_%s_forSW.set', subject_id, session_id));
        
        try
            if exist(eeg_file, 'file')
                % Load EEG file to get urchanlocs (original unmodified channel locations)
                fprintf('Loading EEG data for channel info from: %s\n', eeg_file);
                EEG = pop_loadset(eeg_file);
                
                % Check if urchanlocs exists and use it instead of chanlocs
                if isfield(EEG, 'urchanlocs') && ~isempty(EEG.urchanlocs)
                    fprintf('Using urchanlocs for more accurate electrode positions\n');
                    % Create a copy of EEG with urchanlocs as chanlocs
                    EEG_with_ur = EEG;
                    EEG_with_ur.chanlocs = EEG.urchanlocs;
                    createGroupTopoplots(EEG_with_ur, group_topos, group_diffs, condition, group_results_dir, excluded_channels, count);
                else
                    % Fallback to regular chanlocs
                    fprintf('urchanlocs not found, using regular chanlocs\n');
                    createGroupTopoplots(EEG, group_topos, group_diffs, condition, group_results_dir, excluded_channels, count);
                end
            else
                warning('Could not find EEG file for subject %s. Skipping topoplots.', subject_id);
            end
        catch err
            warning('Error loading EEG data: %s\nTrying alternative method...', err.message);
            try
                % Try with electrode file instead
                fprintf('Attempting to use electrode file directly: %s\n', electrode_file);
                chanlocs = readlocs(electrode_file);
                EEG = struct('chanlocs', chanlocs, 'nbchan', length(chanlocs));
                
                % Create both types of visualizations
                try
                    % First try normal topoplots
                    createGroupTopoplots(EEG, group_topos, group_diffs, condition, group_results_dir, excluded_channels, count);
                catch err2
                    fprintf('Error creating standard topoplots: %s\n', err2.message);
                    % Fallback to simple visualization
                    createSimpleTopoplots(EEG, group_topos, group_diffs, condition, group_results_dir, excluded_channels, count);
                end
            catch err2
                warning('Failed to create topoplots: %s', err2.message);
            end
        end
    end
end

function createGroupTopoplots(EEG, group_topos, group_diffs, condition, group_results_dir, excluded_channels, subject_count)
    % Create topoplots for group average data (using EEG data from a subject)
    fprintf('Creating topoplots for %s group average using subject EEG data...\n', condition);
    
    % Define standard segment types and their display names
    standard_segments = {'pre_stim', 'early_stim', 'late_stim', 'post_stim'};
    display_names = {'Pre-Stim', 'Early-Stim', 'Late-Stim', 'Post-Stim'};
    
    % Find min/max values for consistent color scales
    segment_values = [];
    for i = 1:length(standard_segments)
        segment_type = standard_segments{i};
        if isfield(group_topos, segment_type)
            segment_values = [segment_values; group_topos.(segment_type)];
        end
    end
    
    segment_min = min(segment_values);
    segment_max = max(segment_values);
    
    % Find min/max values for difference maps
    diff_values = [];
    diff_fields = fieldnames(group_diffs);
    for i = 1:length(diff_fields)
        diff_values = [diff_values; group_diffs.(diff_fields{i})];
    end
    
    diff_abs_max = max(abs(diff_values));
    diff_limit = [-diff_abs_max diff_abs_max]; % Symmetric colormap for differences
    
    % Create comprehensive figure with segments and differences
    figure('Name', sprintf('%s Group Average', condition), 'Position', [50, 50, 1200, 800], 'visible', 'off');
    
    % Calculate total number of plots
    num_segments = sum(isfield(group_topos, standard_segments));
    num_diffs = length(diff_fields);
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
            
            % Make sure data length matches number of channels
            if length(topo_data) ~= length(EEG.chanlocs)
                fprintf('Data length (%d) does not match channel count (%d) for %s\n', ...
                    length(topo_data), length(EEG.chanlocs), segment_type);
                
                % If there's a stim channel, remove it (usually the last channel)
                if length(topo_data) == length(EEG.chanlocs) + 1
                    topo_data = topo_data(1:end-1);
                end
            end
            
            try
                % Create topoplot with only valid channels
                topoplot(topo_data, EEG.chanlocs, 'maplimits', [segment_min segment_max], ...
                    'electrodes', 'on', 'plotchans', included_indices);
                
                title(display_names{i});
                colorbar;
            catch err
                fprintf('Error in topoplot for %s: %s\n', segment_type, err.message);
            end
            
            plot_idx = plot_idx + 1;
        end
    end
    
    % Then plot the differences
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
        
        % Make sure data length matches number of channels
        if length(topo_data) ~= length(EEG.chanlocs)
            fprintf('Data length (%d) does not match channel count (%d) for %s\n', ...
                length(topo_data), length(EEG.chanlocs), diff_type);
            
            % If there's a stim channel, remove it (usually the last channel)
            if length(topo_data) == length(EEG.chanlocs) + 1
                topo_data = topo_data(1:end-1);
            end
        end
        
        try
            % Create topoplot with only valid channels
            topoplot(topo_data, EEG.chanlocs, 'maplimits', diff_limit, ...
                'electrodes', 'on', 'plotchans', included_indices);
            
            title(display_name);
            colorbar;
        catch err
            fprintf('Error in topoplot for diff %s: %s\n', diff_type, err.message);
        end
        
        plot_idx = plot_idx + 1;
    end
    
    % Update the title to show correct subject count
    sgtitle(sprintf('%s Group Average (n=%d)', condition, subject_count), 'FontSize', 16);
    
    % Save the comprehensive figure
    saveas(gcf, fullfile(group_results_dir, sprintf('group_%s_topoplots.png', lower(condition))));
    fprintf('Saved group topoplots to: %s\n', fullfile(group_results_dir, sprintf('group_%s_topoplots.png', lower(condition))));
    
    % Close the figure
    close(gcf);
end

function createSimpleTopoplots(EEG, group_topos, group_diffs, condition, group_results_dir, excluded_channels, subject_count)
    % Create simple topoplots without relying on EEGLAB's topoplot function
    fprintf('Creating simple topoplots for %s group average...\n', condition);
    
    % Define standard segment types and their display names
    standard_segments = {'pre_stim', 'early_stim', 'late_stim', 'post_stim'};
    display_names = {'Pre-Stim', 'Early-Stim', 'Late-Stim', 'Post-Stim'};
    
    % Find min/max values for consistent color scales
    segment_values = [];
    for i = 1:length(standard_segments)
        segment_type = standard_segments{i};
        if isfield(group_topos, segment_type)
            segment_values = [segment_values; group_topos.(segment_type)];
        end
    end
    
    segment_min = min(segment_values);
    segment_max = max(segment_values);
    
    % Find min/max values for difference maps
    diff_values = [];
    diff_fields = fieldnames(group_diffs);
    for i = 1:length(diff_fields)
        diff_values = [diff_values; group_diffs.(diff_fields{i})];
    end
    
    diff_abs_max = max(abs(diff_values));
    diff_limit = [-diff_abs_max diff_abs_max]; % Symmetric colormap for differences
    
    % Create comprehensive figure with segments and differences
    figure('Name', sprintf('%s Group Average', condition), 'Position', [50, 50, 1200, 800], 'visible', 'off');
    
    % Calculate total number of plots
    num_segments = sum(isfield(group_topos, standard_segments));
    num_diffs = length(diff_fields);
    total_plots = num_segments + num_diffs;
    
    % Calculate a suitable grid layout
    grid_cols = ceil(sqrt(total_plots));
    grid_rows = ceil(total_plots / grid_cols);
    
    % Keep track of current subplot index
    plot_idx = 1;
    
    % First plot the segments
    for i = 1:length(standard_segments)
        segment_type = standard_segments{i};
        if isfield(group_topos, segment_type)
            subplot(grid_rows, grid_cols, plot_idx);
            
            % Just plot the data as a bar chart for now
            bar(group_topos.(segment_type));
            
            title([display_names{i} ' (Simple View)']);
            ylabel('ISPC Value');
            xlabel('Channel Index');
            
            plot_idx = plot_idx + 1;
        end
    end
    
    % Then plot the differences
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
        
        % Just plot the data as a bar chart for now
        bar(group_diffs.(diff_type));
        
        title([display_name ' (Simple View)']);
        ylabel('ISPC Value');
        xlabel('Channel Index');
        
        plot_idx = plot_idx + 1;
    end
    
    % Update the title to show correct subject count
    sgtitle(sprintf('%s Group Average (n=%d) - Simple View', condition, subject_count), 'FontSize', 16);
    
    % Save the comprehensive figure
    saveas(gcf, fullfile(group_results_dir, sprintf('group_%s_simple_view.png', lower(condition))));
    fprintf('Saved simple visualization to: %s\n', fullfile(group_results_dir, sprintf('group_%s_simple_view.png', lower(condition))));
    
    % Close the figure
    close(gcf);
end