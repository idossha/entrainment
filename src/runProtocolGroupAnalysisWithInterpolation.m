function runProtocolGroupAnalysisWithInterpolation(experiment_path, session_id, max_stim)
    % RUNPROTOCOLGROUPANALYSISWITHINTERPOLATION - Creates interpolated group averages for specific protocols
    %
    % This function creates interpolated group averages for ACTIVE and SHAM conditions
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
    
    fprintf('Running protocol-specific group analysis with interpolated data for session %s (up to Stim_%d)...\n', ...
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
    
    % Define electrode file
    electrode_file = fullfile(code_path, 'src/utils/egi256_GSN_HydroCel.sfp');
    fprintf('Using electrode locations from: %s\n', electrode_file);
    
    % Define base group results directory
    base_group_dir = fullfile(experiment_path, 'group_results', session_id, 'protocols_normalized_interp');
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
    
    % Process each stimulation protocol
    for stim_num = 1:max_stim
        stim_id = sprintf('stim_%d', stim_num);
        fprintf('\n==== Processing %s with Interpolated Data ====\n', stim_id);
        
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
        
        % Process data for both conditions
        processProtocolWithInterpolation('ACTIVE', active_subjects, session_id, experiment_path, stim_group_dir, electrode_file, excluded_channels, stim_num, stim_id);
        processProtocolWithInterpolation('SHAM', sham_subjects, session_id, experiment_path, stim_group_dir, electrode_file, excluded_channels, stim_num, stim_id);
    end
    
    fprintf('\nProtocol-specific group analysis with interpolated data completed successfully.\n');
end

function processProtocolWithInterpolation(condition, subjects, session_id, experiment_path, stim_group_dir, electrode_file, excluded_channels, stim_num, stim_id)
    % Process protocol-specific data with interpolation for a condition
    fprintf('Processing interpolated data for %s condition, %s (%d subjects)...\n', condition, stim_id, length(subjects));
    
    % Initialize data structures
    all_segment_means = struct();
    all_pct_changes = struct();
    count = 0;
    
    % Define standard segment types
    standard_segments = {'pre_stim', 'early_stim', 'late_stim', 'post_stim'};
    
    % Initialize segment data
    for i = 1:length(standard_segments)
        all_segment_means.(standard_segments{i}) = [];
    end
    
    % Define standard percent changes
    change_types = {'early_vs_pre', 'late_vs_pre', 'post_vs_pre', 'late_vs_early'};
    
    % Initialize percent change data
    for i = 1:length(change_types)
        all_pct_changes.(change_types{i}) = [];
    end
    
    % Load data from each subject
    for i = 1:length(subjects)
        subject_id = subjects{i};
        stim_dir = fullfile(experiment_path, subject_id, session_id, 'output', ...
            sprintf('entrainment_%s_%s', subject_id, session_id), stim_id);
        interp_file = fullfile(stim_dir, 'interp_results.mat');
        
        if exist(interp_file, 'file')
            fprintf('Loading interpolated data for subject %s, %s...\n', subject_id, stim_id);
            
            try
                % Load subject's interpolated stimulation data
                data = load(interp_file);
                
                % Ensure the file contains required fields
                if ~isfield(data, 'interp_segment_means') || ~isfield(data, 'interp_pct_change')
                    warning('Interpolated file for %s, %s does not contain required fields. Skipping.', subject_id, stim_id);
                    continue;
                end
                
                % Add segment means to accumulated data
                for s = 1:length(standard_segments)
                    segment_type = standard_segments{s};
                    if isfield(data.interp_segment_means, segment_type)
                        all_segment_means.(segment_type) = [all_segment_means.(segment_type), data.interp_segment_means.(segment_type)];
                    end
                end
                
                % Add percent changes to accumulated data
                for c = 1:length(change_types)
                    change_type = change_types{c};
                    if isfield(data.interp_pct_change, change_type)
                        all_pct_changes.(change_type) = [all_pct_changes.(change_type), data.interp_pct_change.(change_type)];
                    end
                end
                
                count = count + 1;
                
            catch err
                warning('Error loading data for subject %s, %s: %s', subject_id, stim_id, err.message);
            end
        else
            fprintf('Interpolated file not found for subject %s, %s: %s\n', subject_id, stim_id, interp_file);
        end
    end
    
    % If no data loaded, exit
    if count == 0
        warning('No valid interpolated data found for any subjects in %s condition for %s. Skipping group average.', ...
            condition, stim_id);
        return;
    end
    
    fprintf('Successfully loaded interpolated data from %d/%d subjects for %s condition, %s\n', ...
        count, length(subjects), condition, stim_id);
    
    % Calculate group averages
    group_segment_means = struct();
    group_pct_changes = struct();
    
    segment_types = fieldnames(all_segment_means);
    for s = 1:length(segment_types)
        segment = segment_types{s};
        if ~isempty(all_segment_means.(segment))
            group_segment_means.(segment) = mean(all_segment_means.(segment), 2, 'omitnan');
        end
    end
    
    change_types = fieldnames(all_pct_changes);
    for c = 1:length(change_types)
        change = change_types{c};
        if ~isempty(all_pct_changes.(change))
            group_pct_changes.(change) = mean(all_pct_changes.(change), 2, 'omitnan');
        end
    end
    
    % Save group average data
    group_data_file = fullfile(stim_group_dir, sprintf('group_%s_interpolated.mat', lower(condition)));
    save(group_data_file, 'group_segment_means', 'group_pct_changes', 'count');
    fprintf('Saved interpolated group data to: %s\n', group_data_file);
    
    % Create visualizations
    try
        % Read electrode locations
        chanlocs = readlocs(electrode_file);
        EEG = struct('chanlocs', chanlocs, 'nbchan', length(chanlocs));
        
        createProtocolInterpolatedTopoplots(EEG, group_segment_means, group_pct_changes, condition, stim_group_dir, excluded_channels, count, stim_id);
    catch err
        warning('Failed to create topoplots: %s', err.message);
    end
end

function createProtocolInterpolatedTopoplots(EEG, group_segment_means, group_pct_changes, condition, group_results_dir, excluded_channels, subject_count, stim_id)
    % Create topoplots for interpolated protocol-specific group average data
    fprintf('Creating interpolated topoplots for %s group average (%s)...\n', condition, stim_id);
    
    % Define standard segment types and their display names
    standard_segments = {'pre_stim', 'early_stim', 'late_stim', 'post_stim'};
    display_names = {'Pre-Stim', 'Early-Stim', 'Late-Stim', 'Post-Stim'};
    
    % Define color limits for normalized topoplots
    norm_limits = [0.5 1.5];
    pct_change_limits = [-50 50];
    
    % Create comprehensive figure with segments and percent changes
    figure('Name', sprintf('%s Interpolated Group Average - %s', condition, stim_id), 'Position', [50, 50, 1200, 800], 'visible', 'off');
    
    % Calculate total number of plots
    num_segments = sum(isfield(group_segment_means, standard_segments));
    num_changes = length(fieldnames(group_pct_changes));
    total_plots = num_segments + num_changes;
    
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
    
    % First plot the normalized segments
    for i = 1:length(standard_segments)
        segment_type = standard_segments{i};
        if isfield(group_segment_means, segment_type)
            subplot(grid_rows, grid_cols, plot_idx);
            
            % Prepare data for topoplot
            topo_data = group_segment_means.(segment_type);
            
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
    sgtitle(sprintf('%s Interpolated Group Average - %s (n=%d)', condition, stim_id, subject_count), 'FontSize', 16);
    
    % Save the comprehensive figure
    saveas(gcf, fullfile(group_results_dir, sprintf('group_%s_interpolated_topoplots.png', lower(condition))));
    fprintf('Saved interpolated group topoplots to: %s\n', fullfile(group_results_dir, sprintf('group_%s_interpolated_topoplots.png', lower(condition))));
    
    % Close the figure
    close(gcf);
end
