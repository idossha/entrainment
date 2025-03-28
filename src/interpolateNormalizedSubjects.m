function interpolateNormalizedSubjects(experiment_path, session_id)
    % INTERPOLATENORMALIZEDSUBJECTS - Interpolates missing channels in normalized ISPC data
    %
    % This function loads normalized ISPC data for each subject, interpolates missing
    % channels using EEGLAB's pop_interp function, and saves the interpolated data.
    %
    % Inputs:
    %   experiment_path - Base path for the experiment data
    %   session_id - Session ID (e.g., 'N1')
    
    % Set default session if not provided
    if nargin < 2
        session_id = 'N1';
    end
    
    fprintf('Interpolating normalized ISPC data for session %s...\n', session_id);
    
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
    
    % Load standard electrode locations
    std_chanlocs = readlocs(electrode_file);
    fprintf('Loaded %d channel locations from standard electrode file\n', length(std_chanlocs));
    
    % Create a list of channels to exclude
    excluded_indices = [];
    for i = 1:length(excluded_channels)
        for c = 1:length(std_chanlocs)
            if strcmpi(std_chanlocs(c).labels, excluded_channels{i})
                excluded_indices = [excluded_indices, c];
                break;
            end
        end
    end
    fprintf('Found %d excluded channels in standard electrode file\n', length(excluded_indices));
    
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
    
    % Process each subject
    subjects_ids = {subjects.id};
    for i = 1:length(subjects_ids)
        subject_id = subjects_ids{i};
        fprintf('\n\n========== Processing Subject %s, Session %s ==========\n\n', subject_id, session_id);
        
        % Define paths
        subject_base_dir = fullfile(experiment_path, subject_id, session_id);
        output_dir = fullfile(subject_base_dir, 'output', sprintf('entrainment_%s_%s', subject_id, session_id));
        
        if ~exist(output_dir, 'dir')
            warning('Output directory not found for subject %s: %s. Skipping.', subject_id, output_dir);
            continue;
        end
        
        % Load EEG file first to get proper channel and chainfo structure
        eeg_file = fullfile(subject_base_dir, sprintf('Strength_%s_%s_forSW.set', subject_id, session_id));
        
        if ~exist(eeg_file, 'file')
            warning('EEG file not found for subject %s: %s. Cannot interpolate without channel information.', subject_id, eeg_file);
            continue;
        end
        
        % Load EEG with channel info
        fprintf('Loading EEG channel info for subject %s...\n', subject_id);
        try
            EEG = pop_loadset('filename', eeg_file, 'loadmode', 'info');
            fprintf('Loaded EEG with %d channels and chaninfo structure\n', EEG.nbchan);
        catch err
            warning('Error loading EEG file: %s. Skipping subject.', err.message);
            continue;
        end
        
        % Look for subject-level normalized data file
        subject_data_file = fullfile(output_dir, sprintf('subject_%s_normalized_average.mat', subject_id));
        
        if exist(subject_data_file, 'file')
            fprintf('Found normalized data file: %s\n', subject_data_file);
            
            try
                % Load subject's normalized data
                data = load(subject_data_file);
                
                % Create interpolated versions of the normalized data
                [interp_norm_topos, interp_pct_topos] = interpolateTopoMaps(data.norm_final_topos, data.pct_change_topos, std_chanlocs, excluded_indices, EEG);
                
                % Save interpolated data
                interp_file = fullfile(output_dir, sprintf('subject_%s_normalized_interp.mat', subject_id));
                save(interp_file, 'interp_norm_topos', 'interp_pct_topos');
                fprintf('Saved interpolated normalized data to: %s\n', interp_file);
            catch err
                warning('Error processing subject %s: %s', subject_id, err.message);
            end
        else
            fprintf('Subject-level normalized data file not found. Checking for individual stimulation results...\n');
        end
        
        % Process individual stimulation protocols
        stim_dirs = dir(fullfile(output_dir, 'stim_*'));
        
        if ~isempty(stim_dirs)
            fprintf('Found %d stimulation directories. Processing each one...\n', length(stim_dirs));
            
            % Look for global ISPC file
            global_file = fullfile(output_dir, 'global_ispc.mat');
            
            if ~exist(global_file, 'file')
                warning('Global ISPC file not found for subject %s. Cannot continue without global normalization value.', subject_id);
                continue;
            end
            
            try
                % Load global ISPC values
                global_data = load(global_file);
                globalISPC = global_data.globalISPC;
                globalStd = global_data.globalStd;
                
                fprintf('Loaded global ISPC: %.4f, std: %.4f\n', globalISPC, globalStd);
                
                % Process each stimulation directory
                for s = 1:length(stim_dirs)
                    stim_dir = fullfile(output_dir, stim_dirs(s).name);
                    results_file = fullfile(stim_dir, 'results.mat');
                    
                    if exist(results_file, 'file')
                        fprintf('Processing %s for subject %s...\n', stim_dirs(s).name, subject_id);
                        
                        try
                            % Load results file
                            stim_data = load(results_file);
                            
                            % Ensure the file contains ISPC results
                            if ~isfield(stim_data, 'ispc_results') || ~isfield(stim_data, 'stim_segments')
                                warning('Results file for %s does not contain required fields. Skipping.', stim_dirs(s).name);
                                continue;
                            end
                            
                            % If normalized data already exists, use it
                            if isfield(stim_data, 'norm_ispc_results')
                                fprintf('Normalized data found in results file. Using existing normalization.\n');
                                norm_ispc_results = stim_data.norm_ispc_results;
                            else
                                % Calculate normalized data
                                fprintf('Calculating normalized data...\n');
                                norm_ispc_results = stim_data.ispc_results / globalISPC;
                            end
                            
                            % For each segment, calculate average and interpolate
                            segment_names = fieldnames(stim_data.stim_segments);
                            interp_segment_means = struct();
                            
                            for seg_idx = 1:length(segment_names)
                                segment = segment_names{seg_idx};
                                segment_range = stim_data.stim_segments.(segment);
                                
                                % Calculate mean normalized ISPC for this segment
                                segment_mean = mean(norm_ispc_results(:, segment_range(1):segment_range(2)), 2, 'omitnan');
                                
                                % Create a proper temporary EEG structure for interpolation
                                temp_EEG = createTempEEG(EEG, segment_mean);
                                
                                % Interpolate to standard locations using EEGLAB
                                interp_EEG = pop_interp(temp_EEG, std_chanlocs, 'spherical');
                                
                                % Set excluded channels to NaN
                                interp_EEG.data(excluded_indices, 1) = NaN;
                                
                                % Store interpolated data
                                interp_segment_means.(segment) = interp_EEG.data(:,1);
                            end
                            
                            % Calculate percent changes
                            interp_pct_change = struct();
                            
                            if isfield(interp_segment_means, 'pre_stim') && isfield(interp_segment_means, 'late_stim')
                                % Calculate late vs pre
                                interp_pct_change.late_vs_pre = ((interp_segment_means.late_stim - interp_segment_means.pre_stim) ./ interp_segment_means.pre_stim) * 100;
                                
                                % Calculate early vs pre (if available)
                                if isfield(interp_segment_means, 'early_stim')
                                    interp_pct_change.early_vs_pre = ((interp_segment_means.early_stim - interp_segment_means.pre_stim) ./ interp_segment_means.pre_stim) * 100;
                                    
                                    % Calculate late vs early
                                    interp_pct_change.late_vs_early = ((interp_segment_means.late_stim - interp_segment_means.early_stim) ./ interp_segment_means.early_stim) * 100;
                                end
                                
                                % Calculate post vs pre (if available)
                                if isfield(interp_segment_means, 'post_stim')
                                    interp_pct_change.post_vs_pre = ((interp_segment_means.post_stim - interp_segment_means.pre_stim) ./ interp_segment_means.pre_stim) * 100;
                                end
                                
                                % Save interpolated data for this stimulation
                                interp_file = fullfile(stim_dir, 'interp_results.mat');
                                save(interp_file, 'interp_segment_means', 'interp_pct_change', 'globalISPC', 'globalStd');
                                fprintf('Saved interpolated data for %s to: %s\n', stim_dirs(s).name, interp_file);
                            else
                                warning('Cannot calculate percent changes: missing pre_stim or late_stim segment.');
                            end
                        catch err
                            warning('Error processing %s for subject %s: %s', stim_dirs(s).name, subject_id, err.message);
                        end
                    else
                        warning('Results file not found for %s, subject %s.', stim_dirs(s).name, subject_id);
                    end
                end
            catch err
                warning('Error processing stimulation protocols for subject %s: %s', subject_id, err.message);
            end
        else
            fprintf('No stimulation directories found for subject %s.\n', subject_id);
        end
    end
    
    fprintf('Interpolation of normalized ISPC data completed.\n');
end

function [interp_norm_topos, interp_pct_topos] = interpolateTopoMaps(norm_topos, pct_topos, std_chanlocs, excluded_indices, EEG)
    % Load standard electrode locations and interpolate topo maps
    fprintf('Interpolating topo maps to standard electrode layout...\n');
    
    % Initialize output structures
    interp_norm_topos = struct();
    interp_pct_topos = struct();
    
    % Get first field to determine channel structure
    norm_fields = fieldnames(norm_topos);
    
    if isempty(norm_fields)
        warning('No fields found in norm_topos. Nothing to interpolate.');
        return;
    end
    
    % Interpolate each field in norm_topos
    for i = 1:length(norm_fields)
        field = norm_fields{i};
        data = norm_topos.(field);
        
        % Create a proper temporary EEG structure for interpolation
        temp_EEG = createTempEEG(EEG, data);
        
        % Interpolate to standard locations
        try
            interp_EEG = pop_interp(temp_EEG, std_chanlocs, 'spherical');
            
            % Set excluded channels to NaN
            interp_EEG.data(excluded_indices, 1) = NaN;
            
            % Store interpolated data
            interp_norm_topos.(field) = interp_EEG.data(:,1);
            
            fprintf('Successfully interpolated %s\n', field);
        catch err
            warning('Error interpolating %s: %s', field, err.message);
            % Fall back to simple approach if interpolation fails
            new_data = zeros(length(std_chanlocs), 1);
            new_data(1:min(length(data), length(std_chanlocs))) = data(1:min(length(data), length(std_chanlocs)));
            new_data(excluded_indices) = NaN;
            interp_norm_topos.(field) = new_data;
        end
    end
    
    % Interpolate each field in pct_topos
    pct_fields = fieldnames(pct_topos);
    
    for i = 1:length(pct_fields)
        field = pct_fields{i};
        data = pct_topos.(field);
        
        % Create a proper temporary EEG structure for interpolation
        temp_EEG = createTempEEG(EEG, data);
        
        % Interpolate to standard locations
        try
            interp_EEG = pop_interp(temp_EEG, std_chanlocs, 'spherical');
            
            % Set excluded channels to NaN
            interp_EEG.data(excluded_indices, 1) = NaN;
            
            % Store interpolated data
            interp_pct_topos.(field) = interp_EEG.data(:,1);
            
            fprintf('Successfully interpolated %s\n', field);
        catch err
            warning('Error interpolating %s: %s', field, err.message);
            % Fall back to simple approach if interpolation fails
            new_data = zeros(length(std_chanlocs), 1);
            new_data(1:min(length(data), length(std_chanlocs))) = data(1:min(length(data), length(std_chanlocs)));
            new_data(excluded_indices) = NaN;
            interp_pct_topos.(field) = new_data;
        end
    end
    
    fprintf('Topo map interpolation complete.\n');
end

function temp_EEG = createTempEEG(original_EEG, data)
    % Create a proper temporary EEG structure for interpolation
    % that includes the chaninfo field and other required fields
    
    temp_EEG = struct();
    temp_EEG.nbchan = length(data);
    
    % Use original EEG for channel locations, trimming if needed
    if temp_EEG.nbchan <= original_EEG.nbchan
        temp_EEG.chanlocs = original_EEG.chanlocs(1:temp_EEG.nbchan);
    else
        % If data has more channels than original EEG (shouldn't happen)
        temp_EEG.chanlocs = original_EEG.chanlocs;
        data = data(1:original_EEG.nbchan); % Trim data to match
        temp_EEG.nbchan = length(temp_EEG.chanlocs);
    end
    
    % Copy the chaninfo field from the original EEG
    if isfield(original_EEG, 'chaninfo')
        temp_EEG.chaninfo = original_EEG.chaninfo;
    else
        % Create a basic chaninfo structure if the original doesn't have one
        temp_EEG.chaninfo = struct();
        temp_EEG.chaninfo.nosedir = '+X';
        temp_EEG.chaninfo.plotrad = [];
        temp_EEG.chaninfo.shrink = [];
    end
    
    % Add data field (required for interpolation)
    temp_EEG.data = zeros(temp_EEG.nbchan, 1);
    temp_EEG.data(:,1) = data;
    
    % Add other required fields
    temp_EEG.pnts = 1;
    temp_EEG.trials = 1;
    temp_EEG.srate = original_EEG.srate;
    temp_EEG.xmin = 0;
    temp_EEG.xmax = 0;
    
    % Add minimal event structure
    temp_EEG.event = struct('type', {}, 'latency', {}, 'urevent', {});
    temp_EEG.urevent = struct('type', {}, 'latency', {});
    
    % Add other fields that might be expected
    temp_EEG.icawinv = [];
    temp_EEG.icasphere = [];
    temp_EEG.icaweights = [];
    temp_EEG.icaact = [];
    temp_EEG.setname = 'temp';
    temp_EEG.filename = '';
    temp_EEG.filepath = '';
    temp_EEG.subject = '';
    temp_EEG.group = '';
    temp_EEG.condition = '';
    temp_EEG.session = [];
    temp_EEG.comments = 'Temporary EEG for interpolation';
    temp_EEG.ref = 'common';
    
    % Make sure all required fields are present
    temp_EEG = eeg_checkset(temp_EEG);
end
