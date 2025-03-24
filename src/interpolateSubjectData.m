function interpolateSubjectData(experiment_path, session_id)
    % INTERPOLATESUBJECTDATA - Interpolate bad channels for subject average data
    %
    % This function interpolates bad channels in subject_*_average.mat files
    % and creates new files with interpolated data that all have the same
    % channel count/structure for easier group analysis.
    %
    % Inputs:
    %   experiment_path - Base path for the experiment data
    %   session_id - Session ID (e.g., 'N1')
    
    % Set default session if not provided
    if nargin < 2
        session_id = 'N1';
    end
    
    fprintf('Interpolating bad channels for subjects in session %s...\n', session_id);
    
    % Get base paths
    code_path = '/Users/idohaber/Git-Projects/entrainment';
    src_path = fullfile(code_path, 'src');
    assets_path = fullfile(code_path, 'src/assets');
    eeglab_path = fullfile(code_path, 'src/eeglab');
    
    % Add paths
    addpath('/Users/idohaber/Documents/MATLAB/eeglab2024.0/');
    addpath(src_path);
    addpath(eeglab_path);
    
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
    config.utilities_path = fullfile(code_path, 'src/utils');
    
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
    
    % Extract all subject IDs
    all_subjects = {};
    for i = 1:length(subjects)
        all_subjects{end+1} = subjects(i).id;
    end
    
    % Create EEG template with standard electrode locations
    try
        % Read the standard electrode locations
        chanlocs = readlocs(electrode_file);
        fprintf('Loaded %d channels from electrode file\n', length(chanlocs));
        
        % Create a template EEG structure
        template_EEG = struct();
        template_EEG.chanlocs = chanlocs;
        template_EEG.nbchan = length(chanlocs);
        
        % Process each subject
        for i = 1:length(all_subjects)
            subject_id = all_subjects{i};
            processSubject(subject_id, session_id, experiment_path, template_EEG, excluded_channels);
        end
        
    catch err
        fprintf('Error creating template EEG: %s\n', err.message);
    end
    
    fprintf('Completed interpolation for all subjects.\n');
end

function processSubject(subject_id, session_id, experiment_path, template_EEG, excluded_channels)
    % Process a single subject's data
    subject_results_dir = fullfile(experiment_path, subject_id, session_id, 'output', ...
        sprintf('entrainment_%s_%s', subject_id, session_id));
    subject_data_file = fullfile(subject_results_dir, sprintf('subject_%s_average.mat', subject_id));
    
    if ~exist(subject_data_file, 'file')
        fprintf('Subject %s: No average data file found. Skipping.\n', subject_id);
        return;
    end
    
    fprintf('\n--- Processing Subject %s ---\n', subject_id);
    
    try
        % Load the subject's average data
        fprintf('Loading data from: %s\n', subject_data_file);
        data = load(subject_data_file);
        
        % Check if data has the expected fields
        if ~isfield(data, 'final_topos') || ~isfield(data, 'diff_topos')
            fprintf('Subject %s: Data file does not contain expected fields. Skipping.\n', subject_id);
            return;
        end
        
        % Get fields to interpolate
        segment_fields = fieldnames(data.final_topos);
        diff_fields = fieldnames(data.diff_topos);
        
        % Create dummy EEG structure for this subject's data
        EEG = struct();
        
        % Process segment data
        interp_final_topos = struct();
        for s = 1:length(segment_fields)
            segment = segment_fields{s};
            fprintf('  Interpolating %s segment\n', segment);
            
            % Get the original data
            original_data = data.final_topos.(segment);
            n_channels = length(original_data);
            
            % Create a temporary EEG structure for interpolation
            EEG.data = zeros(n_channels, 2); % Just need non-empty data
            EEG.data(1:n_channels, 1) = original_data; % Use the data for channel 1
            
            % Try to load actual EEG data for channel locations
            eeg_file = fullfile(experiment_path, subject_id, session_id, sprintf('Strength_%s_%s_forSW.set', subject_id, session_id));
            
            if exist(eeg_file, 'file')
                try
                    % Load just the channel info
                    temp = load(eeg_file, '-mat'); % Load the .set file in MATLAB format
                    if isfield(temp, 'chanlocs')
                        EEG.chanlocs = temp.chanlocs(1:n_channels); % Use the same channel count
                    else
                        % Use template channel locations but limit to current data size
                        EEG.chanlocs = template_EEG.chanlocs(1:min(n_channels, length(template_EEG.chanlocs)));
                    end
                catch err
                    fprintf('  Warning: Could not load EEG channel info. Using template. Error: %s\n', err.message);
                    EEG.chanlocs = template_EEG.chanlocs(1:min(n_channels, length(template_EEG.chanlocs)));
                end
            else
                % Use template channel locations but limit to current data size
                EEG.chanlocs = template_EEG.chanlocs(1:min(n_channels, length(template_EEG.chanlocs)));
            end
            
            EEG.nbchan = length(EEG.chanlocs);
            
            % Now interpolate to the full template
            try
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
                
                % Create list of channels to interpolate (those not in original data)
                orig_labels = {EEG.chanlocs.labels};
                to_interp = [];
                
                for c = 1:length(template_EEG.chanlocs)
                    % Skip excluded channels
                    if ismember(c, excluded_indices)
                        continue;
                    end
                    
                    % Check if this channel is in the original data
                    if ~ismember(template_EEG.chanlocs(c).labels, orig_labels)
                        to_interp = [to_interp, c];
                    end
                end
                
                fprintf('  Interpolating %d channels to match template\n', length(to_interp));
                
                % Now do the interpolation using pop_interp
                % This is a simplified version since we have average data, not EEG data
                if ~isempty(to_interp)
                    % Create a temporary copy of template with our data inserted
                    temp_EEG = template_EEG;
                    temp_EEG.data = zeros(length(template_EEG.chanlocs), 1);
                    
                    % Map original data to temp_EEG
                    for c = 1:length(EEG.chanlocs)
                        label = EEG.chanlocs(c).labels;
                        for t = 1:length(template_EEG.chanlocs)
                            if strcmpi(template_EEG.chanlocs(t).labels, label)
                                temp_EEG.data(t) = original_data(c);
                                break;
                            end
                        end
                    end
                    
                    % Simple nearest-neighbor interpolation for missing channels
                    % Find the 5 nearest neighbors for each channel to interpolate
                    interp_data = temp_EEG.data;
                    for i = 1:length(to_interp)
                        ch_idx = to_interp(i);
                        
                        % Skip if this is an excluded channel
                        if ismember(ch_idx, excluded_indices)
                            continue;
                        end
                        
                        % Get the coordinates of this channel
                        x1 = template_EEG.chanlocs(ch_idx).X;
                        y1 = template_EEG.chanlocs(ch_idx).Y;
                        z1 = template_EEG.chanlocs(ch_idx).Z;
                        
                        % Calculate distances to all other channels
                        distances = zeros(length(template_EEG.chanlocs), 1);
                        for j = 1:length(template_EEG.chanlocs)
                            % Skip if this is a channel to interpolate or an excluded channel
                            if ismember(j, to_interp) || ismember(j, excluded_indices)
                                distances(j) = Inf;
                                continue;
                            end
                            
                            % Calculate Euclidean distance
                            x2 = template_EEG.chanlocs(j).X;
                            y2 = template_EEG.chanlocs(j).Y;
                            z2 = template_EEG.chanlocs(j).Z;
                            
                            distances(j) = sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2);
                        end
                        
                        % Find the 5 closest channels
                        [~, nearest_idx] = sort(distances);
                        nearest_idx = nearest_idx(1:min(5, length(nearest_idx)));
                        
                        % Remove channels with infinite distance
                        nearest_idx = nearest_idx(~isinf(distances(nearest_idx)));
                        
                        if isempty(nearest_idx)
                            % No neighbors found, skip this channel
                            continue;
                        end
                        
                        % Calculate inverse-distance-weighted average
                        weights = 1 ./ distances(nearest_idx);
                        weights = weights / sum(weights);
                        
                        interp_data(ch_idx) = sum(weights .* temp_EEG.data(nearest_idx));
                    end
                    
                    % Set excluded channels to NaN
                    interp_data(excluded_indices) = NaN;
                    
                    % Store the interpolated data
                    interp_final_topos.(segment) = interp_data;
                else
                    % No channels to interpolate
                    interp_final_topos.(segment) = original_data;
                end
            catch err
                fprintf('  Error interpolating segment %s: %s\n', segment, err.message);
                % Just use the original data if interpolation fails
                interp_final_topos.(segment) = original_data;
            end
        end
        
        % Process difference data
        interp_diff_topos = struct();
        for d = 1:length(diff_fields)
            diff_type = diff_fields{d};
            fprintf('  Interpolating %s difference\n', diff_type);
            
            % Get the original data
            original_data = data.diff_topos.(diff_type);
            n_channels = length(original_data);
            
            % Create a temporary EEG structure for interpolation
            EEG.data = zeros(n_channels, 2); % Just need non-empty data
            EEG.data(1:n_channels, 1) = original_data; % Use the data for channel 1
            
            % Try to load actual EEG data for channel locations
            eeg_file = fullfile(experiment_path, subject_id, session_id, sprintf('Strength_%s_%s_forSW.set', subject_id, session_id));
            
            if exist(eeg_file, 'file')
                try
                    % Load just the channel info
                    temp = load(eeg_file, '-mat'); % Load the .set file in MATLAB format
                    if isfield(temp, 'chanlocs')
                        EEG.chanlocs = temp.chanlocs(1:n_channels); % Use the same channel count
                    else
                        % Use template channel locations but limit to current data size
                        EEG.chanlocs = template_EEG.chanlocs(1:min(n_channels, length(template_EEG.chanlocs)));
                    end
                catch err
                    fprintf('  Warning: Could not load EEG channel info. Using template. Error: %s\n', err.message);
                    EEG.chanlocs = template_EEG.chanlocs(1:min(n_channels, length(template_EEG.chanlocs)));
                end
            else
                % Use template channel locations but limit to current data size
                EEG.chanlocs = template_EEG.chanlocs(1:min(n_channels, length(template_EEG.chanlocs)));
            end
            
            EEG.nbchan = length(EEG.chanlocs);
            
            % Now interpolate to the full template
            try
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
                
                % Create list of channels to interpolate (those not in original data)
                orig_labels = {EEG.chanlocs.labels};
                to_interp = [];
                
                for c = 1:length(template_EEG.chanlocs)
                    % Skip excluded channels
                    if ismember(c, excluded_indices)
                        continue;
                    end
                    
                    % Check if this channel is in the original data
                    if ~ismember(template_EEG.chanlocs(c).labels, orig_labels)
                        to_interp = [to_interp, c];
                    end
                end
                
                fprintf('  Interpolating %d channels to match template\n', length(to_interp));
                
                % Now do the interpolation using pop_interp
                % This is a simplified version since we have average data, not EEG data
                if ~isempty(to_interp)
                    % Create a temporary copy of template with our data inserted
                    temp_EEG = template_EEG;
                    temp_EEG.data = zeros(length(template_EEG.chanlocs), 1);
                    
                    % Map original data to temp_EEG
                    for c = 1:length(EEG.chanlocs)
                        label = EEG.chanlocs(c).labels;
                        for t = 1:length(template_EEG.chanlocs)
                            if strcmpi(template_EEG.chanlocs(t).labels, label)
                                temp_EEG.data(t) = original_data(c);
                                break;
                            end
                        end
                    end
                    
                    % Simple nearest-neighbor interpolation for missing channels
                    % Find the 5 nearest neighbors for each channel to interpolate
                    interp_data = temp_EEG.data;
                    for i = 1:length(to_interp)
                        ch_idx = to_interp(i);
                        
                        % Skip if this is an excluded channel
                        if ismember(ch_idx, excluded_indices)
                            continue;
                        end
                        
                        % Get the coordinates of this channel
                        x1 = template_EEG.chanlocs(ch_idx).X;
                        y1 = template_EEG.chanlocs(ch_idx).Y;
                        z1 = template_EEG.chanlocs(ch_idx).Z;
                        
                        % Calculate distances to all other channels
                        distances = zeros(length(template_EEG.chanlocs), 1);
                        for j = 1:length(template_EEG.chanlocs)
                            % Skip if this is a channel to interpolate or an excluded channel
                            if ismember(j, to_interp) || ismember(j, excluded_indices)
                                distances(j) = Inf;
                                continue;
                            end
                            
                            % Calculate Euclidean distance
                            x2 = template_EEG.chanlocs(j).X;
                            y2 = template_EEG.chanlocs(j).Y;
                            z2 = template_EEG.chanlocs(j).Z;
                            
                            distances(j) = sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2);
                        end
                        
                        % Find the 5 closest channels
                        [~, nearest_idx] = sort(distances);
                        nearest_idx = nearest_idx(1:min(5, length(nearest_idx)));
                        
                        % Remove channels with infinite distance
                        nearest_idx = nearest_idx(~isinf(distances(nearest_idx)));
                        
                        if isempty(nearest_idx)
                            % No neighbors found, skip this channel
                            continue;
                        end
                        
                        % Calculate inverse-distance-weighted average
                        weights = 1 ./ distances(nearest_idx);
                        weights = weights / sum(weights);
                        
                        interp_data(ch_idx) = sum(weights .* temp_EEG.data(nearest_idx));
                    end
                    
                    % Set excluded channels to NaN
                    interp_data(excluded_indices) = NaN;
                    
                    % Store the interpolated data
                    interp_diff_topos.(diff_type) = interp_data;
                else
                    % No channels to interpolate
                    interp_diff_topos.(diff_type) = original_data;
                end
            catch err
                fprintf('  Error interpolating difference %s: %s\n', diff_type, err.message);
                % Just use the original data if interpolation fails
                interp_diff_topos.(diff_type) = original_data;
            end
        end
        
        % Save the interpolated data
        interp_file = fullfile(subject_results_dir, sprintf('subject_%s_interp.mat', subject_id));
        save(interp_file, 'interp_final_topos', 'interp_diff_topos', 'template_EEG');
        fprintf('Saved interpolated data to: %s\n', interp_file);
        
    catch err
        fprintf('Error processing subject %s: %s\n', subject_id, err.message);
    end
end