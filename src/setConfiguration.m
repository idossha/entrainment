function config = setConfiguration(subject_id, session_id, experiment_path)
    % SETCONFIGURATION Creates a configuration structure for EEG entrainment analysis
    %
    % This function creates a configuration structure for a specific subject and session
    %
    % Inputs:
    %   subject_id - Subject ID (e.g., '101')
    %   session_id - Session ID (e.g., 'N1')
    %   experiment_path - Base path for the experiment data
    %
    % Output:
    %   config - Configuration structure with all required parameters
    
    % Create base paths
    code_path = '/Users/idohaber/Git-Projects/entrainment';
    
    % Define paths
    config.eeglab_path = '/Users/idohaber/Documents/MATLAB/eeglab2024.0/';
    config.src_path = fullfile(code_path, 'src');
    config.data_path = fullfile(experiment_path, subject_id, session_id);
    config.module_path = fullfile(code_path, 'src');
    config.assets_path = fullfile(code_path, 'src/assets');
    config.utilities_path = fullfile(code_path, 'src/utils');
    config.exclude_json = fullfile(config.assets_path, 'exclude.json');
    config.regions_json = fullfile(config.assets_path, 'regions.json');
    config.subject_condition_json = fullfile(config.assets_path, 'subject_condition.json');
    
    % Define results directory
    output_dir = fullfile(experiment_path, subject_id, session_id, 'output');
    config.results_dir = fullfile(output_dir, sprintf('entrainment_%s_%s', subject_id, session_id));
    
    % EEG data settings
    config.eeg_filename = fullfile(config.data_path, sprintf('Strength_%s_%s_forSW.set', subject_id, session_id));
    config.target_srate = 50; % Target sampling rate in Hz
    
    % Stimulus parameters
    config.stim_freq = 1; % 1Hz stimulation
    
    % Protocol definition - finds ALL occurrences of these markers
    config.protocol = struct();
    config.protocol.start_marker = 'stim start';
    config.protocol.end_marker = 'stim end';
    config.min_stim_duration = 60; % Minimum stimulation duration in seconds
    config.protocol.segment_duration = 45; % 45 seconds for each segment
    
    % Analysis parameters
    config.window_size = 20;    % 20-second window for ISPC calculation
    config.sliding_step = 20;   % Step size for sliding window (in samples) - default: 1 sample
    config.wavelet_cycles = 8;  % Number of cycles in the Morlet wavelet
    config.post_stim_exclusion = 45; % Time in seconds to exclude after stim end marker


    % Process all stimulation instances separately
    config.process_all_instances = true;
    
    % Electrode positions file
    % Prioritize the assets version of the file, then fall back to utils if needed
    config.electrode_file = fullfile(config.assets_path, 'GSN-HydroCel-256.sfp');
    if ~exist(config.electrode_file, 'file')
        config.electrode_file = fullfile(config.utilities_path, 'egi256_GSN_HydroCel.sfp');
        if ~exist(config.electrode_file, 'file')
            warning('Neither electrode position file was found. This may cause issues with topoplot visualizations.');
        end
    end
    
    % Create necessary directories if they don't exist
    if ~exist(output_dir, 'dir')
        [success, msg] = mkdir(output_dir);
        if ~success
            error('Failed to create output directory %s: %s', output_dir, msg);
        end
        fprintf('Created output directory: %s\n', output_dir);
    end
    
    if ~exist(config.results_dir, 'dir')
        [success, msg] = mkdir(config.results_dir);
        if ~success
            error('Failed to create results directory %s: %s', config.results_dir, msg);
        end
        fprintf('Created results directory: %s\n', config.results_dir);
    end
   
    % Parse JSON configuration files
    config = parseJSONConfig(config);
end
