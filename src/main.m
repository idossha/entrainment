%% Main script for EEG entrainment analysis across multiple subjects
% This script analyzes the entrainment of EEG data to an exogenous electrical stimulus
% for multiple subjects and sessions.

clear all; close all;

%% Define experiment settings
% Base experiment path
experiment_path = '/Volumes/Ido/analyze';

% Define subjects and nights
subjects = {'101', '102','107','108','109','110','111','112','114','115','116','117','119','120','121','122','127','132'};
% subjects = {'112','122','127','132'};
nights = {'N1'};

%% Add EEGLAB path and initialize
eeglab_path = '/Users/idohaber/Documents/MATLAB/eeglab2024.0/';
addpath(eeglab_path);
eeglab nogui;

%% Loop through each subject and night
for subjIdx = 1:length(subjects)
    for nightIdx = 1:length(nights)
        try
            %% Define Current Subject and Night
            whichSubj = subjects{subjIdx};
            whichSess = nights{nightIdx};
            
            % Display current subject/session
            fprintf('\n\n========== Processing Subject %s, Session %s ==========\n\n', whichSubj, whichSess);
            
            %% Check if file exists
            eeg_filename = fullfile(experiment_path, whichSubj, whichSess, sprintf('Strength_%s_%s_forSW.set', whichSubj, whichSess));
            if ~exist(eeg_filename, 'file')
                fprintf('File not found: %s. Skipping this subject/session.\n', eeg_filename);
                continue;
            end
            
            %% Set configuration for this subject/session
            config = setConfiguration(whichSubj, whichSess, experiment_path);
            
            % Add required paths (from the config)
            addpath(config.src_path);
            addpath(config.module_path);
            addpath(config.utilities_path);
            
            %% Execute analysis pipeline for this subject/session
            
            % Load EEG data and possibly resample
            fprintf('Loading EEG data for subject %s, session %s\n', whichSubj, whichSess);
            EEG = loadEEGData(config.eeg_filename, config.target_srate);
            
            % Find all stimulation events and define segments
            fprintf('Finding stimulation events for subject %s, session %s\n', whichSubj, whichSess);
            [segments, stim_samples] = findEvents(EEG, config);
            
            % Generate stimulus signal and add it to the EEG data
            fprintf('Adding stimulus signal for subject %s, session %s\n', whichSubj, whichSess);
            [EEG, stim_channel_idx] = addStimulusSignal(EEG, config.stim_freq);
            
            % Determine how many stimulation periods we found
            num_stims = size(stim_samples, 1);
            
            % If we found multiple stimulation instances and are processing them separately
            if num_stims > 1 && config.process_all_instances
                fprintf('Processing %d stimulation instances separately for subject %s, session %s\n', ...
                    num_stims, whichSubj, whichSess);
                
                % Get all segment names
                all_segment_names = fieldnames(segments);
                
                % Store results for each stimulation
                all_ispc_results = cell(num_stims, 1);
                all_stim_segments = cell(num_stims, 1);
                
                for stim_idx = 1:num_stims
                    fprintf('\n===== Processing Stimulation %d/%d for Subject %s, Session %s =====\n', ...
                        stim_idx, num_stims, whichSubj, whichSess);
                    
                    % Create directory for this stimulation instance
                    stim_dir = fullfile(config.results_dir, ['stim_' num2str(stim_idx)]);
                    if ~exist(stim_dir, 'dir')
                        [success, msg] = mkdir(stim_dir);
                        if ~success
                            error('Failed to create directory %s: %s', stim_dir, msg);
                        end
                        fprintf('Created directory: %s\n', stim_dir);
                    end
                    
                    % Find segments for this stimulation
                    stim_segments = struct();
                    for i = 1:length(all_segment_names)
                        segment_name = all_segment_names{i};
                        if endsWith(segment_name, ['_' num2str(stim_idx)]) || ...
                           (~contains(segment_name, '_') && stim_idx == 1)
                            % Extract base segment name (without numerical suffix)
                            if contains(segment_name, '_')
                                parts = strsplit(segment_name, '_');
                                base_name = strjoin(parts(1:end-1), '_');
                            else
                                base_name = segment_name;
                            end
                            
                            % Store with the base name for cleaner visualization
                            stim_segments.(base_name) = segments.(segment_name);
                        end
                    end
                    
                    % Extract samples for just this stimulation
                    this_stim_samples = stim_samples(stim_idx, :);
                    
                    % Temporarily update config to store results in stimulation subdirectory
                    stim_config = config;
                    stim_config.results_dir = stim_dir;
                    
                    % Define the time range for this stimulation (from pre-stim start to post-stim end)
                    all_segment_bounds = zeros(length(fieldnames(stim_segments)), 2);
                    segment_names = fieldnames(stim_segments);
                    for seg_idx = 1:length(segment_names)
                        all_segment_bounds(seg_idx, :) = stim_segments.(segment_names{seg_idx});
                    end
                    
                    % Find the overall time range for this stimulation
                    time_range = [min(all_segment_bounds(:,1)), max(all_segment_bounds(:,2))];
                    
                    % Calculate ISPC between each EEG channel and the stimulus (only for relevant time range)
                    [ispc_results, phases, filtered_data] = calculateISPC(EEG, stim_channel_idx, time_range, stim_config);
                    
                    % Store results for this stimulation
                    all_ispc_results{stim_idx} = ispc_results;
                    all_stim_segments{stim_idx} = stim_segments;
                    
                    % Create visualizations for this stimulation
                    visualizeISPCOverview(EEG, ispc_results, phases, filtered_data, stim_segments, this_stim_samples, stim_channel_idx, stim_config);
                    visualizeISPCTopoplots(EEG, ispc_results, stim_segments, stim_config);
                    visualizeISPCTimecourses(EEG, ispc_results, stim_segments, stim_config);
                    
                    % Save results for this stimulation
                    results_file = fullfile(stim_dir, 'results.mat');
                    try
                        save(results_file, 'ispc_results', 'stim_segments', 'this_stim_samples', '-v7.3');
                        fprintf('Saved results to: %s\n', results_file);
                    catch err
                        warning('Failed to save results: %s', err.message);
                    end
                    
                    % Clean up figure windows
                    close all;
                end
                
                % Create summary visualizations comparing all stimulations
                createStimulationComparison(EEG, all_ispc_results, all_stim_segments, config);
                
                % Create aggregate topoplots showing conditions and differences collapsed across protocols
                % Save subject-level average data for group analysis
                config.subject_id = whichSubj;
                [final_topos, diff_topos] = createAggregateTopoplots(EEG, all_ispc_results, all_stim_segments, config, true);
                
                fprintf('\nAnalysis complete for subject %s, session %s. Results saved in %s\n', ...
                    whichSubj, whichSess, config.results_dir);
            else
                % Process all stimulations together (original behavior)
                fprintf('Processing all stimulation instances together for subject %s, session %s\n', ...
                    whichSubj, whichSess);
                
                % Calculate ISPC between each EEG channel and the stimulus
                [ispc_results, phases, filtered_data] = calculateISPC(EEG, stim_channel_idx, config);
                
                % Create visualizations
                visualizeISPCOverview(EEG, ispc_results, phases, filtered_data, segments, stim_samples, stim_channel_idx, config);
                visualizeISPCTopoplots(EEG, ispc_results, segments, config);
                visualizeISPCTimecourses(EEG, ispc_results, segments, config);
                
                % Save results
                save(fullfile(config.results_dir, 'results.mat'), 'ispc_results', 'segments', 'stim_samples', '-v7.3');
                
                fprintf('\nAnalysis complete for subject %s, session %s. Results saved to %s\n', ...
                    whichSubj, whichSess, config.results_dir);
            end
            
        catch err
            % Log any errors but continue with next subject/session
            fprintf('\nERROR processing subject %s, session %s: %s\n', ...
                whichSubj, whichSess, err.message);
            fprintf('Stack trace:\n');
            disp(err.stack);
        end
    end
end

fprintf('\n\n========== Analysis complete for all subjects ==========\n\n');

% Display information about utility functions
fprintf('\nUtility Functions Available:\n');
fprintf('1. visualizeISPCSteps(EEG, timeRange, config)\n');
fprintf('   - Creates a detailed visualization of ISPC calculation steps\n');
fprintf('   - Example: visualizeISPCSteps(EEG, [30, 40], config);\n\n');
fprintf('2. createISPCAnimation(EEG, timeRange, outputFile, config)\n');
fprintf('   - Creates an animation of phase synchronization\n');
fprintf('   - Example: createISPCAnimation(EEG, [30, 40], ''my_animation.gif'', config);\n\n');
fprintf('3. visualizeWindowSliding(EEG, timeRange, config)\n');
fprintf('   - Visualizes how the sliding window affects ISPC calculation\n');
fprintf('   - Example: visualizeWindowSliding(EEG, [30, 35], config);\n\n');
