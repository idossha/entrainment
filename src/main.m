%% Main script for EEG entrainment analysis across multiple subjects
% This script analyzes the entrainment of EEG data to an exogenous electrical stimulus
% for multiple subjects and sessions using normalized ISPC calculation approach

clear all; close all;

%% Define experiment settings
experiment_path = '/Volumes/Ido/analyze';
subjects = { '102','107','108','109','110','111','112','114','115','116','117','119','120','121','122','127','132'};
% subjects = {'101'};
nights = {'N1'};

%% Add EEGLAB path and initialize
eeglab_path = '/Users/idohaber/Documents/MATLAB/eeglab2024.0/';
addpath(eeglab_path);
eeglab nogui;

%% Process each subject
for subjIdx = 1:length(subjects)
    for nightIdx = 1:length(nights)
        try
            whichSubj = subjects{subjIdx};
            whichSess = nights{nightIdx};
            fprintf('\n\n========== Processing Subject %s, Session %s ==========\n\n', whichSubj, whichSess);
            
            eeg_filename = fullfile(experiment_path, whichSubj, whichSess, sprintf('Strength_%s_%s_forSW.set', whichSubj, whichSess));

            if ~exist(eeg_filename, 'file')
                fprintf('File not found: %s. Skipping this subject/session.\n', eeg_filename);
                continue;
            end
            
            config = setConfiguration(whichSubj, whichSess, experiment_path);
            addpath(config.src_path, config.module_path, config.utilities_path);
            
            % load data, downsample if needed and add stimulus signal
            EEG = loadEEGData(config.eeg_filename, config.target_srate);
            [segments, stim_samples] = findEvents(EEG, config);
            [EEG, stim_channel_idx] = addStimulusSignal(EEG, config.stim_freq);
            
             % Calculate global ISPC, channel ISPC, time ISPC, and save with histograms
            [globalISPC, channelISPC, timeISPC, globalStd, nonStimGlobalISPC, nonStimChannelISPC, nonStimTimeISPC, nonStimGlobalStd] = calculateGlobalISPC(EEG, stim_channel_idx, config);
            
            num_stims = size(stim_samples, 1);
            fprintf('Processing %d stimulation instances for subject %s, session %s\n', num_stims, whichSubj, whichSess);
            
            % Get all segment names and initialize storage
            all_segment_names = fieldnames(segments);
            all_ispc_results = cell(num_stims, 1);
            all_norm_ispc_results = cell(num_stims, 1);
            all_z_scored_ispc = cell(num_stims, 1);
            all_stim_segments = cell(num_stims, 1);
            
            % Keep track of valid stimulation indices
            valid_stim_indices = [];
            
            % Process each stimulation separately
            for stim_idx = 1:num_stims
                fprintf('\n===== Processing Stimulation %d/%d =====\n', stim_idx, num_stims);
                
                % Check if the stimulation is long enough
                stim_start = stim_samples(stim_idx, 1);
                stim_end = stim_samples(stim_idx, 2);
                stim_duration_samples = stim_end - stim_start;
                stim_duration_seconds = stim_duration_samples / EEG.srate;
                
                % Skip this stimulation if it's too short
                if stim_duration_seconds < config.min_stim_duration
                    fprintf('Stimulation %d is too short (%.1f seconds). Minimum required: %.1f seconds. Skipping.\n', ...
                        stim_idx, stim_duration_seconds, config.min_stim_duration);
                    continue;
                end
                
                fprintf('Stimulation duration: %.1f seconds (adequate for analysis)\n', stim_duration_seconds);
                
                % Create directory for this stimulation
                stim_dir = fullfile(config.results_dir, ['stim_' num2str(stim_idx)]);
                if ~exist(stim_dir, 'dir'), mkdir(stim_dir); end
                
                % Extract segments for this stimulation
                stim_segments = extractSegmentsForStimulation(segments, all_segment_names, stim_idx);
                this_stim_samples = stim_samples(stim_idx, :);
                
                % Update config for this stimulation
                stim_config = config;
                stim_config.results_dir = stim_dir;
                
                % Calculate time range for this stimulation
                time_range = calculateTimeRange(stim_segments);
                
                % Calculate and normalize ISPC
                [ispc_results, phases, filtered_data] = calculateISPC(EEG, stim_channel_idx, time_range, stim_config);
                [norm_ispc_results, z_scored_ispc] = normalizeISPC(ispc_results, stim_segments, globalISPC, globalStd);
                
                % Store results
                all_ispc_results{stim_idx} = ispc_results;
                all_norm_ispc_results{stim_idx} = norm_ispc_results;
                all_z_scored_ispc{stim_idx} = z_scored_ispc;
                all_stim_segments{stim_idx} = stim_segments;
                
                % Add to valid stimulation indices
                valid_stim_indices = [valid_stim_indices, stim_idx];
                
                % Visualize results
                visualizeISPCOverview(EEG, ispc_results, phases, filtered_data, stim_segments, this_stim_samples, stim_channel_idx, stim_config);
                visualizeNormalizedISPC(EEG, ispc_results, norm_ispc_results, z_scored_ispc, stim_segments, stim_config);
                
                % Save results
                save(fullfile(stim_dir, 'results.mat'), 'ispc_results', 'norm_ispc_results', 'z_scored_ispc', ...
                    'globalISPC', 'globalStd', 'stim_segments', 'this_stim_samples', '-v7.3');
                
                close all;
            end
            
            % Check if we have any valid stimulations
            if isempty(valid_stim_indices)
                fprintf('No valid stimulations found for subject %s, session %s. Skipping comparison visualizations.\n', ...
                    whichSubj, whichSess);
                continue;
            end
            
            % Filter out empty cells from results (for skipped stimulations)
            valid_ispc_results = all_ispc_results(valid_stim_indices);
            valid_norm_ispc_results = all_norm_ispc_results(valid_stim_indices);
            valid_stim_segments = all_stim_segments(valid_stim_indices);
            
            % Create comparison visualizations with valid stimulations only
            createStimulationComparison(EEG, valid_ispc_results, valid_stim_segments, config);
            
            % Create aggregate topoplots
            config.subject_id = whichSubj;
            [norm_final_topos, pct_change_topos] = createNormalizedAggregateTopoplots(EEG, valid_ispc_results, valid_norm_ispc_results, valid_stim_segments, config);
            
            % Save information about valid stimulations
            save(fullfile(config.results_dir, 'valid_stimulations.mat'), 'valid_stim_indices');
            
            fprintf('Analysis complete for subject %s, session %s\n', whichSubj, whichSess);
            
        catch err
            fprintf('\nERROR processing subject %s, session %s: %s\n', whichSubj, whichSess, err.message);
            disp(err.stack);
        end
    end
end

%% Run group-level analysis
fprintf('\n\n========== Group Analysis Pipeline ==========\n\n');

try
    % Interpolate subjects to standard electrode layout
    interpolateNormalizedSubjects(experiment_path, nights{1});
    
    % Run group-level analyses
    runGroupAnalysisWithInterpolatedData(experiment_path, nights{1});
    runProtocolGroupAnalysisWithInterpolation(experiment_path, nights{1}, 5);
    compareActiveAndShamTopoplots(experiment_path, nights{1});
    
    fprintf('Group-level analysis complete.\n');
catch err
    fprintf('Error in group analysis: %s\n', err.message);
    disp(err.stack);
end
