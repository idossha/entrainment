function [segments, stim_samples] = findEvents(EEG, config)
    % Find all instances of stim start/end events and define segments
    
    % Get marker settings
    start_marker = config.protocol.start_marker;
    end_marker = config.protocol.end_marker;
    segment_duration = config.protocol.segment_duration;
    
    % Find event indices for start and end markers
    start_marker_events = strcmpi({EEG.event.type}, start_marker);
    end_marker_events = strcmpi({EEG.event.type}, end_marker);
    
    if ~any(start_marker_events) || ~any(end_marker_events)
        error('Start or end markers not found: %s/%s', start_marker, end_marker);
    end
    
    % Get all occurrences of start and end events
    start_indices = find(start_marker_events);
    end_indices = find(end_marker_events);
    
    fprintf('Found %d %s events and %d %s events\n', ...
        length(start_indices), start_marker, length(end_indices), end_marker);
    
    % Match starts and ends - simple pairing approach
    stim_pairs = [];
    for i = 1:length(start_indices)
        start_idx = start_indices(i);
        start_sample = round(EEG.event(start_idx).latency);
        
        % Find the next end marker after this start
        valid_ends = find(end_indices > start_idx, 1, 'first');
        
        if ~isempty(valid_ends)
            end_idx = end_indices(valid_ends);
            end_sample = round(EEG.event(end_idx).latency);
            
            % Only add if the stimulation duration makes sense (should be positive and not too short)
            if end_sample > start_sample && (end_sample - start_sample) > EEG.srate
                stim_pairs(end+1, :) = [start_idx, end_idx, start_sample, end_sample];
            end
        end
    end
    
    if isempty(stim_pairs)
        error('No valid start-end event pairs found');
    end
    
    % Number of valid stimulation pairs
    num_stims = size(stim_pairs, 1);
    fprintf('Found %d valid stimulation periods\n', num_stims);
    
    % Initialize outputs
    segments = struct();
    stim_samples = zeros(num_stims, 2);
    
    % Calculate time segments for each stimulation
    srate = EEG.srate;
    segment_samples = round(segment_duration * srate);
    
    for i = 1:num_stims
        % Get start and end samples
        stim_start_sample = stim_pairs(i, 3);
        stim_end_sample = stim_pairs(i, 4);
        stim_samples(i, :) = [stim_start_sample, stim_end_sample];
        
        % Create segment name with index if multiple stimulations
        suffix = '';
        if num_stims > 1
            suffix = ['_' num2str(i)];
        end
        
        % Define the four segments for this stimulation
        % Pre-stim segment: 45 seconds before stimulation
        segments.(['pre_stim' suffix]) = [max(1, stim_start_sample - segment_samples), stim_start_sample - 1];
        
        % Early-stim segment: first 45 seconds of stimulation (or all if shorter)
        early_stim_end = min(stim_start_sample + segment_samples - 1, stim_end_sample);
        segments.(['early_stim' suffix]) = [stim_start_sample, early_stim_end];
        
        % Late-stim segment: last 45 seconds of stimulation (or all if shorter)
        late_stim_start = max(stim_start_sample, stim_end_sample - segment_samples + 1);
        segments.(['late_stim' suffix]) = [late_stim_start, stim_end_sample];
        
        % Post-stim segment: 45 seconds after stimulation
        segments.(['post_stim' suffix]) = [stim_end_sample + 1, min(stim_end_sample + segment_samples, EEG.pnts)];
        
        fprintf('Stimulation %d: starts at %.2f s and ends at %.2f s (duration: %.2f s)\n', ...
            i, stim_start_sample/srate, stim_end_sample/srate, (stim_end_sample-stim_start_sample)/srate);
    end
end
