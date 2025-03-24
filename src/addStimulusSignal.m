function [EEG, stim_channel_idx] = addStimulusSignal(EEG, stim_freq)
    % Generate 1Hz sine wave stimulus for the entire duration
    fprintf('Generating %.1f Hz stimulus signal...\n', stim_freq);
    time = (0:EEG.pnts-1) / EEG.srate;
    stimulus = sin(2 * pi * stim_freq * time);
    
    % Add the stimulus as a new channel to the EEG data
    EEG.data(end+1,:) = stimulus;
    EEG.nbchan = EEG.nbchan + 1;
    if isfield(EEG, 'chanlocs') && isstruct(EEG.chanlocs)
        EEG.chanlocs(end+1).labels = 'STIM';
    end
    stim_channel_idx = EEG.nbchan;
    
    fprintf('Added stimulus channel at index %d\n', stim_channel_idx);
end
