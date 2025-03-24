function EEG = loadEEGData(eeg_filename, target_srate)
    % Load EEG data and resample if needed
    
    % Load EEG data
    EEG = pop_loadset(eeg_filename);
    
    fprintf('Loaded EEG data: %d channels, %d time points, %.2f Hz sampling rate\n', ...
        EEG.nbchan, EEG.pnts, EEG.srate);
    
    % Resample data if needed
    if EEG.srate > target_srate
        fprintf('Resampling data from %.2f Hz to %d Hz for faster processing...\n', EEG.srate, target_srate);
        EEG = pop_resample(EEG, target_srate);
        fprintf('After resampling: %d channels, %d time points, %.2f Hz sampling rate\n', ...
            EEG.nbchan, EEG.pnts, EEG.srate);
    end
end
