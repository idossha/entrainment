function [ispc_results, phases, filtered_data] = calculateISPC(EEG, stim_channel_idx, time_range, config)
    % Calculates ISPC between each channel and the stimulus for the specified time range
    % time_range: [start_sample, end_sample] - the range of samples to analyze
    
    if nargin < 3 || isempty(time_range)
        % If no time range specified, process entire dataset
        time_range = [1, EEG.pnts];
        fprintf('Calculating ISPC for entire dataset (%d samples)...\n', EEG.pnts);
    else
        fprintf('Calculating ISPC for time range [%d-%d] (%d samples)...\n', ...
            time_range(1), time_range(2), time_range(2)-time_range(1)+1);
    end
    
    % Extract the relevant data segment
    start_idx = time_range(1);
    end_idx = time_range(2);
    
    % Number of samples in the selected range
    n_samples = end_idx - start_idx + 1;
    
    % Initialize output variables for full dataset (we'll fill only the relevant part)
    ispc_results = zeros(EEG.nbchan-1, EEG.pnts); % ISPC time courses
    phases = zeros(EEG.nbchan, EEG.pnts); % Phase angles
    filtered_data = zeros(EEG.nbchan, EEG.pnts); % Filtered signals
    
    % Wavelet parameters
    center_freq = config.stim_freq; % Analyze at the stimulation frequency
    cycles = 8; % Number of cycles in the wavelet
    wavelet_time = -2:1/EEG.srate:2; % Time for wavelet (4 seconds)
    s = cycles/(2*pi*center_freq); % Gaussian width
    
    % Create Morlet wavelet
    wavelet = exp(2*1i*pi*center_freq.*wavelet_time) .* exp(-wavelet_time.^2./(2*s^2))/center_freq;
    half_wavelet_size = (length(wavelet_time)-1)/2;
    
    % FFT parameters
    n_wavelet = length(wavelet_time);
    n_data = n_samples; % Only use the specified range
    n_convolution = n_wavelet + n_data - 1;
    
    % FFT of wavelet
    fft_wavelet = fft(wavelet, n_convolution);
    
    % First get the stimulus phase (last channel)
    fft_stim = fft(EEG.data(stim_channel_idx, start_idx:end_idx), n_convolution);
    convolution_result_fft = ifft(fft_wavelet .* fft_stim, n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet_size+1:end-half_wavelet_size);
    phases(stim_channel_idx, start_idx:end_idx) = angle(convolution_result_fft);
    filtered_data(stim_channel_idx, start_idx:end_idx) = real(convolution_result_fft);
    
    % Sliding window parameters for ISPC calculation
    window_size = round(EEG.srate * config.window_size); % Default: 2-second window
    half_window = floor(window_size/2);
    
    % Calculate ISPC for each channel
    for chan_idx = 1:EEG.nbchan-1
        % FFT of channel data (only for the specified range)
        fft_data = fft(EEG.data(chan_idx, start_idx:end_idx), n_convolution);
        
        % Convolution with wavelet
        convolution_result_fft = ifft(fft_wavelet .* fft_data, n_convolution);
        convolution_result_fft = convolution_result_fft(half_wavelet_size+1:end-half_wavelet_size);
        
        % Extract phase and filtered signal
        phases(chan_idx, start_idx:end_idx) = angle(convolution_result_fft);
        filtered_data(chan_idx, start_idx:end_idx) = real(convolution_result_fft);
        
        % Calculate phase difference with stimulus
        phase_diff = phases(chan_idx, start_idx:end_idx) - phases(stim_channel_idx, start_idx:end_idx);
        
        % Calculate ISPC using a sliding window
        for t_offset = 1:n_samples
            t = start_idx + t_offset - 1;
            
            % Define window boundaries
            win_start_offset = max(1, t_offset - half_window);
            win_end_offset = min(n_samples, t_offset + half_window);
            
            % Calculate ISPC for this window
            ispc_results(chan_idx, t) = abs(mean(exp(1i * phase_diff(win_start_offset:win_end_offset))));
        end
        
        if mod(chan_idx, 20) == 0
            fprintf('Processed %d of %d channels\n', chan_idx, EEG.nbchan-1);
        end
    end
    
    fprintf('ISPC calculation complete for time range [%d-%d].\n', start_idx, end_idx);
end
