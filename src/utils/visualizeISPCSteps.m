function visualizeISPCSteps(EEG, timeRange, config)
    % VISUALIZEIISPCSTEPS Creates visualization of ISPC calculation steps for a specific time range
    %
    % This is a utility function that can be used to visualize the calculation steps
    % for a specific time range, independent of the main analysis pipeline.
    %
    % Inputs:
    %   EEG - EEGLAB EEG structure
    %   timeRange - [start_time, end_time] in seconds from the start of the recording
    %   config - Configuration structure with at least stim_freq field
    %
    % Example:
    %   visualizeISPCSteps(EEG, [30, 40], config);
    
    fprintf('Creating visualization of calculation steps for time range [%.2f-%.2f] seconds...\n', ...
        timeRange(1), timeRange(2));
    
    % Convert time range to samples
    sample_range = round(timeRange * EEG.srate);
    sample_range = max(1, min(sample_range, EEG.pnts));
    
    % Generate stimulus signal if needed
    if ~any(strcmp({EEG.chanlocs.labels}, 'STIM'))
        fprintf('Stimulus channel not found. Generating %.1f Hz stimulus signal...\n', config.stim_freq);
        time = (0:EEG.pnts-1) / EEG.srate;
        stimulus = sin(2 * pi * config.stim_freq * time);
        
        % Add the stimulus as a new channel to the EEG data (temporarily)
        EEG.data(end+1,:) = stimulus;
        EEG.nbchan = EEG.nbchan + 1;
        if isfield(EEG, 'chanlocs') && isstruct(EEG.chanlocs)
            EEG.chanlocs(end+1).labels = 'STIM';
        end
        stim_channel_idx = EEG.nbchan;
    else
        % Find the existing stimulus channel
        stim_channel_idx = find(strcmp({EEG.chanlocs.labels}, 'STIM'));
    end
    
    % Select a central channel for demonstration
    channel_labels = {EEG.chanlocs(1:end-1).labels};
    central_channels = {'Cz', 'CPz', 'FCz', 'Pz', 'Fz'};
    channel_idx = find(ismember(channel_labels, central_channels), 1);
    if isempty(channel_idx), channel_idx = 1; end
    
    % Calculate phases and ISPC specifically for the requested time range
    
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
    n_data = sample_range(2) - sample_range(1) + 1;
    n_convolution = n_wavelet + n_data - 1;
    
    % FFT of wavelet
    fft_wavelet = fft(wavelet, n_convolution);
    
    % Preallocate
    phases = zeros(EEG.nbchan, n_data);
    filtered_data = zeros(EEG.nbchan, n_data);
    
    % Process stimulus channel
    fft_stim = fft(EEG.data(stim_channel_idx, sample_range(1):sample_range(2)), n_convolution);
    convolution_result_fft = ifft(fft_wavelet .* fft_stim, n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet_size+1:end-half_wavelet_size);
    phases(stim_channel_idx, :) = angle(convolution_result_fft);
    filtered_data(stim_channel_idx, :) = real(convolution_result_fft);
    
    % Process the selected channel
    fft_data = fft(EEG.data(channel_idx, sample_range(1):sample_range(2)), n_convolution);
    convolution_result_fft = ifft(fft_wavelet .* fft_data, n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet_size+1:end-half_wavelet_size);
    phases(channel_idx, :) = angle(convolution_result_fft);
    filtered_data(channel_idx, :) = real(convolution_result_fft);
    
    % Calculate phase difference and ISPC
    phase_diff = phases(channel_idx, :) - phases(stim_channel_idx, :);
    z = exp(1i * phase_diff);
    ispc = abs(mean(z));
    
    % Create time vector for plotting
    abs_time = timeRange(1) + (0:n_data-1) / EEG.srate;
    rel_time = abs_time - abs_time(1);
    
    % Create figure
    fig = figure('Name', 'ISPC Calculation Steps', 'Position', [100, 100, 1200, 800], 'visible', 'off');
    
    % Panel 1: Raw, filtered EEG and stimulus signals overlaid
    subplot(4, 1, 1);
    % Normalize signals for better visualization
    raw_eeg_norm = EEG.data(channel_idx, sample_range(1):sample_range(2)) / max(abs(EEG.data(channel_idx, sample_range(1):sample_range(2))));
    filtered_eeg_norm = filtered_data(channel_idx, :) / max(abs(filtered_data(channel_idx, :)));
    stim_norm = EEG.data(stim_channel_idx, sample_range(1):sample_range(2)) / max(abs(EEG.data(stim_channel_idx, sample_range(1):sample_range(2))));
    
    % Plot all signals
    plot(rel_time, raw_eeg_norm, 'k', 'LineWidth', 1);
    hold on;
    plot(rel_time, filtered_eeg_norm, 'b', 'LineWidth', 1.5);
    plot(rel_time, stim_norm, 'r', 'LineWidth', 1.5);
    title('Overlaid Signals (Normalized)');
    xlabel('Time (s)'); ylabel('Amplitude');
    legend('Raw EEG', ['Filtered EEG (' num2str(config.stim_freq) ' Hz)'], 'Stimulus Signal');
    grid on;
    
    % Panel 2: Phase angles
    subplot(4, 1, 2);
    plot(rel_time, phases(channel_idx, :), 'b', 'LineWidth', 1.5);
    hold on;
    plot(rel_time, phases(stim_channel_idx, :), 'r', 'LineWidth', 1.5);
    title('Phase Angles');
    xlabel('Time (s)'); ylabel('Phase (rad)');
    legend('EEG', 'Stimulus');
    ylim([-pi pi]); % Setting y-limits to [-π, π]
    grid on;
    
    % Panel 3: Phase difference
    subplot(4, 1, 3);
    plot(rel_time, phase_diff, 'k', 'LineWidth', 1.5);
    title('Phase Difference');
    xlabel('Time (s)'); ylabel('Phase Difference (rad)');
    ylim([-pi pi]); % Setting y-limits to [-π, π]
    grid on;
    
    % Panel 4: Phase difference on complex plane
    subplot(4, 1, 4);
    % Convert to polar coordinates for plotting
    scatter(real(z), imag(z), 20, 'filled', 'MarkerFaceAlpha', 0.3);
    hold on;
    % Add mean vector
    mean_z = mean(z);
    quiver(0, 0, real(mean_z), imag(mean_z), 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    % Add unit circle
    theta = linspace(0, 2*pi, 100);
    plot(cos(theta), sin(theta), 'k--');
    axis equal;
    xlim([-1.1 1.1]); ylim([-1.1 1.1]);
    grid on;
    title(['Phase Difference on Complex Plane (ISPC = ' num2str(ispc, '%.3f') ')']);
    xlabel('Real'); ylabel('Imaginary');
    
    sgtitle(['ISPC Calculation for Channel ' EEG.chanlocs(channel_idx).labels ...
        ' (' num2str(timeRange(1), '%.1f') '-' num2str(timeRange(2), '%.1f') ' s)']);
    
    % If a save directory is provided in config, save the figure
    if isfield(config, 'results_dir') && ~isempty(config.results_dir)
        if ~exist(config.results_dir, 'dir')
            mkdir(config.results_dir);
        end
        
        % Create a filename based on the time range
        filename = sprintf('ispc_steps_%.1f_%.1f.png', timeRange(1), timeRange(2));
        save_path = fullfile(config.results_dir, filename);
        saveas(gcf, save_path);
        fprintf('Saved figure to: %s\n', save_path);
    end
end
