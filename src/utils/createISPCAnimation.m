function createISPCAnimation(EEG, timeRange, outputFile, config)
    % CREATEISPCANIMATION Creates an animation of phase synchronization for a specific time range
    %
    % This is a utility function that can be used to create an animation
    % independently of the main analysis pipeline.
    %
    % Inputs:
    %   EEG - EEGLAB EEG structure
    %   timeRange - [start_time, end_time] in seconds from the start of the recording
    %   outputFile - Full path to the output GIF file
    %   config - Configuration structure with at least stim_freq field
    %
    % Example:
    %   createISPCAnimation(EEG, [30, 40], 'my_animation.gif', config);
    
    fprintf('Creating phase animation for time range [%.2f-%.2f] seconds...\n', ...
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
    
    % Calculate phases and filtered data specifically for the requested time range
    n_samples = sample_range(2) - sample_range(1) + 1;
    
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
    n_convolution = n_wavelet + n_samples - 1;
    
    % FFT of wavelet
    fft_wavelet = fft(wavelet, n_convolution);
    
    % Preallocate
    phases = zeros(2, n_samples); % Only for selected channel and stimulus
    filtered_data = zeros(2, n_samples);
    
    % Process stimulus channel
    fft_stim = fft(EEG.data(stim_channel_idx, sample_range(1):sample_range(2)), n_convolution);
    convolution_result_fft = ifft(fft_wavelet .* fft_stim, n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet_size+1:end-half_wavelet_size);
    phases(2, :) = angle(convolution_result_fft);
    filtered_data(2, :) = real(convolution_result_fft);
    
    % Process the selected channel
    fft_data = fft(EEG.data(channel_idx, sample_range(1):sample_range(2)), n_convolution);
    convolution_result_fft = ifft(fft_wavelet .* fft_data, n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet_size+1:end-half_wavelet_size);
    phases(1, :) = angle(convolution_result_fft);
    filtered_data(1, :) = real(convolution_result_fft);
    
    % Create absolute and relative time vectors
    abs_time = timeRange(1) + (0:n_samples-1) / EEG.srate;
    rel_time = abs_time - abs_time(1);
    
    % Animation parameters
    window_size = round(EEG.srate * 2); % 2-second window for visualization
    step_size = round(0.1 * EEG.srate); % advance by 0.1 seconds each frame (smoother)
    frame_delay = 0.1; % delay between frames in seconds
    display_window = 5; % Display 5 seconds of data at a time
    
    % Create figure for animation
    fig = figure('Name', 'Phase Synchronization Animation', 'Position', [100, 100, 1200, 600]);
    
    % Setup signals subplot (top-left)
    subplot(2, 2, 1);
    % Initialize with empty data - will grow as animation proceeds
    h_eeg = plot(0, 0, 'b-', 'LineWidth', 1.5);
    hold on;
    h_stim = plot(0, 0, 'r-', 'LineWidth', 1.5);
    ylim([-1.2 1.2]);
    xlim([0 display_window]);
    title('Filtered Signals');
    xlabel('Time (s)');
    ylabel('Amplitude (Normalized)');
    legend('EEG', 'Stimulus');
    grid on;
    
    % Setup phase angle subplot (top-right)
    subplot(2, 2, 2);
    h_phase_eeg = plot(0, 0, 'b-', 'LineWidth', 1.5);
    hold on;
    h_phase_stim = plot(0, 0, 'r-', 'LineWidth', 1.5);
    h_phase_diff = plot(0, 0, 'k-', 'LineWidth', 1);
    ylim([-pi pi]);
    xlim([0 display_window]);
    title('Phase Angles');
    xlabel('Time (s)');
    ylabel('Phase (rad)');
    legend('EEG Phase', 'Stimulus Phase', 'Phase Difference');
    grid on;
    
    % Setup polar subplot (bottom-left)
    subplot(2, 2, 3);
    % Start with placeholder values
    initial_phase_diff = phases(1, 1) - phases(2, 1);
    z_init = exp(1i * initial_phase_diff);
    h_polar = polarscatter(angle(z_init), 1, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.3);
    hold on;
    h_mean_vector = polarplot([0 angle(z_init)], [0 abs(z_init)], 'r-', 'LineWidth', 2); % mean vector
    title('Phase Differences on Complex Plane');
    rlim([0 1.1]);
    
    % Setup ISPC plot subplot (bottom-right)
    subplot(2, 2, 4);
    h_ispc = plot(0, 0, 'b-', 'LineWidth', 2);
    hold on;
    h_current = plot([0 0], [0 1], 'r--', 'LineWidth', 1);
    ylim([0 1]);
    xlim([0 display_window]);
    xlabel('Time (s)');
    ylabel('ISPC');
    title('ISPC Over Time');
    grid on;
    
    % Initialize arrays to store data for the animation
    eeg_times = [];
    eeg_values = [];
    stim_values = [];
    phase_eeg_values = [];
    phase_stim_values = [];
    phase_diff_values = [];
    ispc_values = [];
    
    % Prepare the GIF file
    %gif_filename = fullfile(config.results_dir, 'phase_animation.gif');
    frame_count = 0;
    
    try
        % Loop through time points, making sure we don't go beyond the end with the window
        for t_idx = 1:step_size:n_samples-window_size
            % Get actual sample indices for this window
            window_end = t_idx + window_size - 1;
            
            % Calculate current time since animation start
            current_rel_time = rel_time(t_idx);
            
            % Get phase differences for this window
            phase_diff = phases(1, t_idx:window_end) - phases(2, t_idx:window_end);
            
            % Convert to complex numbers on the unit circle
            z = exp(1i * phase_diff);
            mean_z = mean(z);
            current_ispc = abs(mean_z);
            
            % Append new data to arrays
            eeg_times = [eeg_times, current_rel_time];
            
            % Normalize signal amplitudes
            eeg_val = filtered_data(1, t_idx) / max(abs(filtered_data(1, :)));
            stim_val = filtered_data(2, t_idx) / max(abs(filtered_data(2, :)));
            
            eeg_values = [eeg_values, eeg_val];
            stim_values = [stim_values, stim_val];
            phase_eeg_values = [phase_eeg_values, phases(1, t_idx)];
            phase_stim_values = [phase_stim_values, phases(2, t_idx)];
            phase_diff_values = [phase_diff_values, phase_diff(1)];
            ispc_values = [ispc_values, current_ispc];
            
            % Determine the xlim range to show a moving window
            if current_rel_time > display_window
                x_min = current_rel_time - display_window;
                x_max = current_rel_time;
            else
                x_min = 0;
                x_max = display_window;
            end
            
            % Update signal plot
            subplot(2, 2, 1);
            set(h_eeg, 'XData', eeg_times, 'YData', eeg_values);
            set(h_stim, 'XData', eeg_times, 'YData', stim_values);
            xlim([x_min x_max]);
            
            % Update phase angle plot
            subplot(2, 2, 2);
            set(h_phase_eeg, 'XData', eeg_times, 'YData', phase_eeg_values);
            set(h_phase_stim, 'XData', eeg_times, 'YData', phase_stim_values);
            set(h_phase_diff, 'XData', eeg_times, 'YData', phase_diff_values);
            xlim([x_min x_max]);
            
            % Update polar plot
            subplot(2, 2, 3);
            set(h_polar, 'ThetaData', angle(z), 'RData', ones(size(z)));
            set(h_mean_vector, 'ThetaData', [0 angle(mean_z)], 'RData', [0 abs(mean_z)]);
            title(sprintf('Phase Differences (ISPC = %.3f)', current_ispc));
            
            % Update ISPC plot
            subplot(2, 2, 4);
            set(h_ispc, 'XData', eeg_times, 'YData', ispc_values);
            set(h_current, 'XData', [current_rel_time current_rel_time]);
            xlim([x_min x_max]);
            
            % Add overall title
            sgtitle(['ISPC Animation: Channel ' EEG.chanlocs(channel_idx).labels ' vs. Stimulus'], 'FontSize', 12);
            
            % Update figure
            drawnow;
            
            % Capture frame for GIF
            frame = getframe(fig);
            im = frame2im(frame);
            [imind, cm] = rgb2ind(im, 256);
            
            % Write to GIF file
            if frame_count == 0
                imwrite(imind, cm, outputFile, 'gif', 'Loopcount', inf, 'DelayTime', frame_delay);
            else
                imwrite(imind, cm, outputFile, 'gif', 'WriteMode', 'append', 'DelayTime', frame_delay);
            end
            
            frame_count = frame_count + 1;
        end
        
        fprintf('Animation saved to %s\n', outputFile);
    catch err
        warning('Error creating animation: %s\n', err.message);
        
        % Try a simplified approach for more compatibility
        try
            fprintf('Attempting to create a simplified animation...\n');
            % Close and recreate figure
            close(fig);
            
            % Create a simpler figure
            fig = figure('Name', 'Phase Synchronization', 'Position', [100, 100, 800, 600]);
            
            % Create a static image for each frame
            % Prepare the GIF file with simpler approach
            [filepath, name, ~] = fileparts(outputFile);
            simplified_gif = fullfile(filepath, [name '_simple.gif']);
            frame_count = 0;
            
            % Use fewer frames for the simplified version
            for t_idx = 1:step_size*2:n_samples-window_size
                % Calculate current time since animation start
                current_rel_time = rel_time(t_idx);
                
                % Get phase differences for this window
                phase_diff = phases(1, t_idx:min(t_idx+window_size-1, n_samples)) - ...
                             phases(2, t_idx:min(t_idx+window_size-1, n_samples));
                
                % Convert to complex numbers on the unit circle
                z = exp(1i * phase_diff);
                mean_z = mean(z);
                current_ispc = abs(mean_z);
                
                % Create a simple figure for this frame
                clf;
                
                % Plot phase points on complex plane
                polarscatter(angle(z), ones(size(z)), 10, 'k', 'filled', 'MarkerFaceAlpha', 0.3);
                hold on;
                polarplot([0 angle(mean_z)], [0 abs(mean_z)], 'r-', 'LineWidth', 2);
                title(sprintf('Phase Differences at %.1f s (ISPC = %.3f)', current_rel_time, current_ispc));
                
                % Update figure
                drawnow;
                
                % Capture frame for GIF
                frame = getframe(fig);
                im = frame2im(frame);
                [imind, cm] = rgb2ind(im, 256);
                
                % Write to GIF file
                if frame_count == 0
                    imwrite(imind, cm, simplified_gif, 'gif', 'Loopcount', inf, 'DelayTime', frame_delay*2);
                else
                    imwrite(imind, cm, simplified_gif, 'gif', 'WriteMode', 'append', 'DelayTime', frame_delay*2);
                end
                
                frame_count = frame_count + 1;
            end
            
            fprintf('Simplified animation saved to %s\n', simplified_gif);
        catch err2
            warning('Failed to create simplified animation: %s\n', err2.message);
            fprintf('Please check if your MATLAB version supports GIF creation.\n');
        end
    end
    
    % Close figure at the end
    if exist('fig', 'var') && ishandle(fig)
        close(fig);
    end
end
