function visualizeWindowSliding(EEG, timeRange, config)
    % VISUALIZEWINDOWSLIDING Creates visualization of how the sliding window
    % affects ISPC calculation over time
    %
    % This utility function demonstrates how the window size parameter affects
    % the ISPC calculation by visualizing:
    % 1. The raw signal with moving window
    % 2. The phase differences within each window on the complex plane
    % 3. How the ISPC vector (mean of complex values) changes with window position
    %
    % Inputs:
    %   EEG - EEGLAB EEG structure
    %   timeRange - [start_time, end_time] in seconds from the start of the recording
    %   config - Configuration structure with at least stim_freq and window_size fields
    %
    % Example:
    %   visualizeWindowSliding(EEG, [30, 35], config);
    
    fprintf('Creating visualization of sliding window ISPC calculation for time range [%.2f-%.2f] seconds...\n', ...
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
    
    % Select a channel for demonstration (use a central channel if available)
    channel_labels = {EEG.chanlocs(1:end-1).labels};
    central_channels = {'Cz', 'CPz', 'FCz', 'Pz', 'Fz'};
    channel_idx = find(ismember(channel_labels, central_channels), 1);
    if isempty(channel_idx), channel_idx = 1; end
    
    % Extract data for the specified time range
    start_idx = sample_range(1);
    end_idx = sample_range(2);
    
    % Get time vector for plotting
    abs_time = (start_idx:end_idx) / EEG.srate;
    rel_time = abs_time - abs_time(1);
    
    % Calculate phases using wavelet convolution (similar to calculateISPC.m)
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
    n_data = end_idx - start_idx + 1;
    n_convolution = n_wavelet + n_data - 1;
    
    % FFT of wavelet
    fft_wavelet = fft(wavelet, n_convolution);
    
    % Process stimulus channel
    fft_stim = fft(EEG.data(stim_channel_idx, start_idx:end_idx), n_convolution);
    convolution_result_fft = ifft(fft_wavelet .* fft_stim, n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet_size+1:end-half_wavelet_size);
    stim_phase = angle(convolution_result_fft);
    stim_filtered = real(convolution_result_fft);
    
    % Process the selected channel
    fft_data = fft(EEG.data(channel_idx, start_idx:end_idx), n_convolution);
    convolution_result_fft = ifft(fft_wavelet .* fft_data, n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet_size+1:end-half_wavelet_size);
    eeg_phase = angle(convolution_result_fft);
    eeg_filtered = real(convolution_result_fft);
    
    % Calculate phase difference
    phase_diff = eeg_phase - stim_phase;
    
    % Convert to complex representation on unit circle
    complex_phase = exp(1i * phase_diff);
    
    % Calculate ISPC for each time point using sliding window
    window_samples = round(config.window_size * EEG.srate);
    half_window = floor(window_samples / 2);
    ispc_values = zeros(1, n_data);
    
    for t = 1:n_data
        % Define window boundaries (ensuring we stay within data limits)
        win_start = max(1, t - half_window);
        win_end = min(n_data, t + half_window);
        
        % Calculate ISPC for this window
        ispc_values(t) = abs(mean(complex_phase(win_start:win_end)));
    end
    
    % Create a figure to visualize the sliding window
    % We'll create a series of subplots for different window positions
    
    % Number of snapshots to show (evenly spaced through the data)
    num_snapshots = 5;
    snapshot_indices = round(linspace(half_window+1, n_data-half_window, num_snapshots));
    
    % Create figure
    fig = figure('Name', 'Sliding Window ISPC Visualization', 'Position', [100, 100, 1200, 800]);
    
    % Create subplots for each snapshot
    for snap_idx = 1:num_snapshots
        t = snapshot_indices(snap_idx); % Current time point
        
        % Define window boundaries
        win_start = max(1, t - half_window);
        win_end = min(n_data, t + half_window);
        
        % Current window phase differences
        window_complex = complex_phase(win_start:win_end);
        current_ispc = abs(mean(window_complex));
        mean_vector = mean(window_complex);
        
        % Create subplots for this snapshot
        
        % 1. Signal with window highlight
        subplot(num_snapshots, 3, (snap_idx-1)*3 + 1);
        
        % Plot filtered signals
        plot(rel_time, eeg_filtered, 'b-', 'LineWidth', 1);
        hold on;
        plot(rel_time, stim_filtered, 'r-', 'LineWidth', 1);
        
        % Highlight the current window
        x_window = [rel_time(win_start), rel_time(win_end), rel_time(win_end), rel_time(win_start)];
        y_limits = ylim();
        y_window = [y_limits(1), y_limits(1), y_limits(2), y_limits(2)];
        patch(x_window, y_window, 'yellow', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        
        % Add vertical line for current time point - FIXED LINE
        plot([rel_time(t), rel_time(t)], y_limits, 'k--', 'LineWidth', 1.5);
        
        title(sprintf('Window at %.1f s', abs_time(t)));
        xlabel('Time (s)');
        ylabel('Amplitude');
        legend('EEG', 'Stimulus', 'Window', 'Current Time');
        
        % 2. Phase differences on complex plane for this window
        subplot(num_snapshots, 3, (snap_idx-1)*3 + 2);
        
        % Plot unit circle
        theta = linspace(0, 2*pi, 100);
        plot(cos(theta), sin(theta), 'k--');
        hold on;
        axis equal;
        xlim([-1.1 1.1]); ylim([-1.1 1.1]);
        
        % Plot phase points
        scatter(real(window_complex), imag(window_complex), 20, 'filled', 'MarkerFaceAlpha', 0.5);
        
        % Plot mean vector
        quiver(0, 0, real(mean_vector), imag(mean_vector), 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
        
        title(sprintf('Phase Differences (ISPC = %.3f)', current_ispc));
        xlabel('Real'); ylabel('Imaginary');
        
        % 3. ISPC time course with current point
        subplot(num_snapshots, 3, (snap_idx-1)*3 + 3);
        
        % Plot ISPC time course
        plot(rel_time, ispc_values, 'b-', 'LineWidth', 1.5);
        hold on;
        
        % Highlight current point
        scatter(rel_time(t), ispc_values(t), 50, 'r', 'filled');
        
        % Highlight the window
        x_window = [rel_time(win_start), rel_time(win_end), rel_time(win_end), rel_time(win_start)];
        y_limits = [0, 1];
        y_window = [y_limits(1), y_limits(1), y_limits(2), y_limits(2)];
        patch(x_window, y_window, 'yellow', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        
        title('ISPC Time Course');
        xlabel('Time (s)');
        ylabel('ISPC');
        ylim([0 1]);
    end
    
    % Overall title
    sgtitle(sprintf('Sliding Window ISPC Calculation (Window Size: %.1f s)', config.window_size), 'FontSize', 14);
    
    % Create animation of the sliding window
    % This will show how the window moves over time and how the ISPC calculation changes
    
    % Create separate figure for animation
    anim_fig = figure('Name', 'Sliding Window Animation', 'Position', [100, 100, 1200, 400]);
    
    % Create subplots for the animation
    subplot(1, 3, 1); % Signal with window
    sig_plot = plot(rel_time, eeg_filtered, 'b-', rel_time, stim_filtered, 'r-');
    hold on;
    window_patch = patch([0 0 0 0], [0 0 0 0], 'yellow', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    time_line_y = ylim;
    time_line = plot([0 0], time_line_y, 'k--', 'LineWidth', 1.5);
    title('Filtered Signals');
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend('EEG', 'Stimulus');
    
    subplot(1, 3, 2); % Complex plane
    % Plot unit circle
    theta = linspace(0, 2*pi, 100);
    plot(cos(theta), sin(theta), 'k--');
    hold on;
    axis equal;
    xlim([-1.1 1.1]); ylim([-1.1 1.1]);
    phase_scatter = scatter([], [], 20, 'filled', 'MarkerFaceAlpha', 0.5);
    mean_arrow = quiver(0, 0, 0, 0, 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    phase_title = title('Phase Differences');
    xlabel('Real'); ylabel('Imaginary');
    
    subplot(1, 3, 3); % ISPC time course
    ispc_plot = plot(rel_time, ispc_values, 'b-', 'LineWidth', 1.5);
    hold on;
    current_point = scatter(0, 0, 50, 'r', 'filled');
    ispc_window_patch = patch([0 0 0 0], [0 0 0 0], 'yellow', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    title('ISPC Time Course');
    xlabel('Time (s)');
    ylabel('ISPC');
    ylim([0 1]);
    
    % Overall title
    sgtitle(sprintf('Sliding Window ISPC Calculation (Window Size: %.1f s)', config.window_size), 'FontSize', 14);
    
    % Animation parameters
    step = max(1, round(n_data / 100)); % Skip points for smoother animation
    
    % Prepare the GIF file
    if isfield(config, 'results_dir') && ~isempty(config.results_dir)
        gif_filename = fullfile(config.results_dir, 'sliding_window_animation.gif');
    else
        % Use current directory as fallback
        gif_filename = 'sliding_window_animation.gif';
    end
    
    % Initialize GIF
    frame_delay = 0.1; % seconds between frames
    
    % Create animation frames
    for t = (half_window+1):step:(n_data-half_window)
        % Define window boundaries
        win_start = max(1, t - half_window);
        win_end = min(n_data, t + half_window);
        
        % Current window phase differences
        window_complex = complex_phase(win_start:win_end);
        current_ispc = abs(mean(window_complex));
        mean_vector = mean(window_complex);
        
        % Update signal plot with window
        subplot(1, 3, 1);
        x_window = [rel_time(win_start), rel_time(win_end), rel_time(win_end), rel_time(win_start)];
        y_lim = ylim();
        y_window = [y_lim(1), y_lim(1), y_lim(2), y_lim(2)];
        set(window_patch, 'XData', x_window, 'YData', y_window);
        set(time_line, 'XData', [rel_time(t), rel_time(t)], 'YData', y_lim);
        
        % Update complex plane
        subplot(1, 3, 2);
        set(phase_scatter, 'XData', real(window_complex), 'YData', imag(window_complex));
        set(mean_arrow, 'UData', real(mean_vector), 'VData', imag(mean_vector));
        set(phase_title, 'String', sprintf('Phase Differences (ISPC = %.3f)', current_ispc));
        
        % Update ISPC plot
        subplot(1, 3, 3);
        set(current_point, 'XData', rel_time(t), 'YData', ispc_values(t));
        x_window = [rel_time(win_start), rel_time(win_end), rel_time(win_end), rel_time(win_start)];
        y_window = [0, 0, 1, 1];
        set(ispc_window_patch, 'XData', x_window, 'YData', y_window);
        
        % Force update
        drawnow;
        
        % Capture frame for GIF
        frame = getframe(anim_fig);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        
        % Write to GIF file
        if t == (half_window+1)
            imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', frame_delay);
        else
            imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', frame_delay);
        end
    end
    
    fprintf('Sliding window visualization complete.\n');
    fprintf('Saved animation to: %s\n', gif_filename);
    
    % Save static figure
    if isfield(config, 'results_dir') && ~isempty(config.results_dir)
        saveas(fig, fullfile(config.results_dir, 'sliding_window_snapshots.png'));
        fprintf('Saved snapshot figure to: %s\n', fullfile(config.results_dir, 'sliding_window_snapshots.png'));
    end
    
    % Close animation figure
    close(anim_fig);
end
