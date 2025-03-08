% MATLAB Script for EEG Entrainment Analysis
% This script generates 4 minutes of simulated NREM sleep EEG data with 64 channels
% plus a 65th channel with a 1Hz sine wave for entrainment analysis
% First 2 minutes: Random EEG with minimal entrainment
% Second 2 minutes: Clear entrainment to the 1Hz stimulus
%% Initialize EEGLAB
clear; close all; clc;
% Make sure EEGLAB is in your MATLAB path, otherwise you'll need to add it:
addpath('/Users/idohaber/Documents/MATLAB/eeglab2024.0/');
eeglab;
% load data
[EEG] = pop_loadset('../data/nrem_sleep_with_entrainment_contrast.set');

%% Parameters for simulated EEG data
fs = EEG.srate;                   % Sampling frequency (Hz)
total_duration = EEG.pnts/fs;
segment_duration = total_duration/2;    % 2 minutes in seconds (for each segment)
n_channels = 64;            % Number of main EEG channels
sine_freq = 1;              % Frequency of the entrainment sine wave (Hz)

% Create separate time vectors for each segment
t1 = 0:1/fs:segment_duration-1/fs;   % First 2 minutes
t2 = 0:1/fs:segment_duration-1/fs;   % Second 2 minutes
t_full = [t1, t1(end) + t2];         % Full 4 minutes

n_samples_segment = length(t1);      % Number of samples per 2-min segment
n_samples_total = length(t_full);    % Total number of samples


%% Analysis 2: Improved Inter-trial Phase Coherence (ITPC) with Simple Visualization
% We'll calculate coherence between each EEG channel and the reference SI-ENV channel
% Create epochs of 2 seconds each


epoch_length = 2; % in seconds
samples_per_epoch = epoch_length * fs;

% Process each segment separately
% First 2 minutes
n_epochs_segment1 = floor(n_samples_segment / samples_per_epoch);
epoched_data1 = zeros(EEG.nbchan, samples_per_epoch, n_epochs_segment1);
for ep = 1:n_epochs_segment1
    start_idx = (ep-1) * samples_per_epoch + 1;
    end_idx = ep * samples_per_epoch;
    epoched_data1(:, :, ep) = EEG.data(:, start_idx:end_idx);
end

% Second 2 minutes
n_epochs_segment2 = floor(n_samples_segment / samples_per_epoch);
epoched_data2 = zeros(EEG.nbchan, samples_per_epoch, n_epochs_segment2);
for ep = 1:n_epochs_segment2
    start_idx = n_samples_segment + (ep-1) * samples_per_epoch + 1;
    end_idx = n_samples_segment + ep * samples_per_epoch;
    epoched_data2(:, :, ep) = EEG.data(:, start_idx:end_idx);
end

% Create epoched EEGLAB datasets for each segment (we'll still need these)
EEG_epoched1 = EEG;
EEG_epoched1.data = epoched_data1;
EEG_epoched1.pnts = samples_per_epoch;
EEG_epoched1.trials = n_epochs_segment1;

EEG_epoched2 = EEG;
EEG_epoched2.data = epoched_data2;
EEG_epoched2.pnts = samples_per_epoch;
EEG_epoched2.trials = n_epochs_segment2;

% Reference channel is the 65th (SI-ENV)
ref_channel = EEG.nbchan;

% Initialize arrays to store coherence values
itpc_vals_seg1 = zeros(n_channels, 1);
itpc_vals_seg2 = zeros(n_channels, 1);

% Calculate phase coherence for each EEG channel with respect to the reference for segment 1
for ch = 1:n_channels
    % Initialize complex phase difference array
    phase_diff = zeros(1, n_epochs_segment1);
    
    % For each epoch
    for ep = 1:n_epochs_segment1
        % Get data for current epoch
        ch_data = squeeze(epoched_data1(ch, :, ep));
        ref_data = squeeze(epoched_data1(ref_channel, :, ep));
        
        % Apply FFT 
        fft_len = 2^nextpow2(samples_per_epoch);
        fft_ch = fft(ch_data, fft_len);
        fft_ref = fft(ref_data, fft_len);
        
        % Calculate frequency resolution
        freq_res = fs / fft_len;
        
        % Find bin index for 1 Hz (the entrainment frequency)
        freq_idx = round(sine_freq / freq_res) + 1;
        
        % Get phase difference at entrainment frequency
        phase_ch = angle(fft_ch(freq_idx));
        phase_ref = angle(fft_ref(freq_idx));
        
        % Store phase difference as complex unit vector
        phase_diff(ep) = exp(1i * (phase_ch - phase_ref));
    end
    
    % ITPC = magnitude of the mean of the complex phase differences
    itpc_vals_seg1(ch) = abs(mean(phase_diff));
end

% Calculate phase coherence for segment 2
for ch = 1:n_channels
    % Initialize complex phase difference array
    phase_diff = zeros(1, n_epochs_segment2);
    
    % For each epoch
    for ep = 1:n_epochs_segment2
        % Get data for current epoch
        ch_data = squeeze(epoched_data2(ch, :, ep));
        ref_data = squeeze(epoched_data2(ref_channel, :, ep));
        
        % Apply FFT
        fft_len = 2^nextpow2(samples_per_epoch);
        fft_ch = fft(ch_data, fft_len);
        fft_ref = fft(ref_data, fft_len);
        
        % Calculate frequency resolution
        freq_res = fs / fft_len;
        
        % Find bin index for 1 Hz
        freq_idx = round(sine_freq / freq_res) + 1;
        
        % Get phase difference at entrainment frequency
        phase_ch = angle(fft_ch(freq_idx));
        phase_ref = angle(fft_ref(freq_idx));
        
        % Store phase difference as complex unit vector
        phase_diff(ep) = exp(1i * (phase_ch - phase_ref));
    end
    
    % ITPC = magnitude of the mean of the complex phase differences
    itpc_vals_seg2(ch) = abs(mean(phase_diff));
end

%% Simple Visualization of ITPC across channels
% 1. Select representative channels for visualization
% Choose a few channels across different regions (e.g., frontal, central, parietal)
% For this example, let's use channels 1, 20, 40, and 60

selected_channels = [1, 20, 40, 60]; % Four representative channels
channel_labels = {'Chan1', 'Chan20', 'Chan40', 'Chan60'};

% 2. Plot ITPC values for selected channels as a bar chart
figure('Position', [100 100 600 400]);
X = categorical(channel_labels);
Y = [itpc_vals_seg1(selected_channels), itpc_vals_seg2(selected_channels)];

b = bar(X, Y);
b(1).FaceColor = [0.8 0.8 0.8]; % Light gray for segment 1
b(2).FaceColor = [0.2 0.6 0.8]; % Blue for segment 2

% Add labels and title
xlabel('EEG Channels');
ylabel('ITPC with Sine Wave (1 Hz)');
title('Inter-Trial Phase Coherence Comparison at 1 Hz');
legend({'Segment 1 (Random)', 'Segment 2 (Entrained)'}, 'Location', 'northwest');
grid on;

% Save figure
saveas(gcf, 'itpc_selected_channels.png');

%% 3. Plot ITPC values across all channels
figure('Position', [100 100 800 400]);
subplot(1, 2, 1);
plot(1:n_channels, itpc_vals_seg1, 'bo-', 'LineWidth', 1.5);
hold on;
plot(1:n_channels, itpc_vals_seg2, 'ro-', 'LineWidth', 1.5);
xlabel('Channel Number');
ylabel('ITPC with Sine Wave (1 Hz)');
title('ITPC Comparison Across All Channels');
legend({'Segment 1 (Random)', 'Segment 2 (Entrained)'});
grid on;

% Also plot the difference to highlight the entrainment effect
subplot(1, 2, 2);
plot(1:n_channels, itpc_vals_seg2 - itpc_vals_seg1, 'ko-', 'LineWidth', 1.5);
xlabel('Channel Number');
ylabel('ITPC Difference (Seg2 - Seg1)');
title('Entrainment Effect by Channel');
grid on;

% Save figure
saveas(gcf, 'itpc_all_channels.png');

%% Analysis 3 (Improved): Topographic Visualization of Entrainment for both segments
% Define frequency bands
entrainment_band = [sine_freq-0.2 sine_freq+0.2]; % Around the entrainment frequency

% First 2 minutes
itpc_values_segment1 = zeros(1, n_channels);
% Calculate ITPC with respect to the SI-ENV channel for first segment
for ch = 1:n_channels
    % Compute cross-spectral density between channel and SI-ENV
    [Cxy, f] = mscohere(EEG.data(ch, 1:n_samples_segment), EEG.data(end, 1:n_samples_segment), hamming(fs*2), fs, [], fs);
    % Get ITPC at the entrainment frequency
    entrainment_f_idx = find(f >= entrainment_band(1) & f <= entrainment_band(2));
    itpc_values_segment1(ch) = mean(Cxy(entrainment_f_idx));
end

% Second 2 minutes
itpc_values_segment2 = zeros(1, n_channels);
% Calculate ITPC with respect to the SI-ENV channel for second segment
for ch = 1:n_channels
    % Compute cross-spectral density between channel and SI-ENV
    [Cxy, f] = mscohere(EEG.data(ch, n_samples_segment+1:end), EEG.data(end, n_samples_segment+1:end), hamming(fs*2), fs, [], fs);
    % Get ITPC at the entrainment frequency
    entrainment_f_idx = find(f >= entrainment_band(1) & f <= entrainment_band(2));
    itpc_values_segment2(ch) = mean(Cxy(entrainment_f_idx));
end

% Ensure ITPC values are within 0-1 range (coherence should already be in this range)
itpc_values_segment1 = min(max(itpc_values_segment1, 0), 1);
itpc_values_segment2 = min(max(itpc_values_segment2, 0), 1);

% Create even electrode distribution using a standard layout approach
% Number of channels
n_channels_viz = 64;  % Standard number for visualization

% Initialize channel structure 
chanlocs = [];

% Create center electrode
chanlocs(1).labels = 'Cz';
chanlocs(1).theta = 0;
chanlocs(1).radius = 0;
chanlocs(1).X = 0;
chanlocs(1).Y = 0;
chanlocs(1).Z = 1;

% Create 3 concentric rings with increasing number of electrodes
ring_radii = [0.4, 0.7, 0.9];  % Normalized radii for the rings
ring_channels = [8, 16, n_channels_viz-1-8-16];  % Number of channels in each ring

ch_idx = 2;  % Start after center
for ring = 1:3
    radius = ring_radii(ring);
    n_chans = ring_channels(ring);
    
    for i = 1:n_chans
        angle_deg = (i-1) * (360 / n_chans);
        angle_rad = angle_deg * pi/180;
        
        chanlocs(ch_idx).labels = sprintf('Ch%d', ch_idx);
        chanlocs(ch_idx).theta = angle_deg;
        chanlocs(ch_idx).radius = radius;
        
        % 3D coordinates
        chanlocs(ch_idx).X = radius * cos(angle_rad);
        chanlocs(ch_idx).Y = radius * sin(angle_rad);
        
        % Z coordinate (height) based on distance from center
        % Further from center = lower height
        chanlocs(ch_idx).Z = sqrt(1 - min(1, chanlocs(ch_idx).X^2 + chanlocs(ch_idx).Y^2));
        
        ch_idx = ch_idx + 1;
    end
end

% Trim or extend channel locations to match our data
n_viz = length(chanlocs);
n_data = length(itpc_values_segment1);

% Map our data values to visualization channels
if n_viz ~= n_data
    % Simple linear mapping 
    itpc_viz_segment1 = zeros(1, n_viz);
    itpc_viz_segment2 = zeros(1, n_viz);
    
    for i = 1:n_viz
        % Map to closest data channel
        data_idx = max(1, min(n_data, round((i-1) * n_data / n_viz) + 1));
        itpc_viz_segment1(i) = itpc_values_segment1(data_idx);
        itpc_viz_segment2(i) = itpc_values_segment2(data_idx);
    end
else
    itpc_viz_segment1 = itpc_values_segment1;
    itpc_viz_segment2 = itpc_values_segment2;
end

% Create topoplots comparing both segments
figure('Position', [100 100 800 400]);

% Common colormap and scaling
cmap = jet(64);  % Could use other colormaps like parula, viridis, etc.

% Topoplot for ITPC with respect to SI-ENV - First 2 minutes
subplot(1, 2, 1);
topoplot(itpc_viz_segment1, chanlocs, 'electrodes', 'on', 'maplimits', [0 1], ...
         'emarker', {'.', 'k', 10, 1}, 'headrad', 'rim', 'colormap', cmap, 'style', 'map');
h = colorbar;
ylabel(h, 'Phase Coherence');
title('ITPC with SI-ENV - First 2 Minutes (Random)');

% Topoplot for ITPC with respect to SI-ENV - Second 2 minutes
subplot(1, 2, 2);
topoplot(itpc_viz_segment2, chanlocs, 'electrodes', 'on', 'maplimits', [0 1], ...
         'emarker', {'.', 'k', 10, 1}, 'headrad', 'rim', 'colormap', cmap, 'style', 'map');
h = colorbar;
ylabel(h, 'Phase Coherence');
title('ITPC with SI-ENV - Second 2 Minutes (Entrained)');

% Save the figure
saveas(gcf, 'topoplots_itpc_comparison.png')


%% 1. Circular Phase Histograms (Rose Plots)
% This visualizes the distribution of phase differences on a unit circle

% We'll create rose plots for selected channels and compare segments
% First, let's collect phase differences for a few representative channels
selected_channels = [1, 20, 40, 60]; % Frontal to posterior gradient
channel_names = {'Frontal (Ch 1)', 'Central (Ch 20)', 'Parietal (Ch 40)', 'Posterior (Ch 60)'};

% Create a new figure for the rose plots
figure('Position', [100 100 900 700]);

for i = 1:length(selected_channels)
    ch = selected_channels(i);
    
    % Get phase differences for segment 1
    phase_diff_seg1 = zeros(1, n_epochs_segment1);
    for ep = 1:n_epochs_segment1
        ch_data = squeeze(epoched_data1(ch, :, ep));
        ref_data = squeeze(epoched_data1(ref_channel, :, ep));
        
        % Get FFT
        fft_len = 2^nextpow2(samples_per_epoch);
        fft_ch = fft(ch_data, fft_len);
        fft_ref = fft(ref_data, fft_len);
        
        % Find 1Hz bin
        freq_res = fs / fft_len;
        freq_idx = round(sine_freq / freq_res) + 1;
        
        % Store phase difference (not as unit vector this time)
        phase_ch = angle(fft_ch(freq_idx));
        phase_ref = angle(fft_ref(freq_idx));
        phase_diff_seg1(ep) = mod(phase_ch - phase_ref, 2*pi);
    end
    
    % Get phase differences for segment 2
    phase_diff_seg2 = zeros(1, n_epochs_segment2);
    for ep = 1:n_epochs_segment2
        ch_data = squeeze(epoched_data2(ch, :, ep));
        ref_data = squeeze(epoched_data2(ref_channel, :, ep));
        
        % Get FFT
        fft_len = 2^nextpow2(samples_per_epoch);
        fft_ch = fft(ch_data, fft_len);
        fft_ref = fft(ref_data, fft_len);
        
        % Find 1Hz bin
        freq_res = fs / fft_len;
        freq_idx = round(sine_freq / freq_res) + 1;
        
        % Store phase difference
        phase_ch = angle(fft_ch(freq_idx));
        phase_ref = angle(fft_ref(freq_idx));
        phase_diff_seg2(ep) = mod(phase_ch - phase_ref, 2*pi);
    end
    
    % Create rose plots (circular histograms)
    % Segment 1
    subplot(2, length(selected_channels), i);
    polarhistogram(phase_diff_seg1, 18, 'Normalization', 'probability', 'FaceColor', [0.8 0.8 0.8]);
    title([channel_names{i} ' - Segment 1 (Random)']);
    
    % Calculate and display circular statistics
    mean_angle1 = circ_mean(phase_diff_seg1');
    r1 = circ_r(phase_diff_seg1'); % Resultant vector length (measure of concentration)
    str1 = sprintf('Mean: %.1f°\nR: %.2f', mean_angle1 * 180/pi, r1);
    text(-1.5, 1.5, str1, 'FontSize', 9);
    
    % Segment 2
    subplot(2, length(selected_channels), i + length(selected_channels));
    polarhistogram(phase_diff_seg2, 18, 'Normalization', 'probability', 'FaceColor', [0.2 0.6 0.8]);
    title([channel_names{i} ' - Segment 2 (Entrained)']);
    
    % Calculate and display circular statistics
    mean_angle2 = circ_mean(phase_diff_seg2');
    r2 = circ_r(phase_diff_seg2');
    str2 = sprintf('Mean: %.1f°\nR: %.2f', mean_angle2 * 180/pi, r2);
    text(-1.5, 1.5, str2, 'FontSize', 9);
end

% Save the figure
saveas(gcf, 'circular_phase_histograms.png');


%% 2. Time-Resolved Phase Analysis (FIXED)
% This shows how phase relationship evolves over time in both segments

% Choose a representative frontal channel (stronger entrainment expected)
ch_to_analyze = 1;

% Create time vectors for continuous segments
t_seg1 = (0:n_samples_segment-1)/fs;
t_seg2 = (0:n_samples_segment-1)/fs;

% Filter data at 1Hz to extract phase information
% Design a narrow bandpass filter centered at 1Hz
[b, a] = butter(4, [0.9 1.1]/(fs/2), 'bandpass');

% Filter segment 1 data with verification
ch_filtered_seg1 = filtfilt(b, a, double(EEG.data(ch_to_analyze, 1:n_samples_segment)));
ref_filtered_seg1 = filtfilt(b, a, double(EEG.data(ref_channel, 1:n_samples_segment)));

% Filter segment 2 data with verification
ch_filtered_seg2 = filtfilt(b, a, double(EEG.data(ch_to_analyze, n_samples_segment+1:end)));
ref_filtered_seg2 = filtfilt(b, a, double(EEG.data(ref_channel, n_samples_segment+1:end)));

% Verify filtered data is not all zeros
disp(['Channel data check: ', num2str(sum(abs(ch_filtered_seg1)))]);
disp(['Reference data check: ', num2str(sum(abs(ref_filtered_seg1)))]);

% Extract instantaneous phase using Hilbert transform
ch_phase_seg1 = angle(hilbert(ch_filtered_seg1));
ref_phase_seg1 = angle(hilbert(ref_filtered_seg1));
phase_diff_seg1 = mod(ch_phase_seg1 - ref_phase_seg1 + 2*pi, 2*pi); % Ensure positive values

ch_phase_seg2 = angle(hilbert(ch_filtered_seg2));
ref_phase_seg2 = angle(hilbert(ref_filtered_seg2));
phase_diff_seg2 = mod(ch_phase_seg2 - ref_phase_seg2 + 2*pi, 2*pi); % Ensure positive values

% Create new figure
figure('Position', [100 100 900 600], 'Color', 'w');

% Plot a short segment (10 seconds) to see detailed behavior
subplot(2,2,1);
t_range = 1:min(10*fs, length(phase_diff_seg1)); % First 10 seconds, but check bounds
if ~isempty(t_range) && length(t_range) > 1
    % Regular phase plot
    plot(t_seg1(t_range), phase_diff_seg1(t_range), 'b-', 'LineWidth', 1.5);
    hold on;
    
    % Calculate and plot unwrapped phase
    unwrapped = unwrap(phase_diff_seg1(t_range));
    % Shift back to the range of the wrapped phase for comparison
    unwrapped = unwrapped - unwrapped(1) + phase_diff_seg1(t_range(1));
    plot(t_seg1(t_range), unwrapped, 'r--', 'LineWidth', 1);
    
    xlabel('Time (s)');
    ylabel('Phase Difference (rad)');
    title('Phase Difference in First 10s - Segment 1 (Random)');
    ylim([0 2*pi]);
    yticks([0 pi 2*pi]);
    yticklabels({'0', '\pi', '2\pi'});
    legend('Wrapped Phase', 'Unwrapped Trend');
    grid on;
else
    text(0.5, 0.5, 'Not enough data points', 'HorizontalAlignment', 'center');
end

subplot(2,2,2);
t_range = 1:min(10*fs, length(phase_diff_seg2)); % First 10 seconds, but check bounds
if ~isempty(t_range) && length(t_range) > 1
    % Regular phase plot
    plot(t_seg2(t_range), phase_diff_seg2(t_range), 'b-', 'LineWidth', 1.5);
    hold on;
    
    % Calculate and plot unwrapped phase
    unwrapped = unwrap(phase_diff_seg2(t_range));
    % Shift back to the range of the wrapped phase for comparison
    unwrapped = unwrapped - unwrapped(1) + phase_diff_seg2(t_range(1));
    plot(t_seg2(t_range), unwrapped, 'r--', 'LineWidth', 1);
    
    xlabel('Time (s)');
    ylabel('Phase Difference (rad)');
    title('Phase Difference in First 10s - Segment 2 (Entrained)');
    ylim([0 2*pi]);
    yticks([0 pi 2*pi]);
    yticklabels({'0', '\pi', '2\pi'});
    legend('Wrapped Phase', 'Unwrapped Trend');
    grid on;
else
    text(0.5, 0.5, 'Not enough data points', 'HorizontalAlignment', 'center');
end

% Calculate and plot phase stability over time (sliding window approach)
window_size = fs * 5; % 5-second sliding window
step_size = fs; % 1-second steps

% Calculate how many windows we can fit in each segment
n_windows_seg1 = max(1, floor((length(phase_diff_seg1) - window_size) / step_size) + 1);
n_windows_seg2 = max(1, floor((length(phase_diff_seg2) - window_size) / step_size) + 1);

% Preallocate arrays
plv_over_time_seg1 = zeros(1, n_windows_seg1);
plv_over_time_seg2 = zeros(1, n_windows_seg2);

% Calculate PLV for segment 1
for w = 1:n_windows_seg1
    start_idx = (w-1) * step_size + 1;
    end_idx = min(start_idx + window_size - 1, length(phase_diff_seg1));
    
    if end_idx > start_idx
        % Calculate PLV for this window
        phase_diff_window = phase_diff_seg1(start_idx:end_idx);
        complex_plv = mean(exp(1i * phase_diff_window));
        plv_over_time_seg1(w) = abs(complex_plv);
    else
        plv_over_time_seg1(w) = NaN;
    end
end

% Calculate PLV for segment 2
for w = 1:n_windows_seg2
    start_idx = (w-1) * step_size + 1;
    end_idx = min(start_idx + window_size - 1, length(phase_diff_seg2));
    
    if end_idx > start_idx
        % Calculate PLV for this window
        phase_diff_window = phase_diff_seg2(start_idx:end_idx);
        complex_plv = mean(exp(1i * phase_diff_window));
        plv_over_time_seg2(w) = abs(complex_plv);
    else
        plv_over_time_seg2(w) = NaN;
    end
end

% Plot Phase Locking Value over time for both segments
subplot(2,1,2);
% Calculate time vector for windows
t_windows1 = (0:n_windows_seg1-1) * step_size / fs;
t_windows2 = (0:n_windows_seg2-1) * step_size / fs;

% Check if we have valid data to plot
has_valid_data = ~isempty(plv_over_time_seg1) && ~all(isnan(plv_over_time_seg1)) && ...
                ~isempty(plv_over_time_seg2) && ~all(isnan(plv_over_time_seg2));

if has_valid_data
    % Clean up any NaN values
    valid_idx1 = ~isnan(plv_over_time_seg1);
    valid_idx2 = ~isnan(plv_over_time_seg2);
    
    plot(t_windows1(valid_idx1), plv_over_time_seg1(valid_idx1), 'b-', 'LineWidth', 2);
    hold on;
    plot(t_windows2(valid_idx2), plv_over_time_seg2(valid_idx2), 'r-', 'LineWidth', 2);
    
    xlabel('Time (s)');
    ylabel('Phase Locking Value');
    title('Phase Stability Over Time (5s sliding window)');
    legend('Segment 1 (Random)', 'Segment 2 (Entrained)');
    ylim([0 1]);
    grid on;
else
    text(0.5, 0.5, 'Not enough valid data for PLV calculation', 'HorizontalAlignment', 'center');
end

% Save the figure
saveas(gcf, 'time_resolved_phase_analysis.png');

% Print debug info
fprintf('Data check complete:\n');
fprintf('Phase diff seg1 range: [%f, %f]\n', min(phase_diff_seg1), max(phase_diff_seg1));
fprintf('Phase diff seg2 range: [%f, %f]\n', min(phase_diff_seg2), max(phase_diff_seg2));
fprintf('PLV seg1 range: [%f, %f]\n', min(plv_over_time_seg1), max(plv_over_time_seg1));
fprintf('PLV seg2 range: [%f, %f]\n', min(plv_over_time_seg2), max(plv_over_time_seg2));


%% Helper functions for circular statistics 
% (normally these would be from the Circular Statistics Toolbox)

function [mu] = circ_mean(alpha)
    % Calculate the mean of circular data
    n = length(alpha);
    if n < 1
        mu = NaN;
        return;
    end
    
    % Convert to complex numbers on unit circle
    z = exp(1i * alpha);
    
    % Compute mean direction
    r = sum(z) / n;
    mu = angle(r);
end

function [r] = circ_r(alpha)
    % Calculate the resultant vector length
    n = length(alpha);
    if n < 1
        r = NaN;
        return;
    end
    
    % Convert to complex numbers on unit circle
    z = exp(1i * alpha);
    
    % Compute resultant vector length
    r = abs(sum(z) / n);
end
