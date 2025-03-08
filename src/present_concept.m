%% EEG Phase Coherence Visualization for Teaching Purposes
% This script creates conceptual visualizations to illustrate phase coherence
% analysis between EEG signals and a 1Hz sine wave stimulus
% It visualizes:
% 1. Continuous EEG with event markers and epochs
% 2. Superimposed epochs aligned to the 1Hz stimulus
% 3. Vector representation of phase relationships and their average


%% Initialize EEGLAB
clear; close all; clc;
% Make sure EEGLAB is in your MATLAB path, otherwise you'll need to add it:
addpath('/Users/idohaber/Documents/MATLAB/eeglab2024.0/');
eeglab;
% load data
[EEG] = pop_loadset('../data/nrem_sleep_with_entrainment_contrast.set');

% Extract basic parameters
fs = EEG.srate;                    % Sampling frequency
total_duration = EEG.pnts/fs;      % Total duration in seconds
segment_duration = total_duration/2; % Each segment duration
n_channels = EEG.nbchan - 1;       % Number of EEG channels (excluding sine wave channel)
sine_channel = EEG.nbchan;         % Index of the sine wave channel
sine_freq = 1;                     % Frequency of entrainment (Hz)

% Time vector for each segment
segment_samples = round(segment_duration * fs);
t_full = (0:EEG.pnts-1)/fs;        % Time vector for full recording
t_seg1 = t_full(1:segment_samples);                % Time for segment 1
t_seg2 = t_full(segment_samples+1:end);    % Time for segment 2

%% PART 1: Visualize Continuous EEG with Event Markers and Epochs
% We'll create simulated events at the peaks of the sine wave for epoching

% 1.1 Find peaks in the sine wave to create event markers
sine_wave = EEG.data(sine_channel, :);
% Use findpeaks to detect positive peaks in the sine wave
[~, peak_locs] = findpeaks(sine_wave, 'MinPeakDistance', fs*0.9); % Minimum distance ~0.9s (for 1Hz)

% 1.2 Create a figure showing continuous EEG with event markers
figure('Position', [100 100 1200 800], 'Color', 'w');

% Choose channels to display (frontal, central and posterior)
display_channels = [1, round(n_channels/2), n_channels];
channel_names = {'Frontal', 'Central', 'Posterior'};

% Define colors for visualization
segment1_color = [0.8 0.8 0.8]; % Light gray for segment 1
segment2_color = [0.2 0.6 0.8]; % Blue for segment 2

% Plot section: Continuous EEG with event markers
subplot(3, 1, 1);
hold on;

% Plot vertical lines separating the two segments
plot([segment_duration segment_duration], [-200 200], 'r--', 'LineWidth', 1.5);
text(segment_duration/2, 180, 'Random EEG', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(segment_duration + segment_duration/2, 180, 'Entrained EEG', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

% Plot selected EEG channels (with offsets for visibility)
offsets = [100, 0, -100]; % Vertical offsets for each channel
for i = 1:length(display_channels)
    ch = display_channels(i);
    % Get channel data with offset
    ch_data = EEG.data(ch, :) + offsets(i);
    
    % Plot first segment
    plot(t_seg1, ch_data(1:segment_samples), 'Color', segment1_color, 'LineWidth', 1);
    
    % Plot second segment
    plot(t_seg2, ch_data(segment_samples+1:end), 'Color', segment2_color, 'LineWidth', 1);
    
    % Add channel label
    text(0, offsets(i), channel_names{i}, 'FontWeight', 'bold');
end

% Plot the sine wave at the bottom
sine_wave_scaled = sine_wave/5 - 175; % Scale and offset for visualization
plot(t_full, sine_wave_scaled, 'k', 'LineWidth', 1.5);
text(0, -175, 'Sine Wave', 'FontWeight', 'bold');

% Plot event markers as vertical lines
for i = 1:length(peak_locs)
    event_time = peak_locs(i)/fs;
    % Different color for events in different segments
    if event_time <= segment_duration
        line_color = 'k';
    else
        line_color = 'r';
    end
    
    % Only plot a selection of events to avoid cluttering
    if mod(i, 5) == 0
        plot([event_time event_time], [-200 150], '--', 'Color', line_color, 'LineWidth', 0.5);
    end
end

% Add some example epoch boxes (green rectangles)
epoch_duration = 2; % 2 seconds per epoch
example_epochs = [30, 100, 170, 220]; % Times in seconds for examples
for ep_time = example_epochs
    rectangle('Position', [ep_time, -200, epoch_duration, 350], ...
              'EdgeColor', 'g', 'LineWidth', 2, 'LineStyle', '-');
end

% Annotate one of the epochs
annotation('textarrow', [0.42 0.38], [0.68 0.65], 'String', 'Epoch', ...
           'FontWeight', 'bold', 'FontSize', 12);

% Formatting
xlim([0 total_duration]);
ylim([-200 200]);
xlabel('Time (s)');
ylabel('Amplitude (μV)');
title('Continuous EEG with Event Markers and Example Epochs');
grid on;

%% PART 2: Visualize Superimposed Epochs Aligned to Sine Wave
% For this visualization, we'll extract epochs around each sine wave peak
% and superimpose them

% 2.1 Extract epochs from both segments
epoch_duration = 2; % seconds
epoch_samples = round(epoch_duration * fs);
half_epoch = round(epoch_samples/2);

% Get peaks from first segment (random EEG)
segment1_peaks = peak_locs(peak_locs <= segment_samples);
% Select a subset of peaks to reduce computation
segment1_peaks = segment1_peaks(1:5:end);

% Get peaks from second segment (entrained EEG)
segment2_peaks = peak_locs(peak_locs > segment_samples);
% Select a subset of peaks to reduce computation
segment2_peaks = segment2_peaks(1:5:end);

% Extract epochs for a chosen EEG channel (e.g., frontal channel 1)
chosen_channel = 1;

% Function to extract epochs around peak locations
function epochs = extract_epochs(data, peak_locs, half_epoch)
    epochs = zeros(length(peak_locs), half_epoch*2+1);
    for i = 1:length(peak_locs)
        peak = peak_locs(i);
        if peak > half_epoch && peak <= length(data) - half_epoch
            epochs(i, :) = data(peak-half_epoch:peak+half_epoch);
        end
    end
end

% Extract epochs around peaks for both segments
epochs_seg1 = extract_epochs(EEG.data(chosen_channel, :), segment1_peaks, half_epoch);
epochs_seg2 = extract_epochs(EEG.data(chosen_channel, :), segment2_peaks, half_epoch);

% Extract sine wave epochs for reference
sine_epochs_seg1 = extract_epochs(EEG.data(sine_channel, :), segment1_peaks, half_epoch);
sine_epochs_seg2 = extract_epochs(EEG.data(sine_channel, :), segment2_peaks, half_epoch);

% Time vector for epoch plotting
epoch_time = (-half_epoch:half_epoch)/fs;

% 2.2 Plot superimposed epochs
subplot(3, 2, 3);
% First segment epochs (random)
hold on;
% Plot individual epochs in light color
for i = 1:size(epochs_seg1, 1)
    plot(epoch_time, epochs_seg1(i, :), 'Color', [0.8 0.8 0.8, 0.3], 'LineWidth', 0.5);
end
% Plot average epoch in bold
plot(epoch_time, mean(epochs_seg1, 1), 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
% Plot average sine wave
plot(epoch_time, mean(sine_epochs_seg1, 1), 'r', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Amplitude (μV)');
title('Superimposed Epochs - Random EEG (Segment 1)');
grid on;

subplot(3, 2, 4);
% Second segment epochs (entrained)
hold on;
% Plot individual epochs in light color
for i = 1:size(epochs_seg2, 1)
    plot(epoch_time, epochs_seg2(i, :), 'Color', [0.2 0.6 0.8, 0.3], 'LineWidth', 0.5);
end
% Plot average epoch in bold
plot(epoch_time, mean(epochs_seg2, 1), 'Color', [0.1 0.4 0.6], 'LineWidth', 2);
% Plot average sine wave
plot(epoch_time, mean(sine_epochs_seg2, 1), 'r', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Amplitude (μV)');
title('Superimposed Epochs - Entrained EEG (Segment 2)');
grid on;

%% PART 3: Vector Representation of Phase Relationships
% This visualization shows phase as vectors on a unit circle

% 3.1 Calculate phase at peak locations (when sine wave = max amplitude)
% For simplicity, we'll use a Hilbert transform to extract the phase

% Filter data at 1Hz (to extract matching frequency components)
[b, a] = butter(4, [0.8 1.2]/(fs/2), 'bandpass');

% Filter the data
eeg_filtered_seg1 = filtfilt(b, a, double(EEG.data(chosen_channel, 1:segment_samples)));
eeg_filtered_seg2 = filtfilt(b, a, double(EEG.data(chosen_channel, segment_samples+1:end)));
sine_filtered_seg1 = filtfilt(b, a, double(EEG.data(sine_channel, 1:segment_samples)));
sine_filtered_seg2 = filtfilt(b, a, double(EEG.data(sine_channel, segment_samples+1:end)));

% Get phase using Hilbert transform
eeg_phase_seg1 = angle(hilbert(eeg_filtered_seg1));
eeg_phase_seg2 = angle(hilbert(eeg_filtered_seg2));
sine_phase_seg1 = angle(hilbert(sine_filtered_seg1));
sine_phase_seg2 = angle(hilbert(sine_filtered_seg2));

% Get phase differences at peak locations
segment1_peak_indices = segment1_peaks(segment1_peaks <= length(eeg_phase_seg1));
phase_diff_seg1 = zeros(1, length(segment1_peak_indices));
for i = 1:length(segment1_peak_indices)
    idx = segment1_peak_indices(i);
    if idx <= length(eeg_phase_seg1)
        phase_diff_seg1(i) = mod(eeg_phase_seg1(idx) - sine_phase_seg1(idx) + 2*pi, 2*pi);
    end
end

segment2_peak_indices = segment2_peaks - segment_samples;
segment2_peak_indices = segment2_peak_indices(segment2_peak_indices > 0 & segment2_peak_indices <= length(eeg_phase_seg2));
phase_diff_seg2 = zeros(1, length(segment2_peak_indices));
for i = 1:length(segment2_peak_indices)
    idx = segment2_peak_indices(i);
    if idx <= length(eeg_phase_seg2)
        phase_diff_seg2(i) = mod(eeg_phase_seg2(idx) - sine_phase_seg2(idx) + 2*pi, 2*pi);
    end
end

% 3.2 Plot phase vectors on a unit circle for both segments
% Segment 1 (Random EEG)
subplot(3, 2, 5);
polarplot([0 0], [0 1.1], 'k-', 'LineWidth', 0.5); % Draw an axis line
hold on;

% Plot individual phase vectors
for i = 1:length(phase_diff_seg1)
    phase = phase_diff_seg1(i);
    % Plot a vector from center to edge of circle
    polarplot([0 phase], [0 1], 'Color', [0.8 0.8 0.8, 0.5], 'LineWidth', 0.5);
end

% Plot mean phase vector
mean_phase_seg1 = circ_mean(phase_diff_seg1');
mean_length_seg1 = circ_r(phase_diff_seg1');

polarplot([0 mean_phase_seg1], [0 mean_length_seg1], 'k-', 'LineWidth', 3);
polarplot(mean_phase_seg1, mean_length_seg1, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);

title(sprintf('Phase Vectors - Random EEG (R = %.2f)', mean_length_seg1));

% Segment 2 (Entrained EEG)
subplot(3, 2, 6);
polarplot([0 0], [0 1.1], 'k-', 'LineWidth', 0.5); % Draw an axis line
hold on;

% Plot individual phase vectors
for i = 1:length(phase_diff_seg2)
    phase = phase_diff_seg2(i);
    % Plot a vector from center to edge of circle
    polarplot([0 phase], [0 1], 'Color', [0.2 0.6 0.8, 0.5], 'LineWidth', 0.5);
end

% Plot mean phase vector
mean_phase_seg2 = circ_mean(phase_diff_seg2');
mean_length_seg2 = circ_r(phase_diff_seg2');

polarplot([0 mean_phase_seg2], [0 mean_length_seg2], 'k-', 'LineWidth', 3);
polarplot(mean_phase_seg2, mean_length_seg2, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);

title(sprintf('Phase Vectors - Entrained EEG (R = %.2f)', mean_length_seg2));

% Adjust the overall figure
sgtitle('Visualization of Phase Coherence Analysis in EEG', 'FontSize', 16, 'FontWeight', 'bold');

% Save the figure for presentations
saveas(gcf, 'eeg_phase_coherence_visualization.png');
saveas(gcf, 'eeg_phase_coherence_visualization.fig');

%% Helper function for circular statistics
function [mu] = circ_mean(alpha)
    % Calculate the mean direction of circular data
    z = exp(1i * alpha);
    mu = angle(mean(z));
end

function [r] = circ_r(alpha)
    % Calculate the resultant vector length (phase consistency)
    z = exp(1i * alpha);
    r = abs(mean(z));
end