% ISPC Transition Region Simulation - Experimental Design
% This script simulates the effect of using a sliding window for ISPC calculation
% across stimulus stage transitions with the actual experimental design timing

%% Parameters
fs = 100;                 % Sampling rate in Hz
duration = 4*60;          % Total duration in seconds (4 minutes)
window_size = 20;         % Window size in seconds (as in the original code)
stim_freq = 1;            % Stimulus frequency in Hz

% Define experimental events (in seconds)
% 45s pre-stim + 2min stimulus + 45s post-stim = 4min total
stim_start = 45;          % Stimulus starts after 45 seconds
stim_end = stim_start + 2*60; % Stimulus lasts for 2 minutes

% Define segment boundaries (in samples)
% Pre-stim: 45s before stim_start
pre_stim_start = 0;
pre_stim_end = stim_start - 1;

% Early-stim: first 45s of stimulation
early_stim_start = stim_start;
early_stim_end = stim_start + 45;

% Late-stim: last 45s of stimulation
late_stim_start = stim_end - 45;
late_stim_end = stim_end;

% Post-stim: 45s after stimulation
post_stim_start = stim_end + 1;
post_stim_end = stim_end + 45;

% Generate time vector
t = (0:fs*duration-1)/fs;

% Create stimulus signal (constant 1 Hz sine wave)
stim_signal = sin(2*pi*stim_freq*t);

% Create EEG signal with different phase synchronization in different segments
eeg_signal = zeros(size(t));
phase_noise = 0.3;    % Phase noise level

% Helper function for random normal distribution
function r = randn_simple()
    % Box-Muller transform
    u = rand();
    v = rand();
    r = sqrt(-2*log(u)) * cos(2*pi*v);
end

% Pre-stim: low synchronization
pre_stim_phase = pi/3;    % Phase offset in pre-stim period
for i = round(pre_stim_start*fs)+1:round(pre_stim_end*fs)
    eeg_signal(i) = sin(2*pi*stim_freq*t(i) + pre_stim_phase + randn_simple()*phase_noise);
end

% Early-stim: medium synchronization
early_stim_phase = pi/6;  % Phase offset in early-stim period
for i = round(early_stim_start*fs)+1:round(early_stim_end*fs)
    eeg_signal(i) = sin(2*pi*stim_freq*t(i) + early_stim_phase + randn_simple()*phase_noise*0.7);
end

% Middle segment (not analyzed): transitional synchronization
middle_start = early_stim_end + 1;
middle_end = late_stim_start - 1;
middle_phase = pi/8;      % Transitional phase
for i = round(middle_start*fs)+1:round(middle_end*fs)
    % Linear transition from early to late stim phase
    progress = (i - round(middle_start*fs)) / (round(middle_end*fs) - round(middle_start*fs));
    current_phase = early_stim_phase * (1-progress) + late_stim_phase * progress;
    current_noise = phase_noise * (0.7 * (1-progress) + 0.3 * progress);
    eeg_signal(i) = sin(2*pi*stim_freq*t(i) + current_phase + randn_simple()*current_noise);
end

% Late-stim: high synchronization
late_stim_phase = pi/12;  % Phase offset in late-stim period
for i = round(late_stim_start*fs)+1:round(late_stim_end*fs)
    eeg_signal(i) = sin(2*pi*stim_freq*t(i) + late_stim_phase + randn_simple()*phase_noise*0.3);
end

% Post-stim: medium to low synchronization
post_stim_phase = pi/4;   % Phase offset in post-stim period
for i = round(post_stim_start*fs)+1:round(post_stim_end*fs)
    eeg_signal(i) = sin(2*pi*stim_freq*t(i) + post_stim_phase + randn_simple()*phase_noise*0.8);
end

%% Calculate phase angles (using Hilbert transform)
% Extract phase angles using Hilbert transform
stim_analytic = hilbert(stim_signal);
eeg_analytic = hilbert(eeg_signal);

stim_phase = angle(stim_analytic);
eeg_phase = angle(eeg_analytic);

% Phase difference
phase_diff = eeg_phase - stim_phase;

%% Calculate ISPC using sliding window (standard approach)
% Window size in samples
half_window = round(window_size * fs / 2);
ispc = zeros(size(t));

% Standard sliding window approach (as in original code)
for i = 1:length(t)
    % Window boundaries
    win_start = max(1, i - half_window);
    win_end = min(length(t), i + half_window);
    
    % Calculate ISPC for this window
    ispc(i) = abs(mean(exp(1i * phase_diff(win_start:win_end))));
end

%% Calculate ISPC with segment-aware windows (proposed solution)
ispc_segment_aware = zeros(size(t));

% Define segment boundaries as a structure (similar to the actual code)
segments = struct();
segments.pre_stim = [round(pre_stim_start*fs)+1, round(pre_stim_end*fs)];
segments.early_stim = [round(early_stim_start*fs)+1, round(early_stim_end*fs)];
segments.late_stim = [round(late_stim_start*fs)+1, round(late_stim_end*fs)];
segments.post_stim = [round(post_stim_start*fs)+1, round(post_stim_end*fs)];

% Calculate ISPC for each segment separately
segment_names = fieldnames(segments);
for s = 1:length(segment_names)
    seg_name = segment_names{s};
    seg_start = segments.(seg_name)(1);
    seg_end = segments.(seg_name)(2);
    
    % For each time point in this segment
    for i = seg_start:seg_end
        % Window boundaries constrained to this segment
        win_start = max(seg_start, i - half_window);
        win_end = min(seg_end, i + half_window);
        
        % Calculate ISPC for this constrained window
        ispc_segment_aware(i) = abs(mean(exp(1i * phase_diff(win_start:win_end))));
    end
end

%% Calculate segment means (center-only approach)
% Define the center portion size (e.g., middle 50% of each segment)
center_portion = 0.5;
segment_means_center = struct();

% Process each segment
for s = 1:length(segment_names)
    seg_name = segment_names{s};
    seg_start = segments.(seg_name)(1);
    seg_end = segments.(seg_name)(2);
    
    % Calculate segment size and center region
    segment_size = seg_end - seg_start + 1;
    margin = round((1 - center_portion) * segment_size / 2);
    
    center_start = seg_start + margin;
    center_end = seg_end - margin;
    
    % Calculate mean ISPC for center portion
    segment_means_center.(seg_name) = mean(ispc(center_start:center_end));
end

%% Create visualization
figure('Position', [100, 100, 1200, 800]);

% Plot signals
subplot(4,1,1);
plot(t, stim_signal, 'b', t, eeg_signal, 'r');
xlim([0, duration]);
title('Simulated Signals');
legend('Stimulus', 'EEG');
ylabel('Amplitude');
hold on;

% Add vertical lines for segment boundaries and event markers
xline(stim_start, '-g', 'Stim Start', 'LineWidth', 2);
xline(stim_end, '-r', 'Stim End', 'LineWidth', 2);
xline(early_stim_end, '--m');
xline(late_stim_start, '--c');
text(stim_start + 22.5, 1.5, 'Early-Stim (45s)', 'HorizontalAlignment', 'center');
text(stim_start + (stim_end-stim_start)/2, 1.5, 'Middle (30s)', 'HorizontalAlignment', 'center');
text(stim_end - 22.5, 1.5, 'Late-Stim (45s)', 'HorizontalAlignment', 'center');
text(pre_stim_start + 22.5, 1.5, 'Pre-Stim (45s)', 'HorizontalAlignment', 'center');
text(post_stim_start + 22.5, 1.5, 'Post-Stim (45s)', 'HorizontalAlignment', 'center');

% Plot phase differences
subplot(4,1,2);
plot(t, phase_diff);
xlim([0, duration]);
title('Phase Difference');
ylabel('Phase (rad)');
hold on;
xline(stim_start, '-g', 'LineWidth', 2);
xline(stim_end, '-r', 'LineWidth', 2);
xline(early_stim_end, '--m');
xline(late_stim_start, '--c');

% Plot original ISPC
subplot(4,1,3);
plot(t, ispc);
xlim([0, duration]);
ylim([0, 1]);
title('ISPC with Standard Sliding Window');
ylabel('ISPC');
hold on;
xline(stim_start, '-g', 'LineWidth', 2);
xline(stim_end, '-r', 'LineWidth', 2);
xline(early_stim_end, '--m');
xline(late_stim_start, '--c');

% Draw window size at a reference point to illustrate the issue
ref_point = round(stim_start*fs) + 1;
% Make sure window bounds don't go out of array bounds
win_start_idx = max(1, ref_point-half_window);
win_end_idx = min(length(t), ref_point+half_window);
plot([t(win_start_idx), t(win_end_idx)], [0.9, 0.9], 'r-', 'LineWidth', 2);
text(t(ref_point), 0.95, 'Window Size', 'HorizontalAlignment', 'center');

% Plot segment-aware ISPC
subplot(4,1,4);
plot(t, ispc_segment_aware);
xlim([0, duration]);
ylim([0, 1]);
title('ISPC with Segment-Aware Windows');
xlabel('Time (s)');
ylabel('ISPC');
hold on;
xline(stim_start, '-g', 'LineWidth', 2);
xline(stim_end, '-r', 'LineWidth', 2);
xline(early_stim_end, '--m');
xline(late_stim_start, '--c');

% Add segment labels
for s = 1:length(segment_names)
    seg_name = segment_names{s};
    seg_center = mean(segments.(seg_name));
    text(t(round(seg_center)), 0.9, seg_name, 'HorizontalAlignment', 'center');
end

%% Display quantitative comparison
% Calculate average ISPC for each segment and method
fprintf('=== ISPC Averages By Segment ===\n');
for s = 1:length(segment_names)
    seg_name = segment_names{s};
    seg_start = segments.(seg_name)(1);
    seg_end = segments.(seg_name)(2);
    
    % Standard ISPC
    std_mean = mean(ispc(seg_start:seg_end));
    % Segment-aware ISPC
    sa_mean = mean(ispc_segment_aware(seg_start:seg_end));
    % Center-only
    center_mean = segment_means_center.(seg_name);
    
    fprintf('Segment: %s\n', seg_name);
    fprintf('  Standard ISPC mean: %.4f\n', std_mean);
    fprintf('  Segment-aware ISPC mean: %.4f\n', sa_mean);
    fprintf('  Center-only ISPC mean: %.4f\n', center_mean);
    fprintf('  Diff (Segment-aware vs Standard): %.4f\n', sa_mean - std_mean);
    fprintf('  Diff (Center-only vs Standard): %.4f\n', center_mean - std_mean);
    fprintf('\n');
end

%% Create comparison figure to highlight transition regions
figure('Position', [100, 500, 1200, 400]);
plot(t, ispc, 'b-', t, ispc_segment_aware, 'r--', 'LineWidth', 1.5);
xlim([0, duration]);
ylim([0, 1]);
title('Comparison of Standard ISPC vs. Segment-Aware ISPC');
xlabel('Time (s)');
ylabel('ISPC');
legend('Standard Window', 'Segment-Aware Window');
hold on;

% Add segment boundaries
xline(stim_start, '-g', 'Stim Start', 'LineWidth', 2);
xline(stim_end, '-r', 'Stim End', 'LineWidth', 2);
xline(early_stim_end, '--m', 'Early-Stim End');
xline(late_stim_start, '--c', 'Late-Stim Start');

% Highlight transition regions - fixing array bound issues
transition_width = half_window;

% Make sure indices are within bounds for pre-stim to early-stim transition
pre_early_center = round(stim_start*fs);
pre_early_start = max(1, pre_early_center - half_window);
pre_early_end = min(length(t), pre_early_center + half_window);
rectangle('Position', [t(pre_early_start), 0, t(pre_early_end)-t(pre_early_start), 1], ...
    'FaceColor', [1 0 0 0.1], 'EdgeColor', 'none');

% Make sure indices are within bounds for early-stim to middle transition
early_mid_center = round(early_stim_end*fs);
early_mid_start = max(1, early_mid_center - half_window);
early_mid_end = min(length(t), early_mid_center + half_window);
rectangle('Position', [t(early_mid_start), 0, t(early_mid_end)-t(early_mid_start), 1], ...
    'FaceColor', [0 1 0 0.1], 'EdgeColor', 'none');

% Make sure indices are within bounds for middle to late-stim transition
mid_late_center = round(late_stim_start*fs);
mid_late_start = max(1, mid_late_center - half_window);
mid_late_end = min(length(t), mid_late_center + half_window);
rectangle('Position', [t(mid_late_start), 0, t(mid_late_end)-t(mid_late_start), 1], ...
    'FaceColor', [0 0 1 0.1], 'EdgeColor', 'none');

% Make sure indices are within bounds for late-stim to post-stim transition
late_post_center = round(stim_end*fs);
late_post_start = max(1, late_post_center - half_window);
late_post_end = min(length(t), late_post_center + half_window);
rectangle('Position', [t(late_post_start), 0, t(late_post_end)-t(late_post_start), 1], ...
    'FaceColor', [1 0 1 0.1], 'EdgeColor', 'none');

% Add labels for transition regions
text(stim_start, 0.2, 'Transition Region', 'HorizontalAlignment', 'center', 'Color', 'r');
text(early_stim_end, 0.2, 'Transition Region', 'HorizontalAlignment', 'center', 'Color', 'g');
text(late_stim_start, 0.2, 'Transition Region', 'HorizontalAlignment', 'center', 'Color', 'b');
text(stim_end, 0.2, 'Transition Region', 'HorizontalAlignment', 'center', 'Color', 'm');
