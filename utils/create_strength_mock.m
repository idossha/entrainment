%% Generate NREM EEG Mock Data with Event Markers and Entrainment
% This script creates 40 minutes of simulated NREM EEG data with
% 'stim start' and 'stim end' event markers and entrainment during stimulation periods

%% Initialize EEGLAB
clear; close all; clc;
% Make sure EEGLAB is in your MATLAB path, otherwise you'll need to add it:
addpath('/Users/idohaber/Documents/MATLAB/eeglab2024.0/');
addpath('/Users/idohaber/Git-Projects/entrainment/src/');
addpath('/Users/idohaber/Git-Projects/entrainment/data/');
eeglab nogui;

%% Parameters for simulated EEG data
fs = 500;                    % Sampling frequency (Hz)
total_duration = 40*60;      % 40 minutes in seconds
n_channels = 64;             % Number of main EEG channels
sine_freq = 1;               % Frequency of the sine wave (Hz)

% Protocol timing parameters
protocol_duration = 200;     % 200 seconds between stim start and stim end
min_between_protocols = 400; % Minimum 400 seconds between protocols
n_protocols = 3;             % Number of protocols
min_start_time = 200;        % Start first protocol at least 3 minutes into the recording
min_end_time = 3*60;         % End last protocol at least 3 minutes before end of recording

% Entrainment parameters
max_entrainment = 0.8;       % Maximum entrainment strength (0-1)
transition_time = 5;         % Seconds to transition entrainment at edges

%% Create time vector and determine sample counts
t = 0:1/fs:total_duration-1/fs;   % Full time vector
n_samples = length(t);            % Total number of samples

%% Generate event timing
% Ensure we have enough time for all protocols with proper spacing
min_total_duration = (n_protocols * protocol_duration) + ...
                    ((n_protocols-1) * min_between_protocols) + ...
                    min_start_time + min_end_time;
if total_duration < min_total_duration
    error('Total duration too short for specified protocols and spacing');
end

% Calculate available space for distributing protocols
available_space = total_duration - (n_protocols * protocol_duration) - min_start_time - min_end_time;

% Generate protocol start times
protocol_starts = zeros(1, n_protocols);
protocol_ends = zeros(1, n_protocols);

% First protocol starts after initial delay
protocol_starts(1) = min_start_time;
protocol_ends(1) = protocol_starts(1) + protocol_duration;

% Distribute remaining protocols with at least min_between_protocols spacing
for i = 2:n_protocols
    min_start = protocol_ends(i-1) + min_between_protocols;
    
    % For the last protocol, ensure it ends at least min_end_time before the end
    if i == n_protocols
        max_start = total_duration - protocol_duration - min_end_time;
    else
        max_start = total_duration - protocol_duration - min_end_time - ...
                   (n_protocols-i)*(min_between_protocols + protocol_duration);
    end
    
    % Random start time within allowed range
    protocol_starts(i) = min_start + rand * (max_start - min_start);
    protocol_ends(i) = protocol_starts(i) + protocol_duration;
end

% Convert to sample indices
protocol_start_samples = round(protocol_starts * fs);
protocol_end_samples = round(protocol_ends * fs);

%% Generate simulated NREM sleep EEG data
% Create the 64 channels of simulated NREM sleep data
eeg_data = zeros(n_channels, n_samples);

% Create 65th channel with 1Hz sine wave for entrainment (present throughout)
sine_wave = 10 * sin(2*pi*sine_freq*t);
sine_wave = sine_wave + 0.5*randn(size(sine_wave)); % Add some noise

for ch = 1:n_channels
    % Generate dominant delta waves (0.5-4 Hz) - NREM characteristic
    delta = 15 * filter_bandpass(randn(1, n_samples), fs, 0.5, 4);
    
    % Add some theta waves (4-8 Hz) - less power than delta
    theta = 5 * filter_bandpass(randn(1, n_samples), fs, 4, 8);
    
    % Add minimal alpha waves (8-12 Hz)
    alpha = 2 * filter_bandpass(randn(1, n_samples), fs, 8, 12);
    
    % Add minimal beta waves (13-30 Hz)
    beta = 1 * filter_bandpass(randn(1, n_samples), fs, 13, 30);
    
    % Add minimal gamma waves (>30 Hz)
    gamma = 0.5 * filter_bandpass(randn(1, n_samples), fs, 30, 45);
    
    % Combine all frequency bands to create realistic NREM EEG
    eeg_data(ch, :) = delta + theta + alpha + beta + gamma;
    
    % Add some spatial correlation between channels (nearby channels are more correlated)
    if ch > 1
        % Add correlation with neighboring channels
        neighbor_contrib = 0.3 * eeg_data(max(1, ch-1), :) + 0.2 * eeg_data(max(1, ch-2), :);
        eeg_data(ch, :) = 0.5 * eeg_data(ch, :) + neighbor_contrib;
    end
    
    % Add some random channel-specific characteristics
    eeg_data(ch, :) = eeg_data(ch, :) * (0.9 + 0.2*rand());
end

% Add occasional K-complexes (large amplitude waveforms) typical in NREM
n_k_complexes = 300; % Number of K-complexes in the whole recording
k_complex_length = round(1.0 * fs); % Approximately 1 second long
k_complex_times = randi([1, n_samples - k_complex_length], 1, n_k_complexes);

for k = 1:n_k_complexes
    % K-complex waveform (sharp negative followed by positive)
    k_wave = 30 * exp(-((1:k_complex_length) - k_complex_length/4).^2 / (2*(k_complex_length/8)^2));
    k_wave(1:round(k_complex_length/3)) = -k_wave(1:round(k_complex_length/3));
    
    % Add to random subset of channels (frontal channels more likely)
    affected_channels = find(rand(1, n_channels) < (0.5 + 0.5*(1:n_channels)/n_channels));
    for ch = affected_channels
        eeg_data(ch, k_complex_times(k):(k_complex_times(k)+k_complex_length-1)) = ...
            eeg_data(ch, k_complex_times(k):(k_complex_times(k)+k_complex_length-1)) + ...
            k_wave * (0.7 + 0.6*rand()); % Vary amplitude slightly
    end
end

% Add occasional sleep spindles (burst of 12-14 Hz activity) typical in NREM
n_spindles = 500; % Number of spindles in the whole recording
spindle_length = round(0.8 * fs); % Approximately 0.8 seconds long
spindle_times = randi([1, n_samples - spindle_length], 1, n_spindles);

for s = 1:n_spindles
    % Create spindle (amplitude-modulated 13Hz oscillation)
    spindle_freq = 13; % Hz
    spindle_time = (0:spindle_length-1)/fs;
    spindle_wave = sin(2*pi*spindle_freq*spindle_time);
    spindle_envelope = sin(pi*spindle_time/max(spindle_time)).^2;
    spindle_wave = 5 * spindle_wave .* spindle_envelope;
    
    % Add to random subset of channels (central and parietal channels more likely)
    channel_weights = sin(pi*(1:n_channels)/n_channels).^2;
    affected_channels = find(rand(1, n_channels) < channel_weights);
    for ch = affected_channels
        eeg_data(ch, spindle_times(s):(spindle_times(s)+spindle_length-1)) = ...
            eeg_data(ch, spindle_times(s):(spindle_times(s)+spindle_length-1)) + ...
            spindle_wave * (0.8 + 0.4*rand()); % Vary amplitude slightly
    end
end

%% Add entrainment to the 1Hz signal during stimulation periods
% Create an entrainment mask (0 = no entrainment, 1 = full entrainment)
entrainment_mask = zeros(1, n_samples);

% Transition sample count (for smooth ramp up/down)
transition_samples = round(transition_time * fs);

for p = 1:n_protocols
    % Stimulation period
    start_idx = protocol_start_samples(p);
    end_idx = protocol_end_samples(p);
    
    % Create linear ramp up/down to avoid sharp transitions
    ramp_up = linspace(0, 1, transition_samples);
    ramp_down = linspace(1, 0, transition_samples);
    
    % Apply ramp up at start of stimulation (ensure we don't go out of bounds)
    ramp_up_end = min(start_idx + transition_samples - 1, n_samples);
    entrainment_mask(start_idx:ramp_up_end) = ramp_up(1:(ramp_up_end-start_idx+1));
    
    % Apply full entrainment for middle of stimulation
    middle_start = start_idx + transition_samples;
    middle_end = end_idx - transition_samples;
    if middle_start <= middle_end
        entrainment_mask(middle_start:middle_end) = 1;
    end
    
    % Apply ramp down at end of stimulation (ensure we don't go out of bounds)
    ramp_down_start = max(1, end_idx - transition_samples + 1);
    ramp_down_end = min(end_idx, n_samples);
    down_range = (ramp_down_start:ramp_down_end) - (end_idx - transition_samples);
    entrainment_mask(ramp_down_start:ramp_down_end) = ramp_down(down_range);
end

% Apply entrainment to each channel with frontal channels (lower numbers) 
% having stronger entrainment
for ch = 1:n_channels
    % Calculate channel-specific entrainment strength
    % Frontal channels (lower numbers) have stronger entrainment
    channel_entrainment = max_entrainment * (1 - min(1, ch/n_channels)*0.7);
    
    % Apply entrainment - add a portion of the sine wave to the EEG channel
    % Scale by the entrainment mask to only apply during stimulation
    eeg_data(ch, :) = eeg_data(ch, :) + sine_wave .* entrainment_mask * channel_entrainment;
end

% Combine all channels including the sine wave channel
all_data = [eeg_data; sine_wave];

%% Create an EEGLAB data structure
EEG = eeg_emptyset();
EEG.data = all_data;
EEG.srate = fs;
EEG.pnts = n_samples;
EEG.trials = 1;
EEG.xmin = 0;
EEG.xmax = total_duration;
EEG.times = t * 1000; % Convert to milliseconds for EEGLAB
EEG.nbchan = size(all_data, 1);

% Create channel labels
for ch = 1:n_channels
    EEG.chanlocs(ch).labels = sprintf('Chan%d', ch);
end
EEG.chanlocs(n_channels+1).labels = 'SI-ENV';

% Create events
event_counter = 1;
for p = 1:n_protocols
    % Start event
    EEG.event(event_counter).type = 'stim start';
    EEG.event(event_counter).latency = protocol_start_samples(p);
    EEG.event(event_counter).duration = 1;
    event_counter = event_counter + 1;
    
    % End event
    EEG.event(event_counter).type = 'stim end';
    EEG.event(event_counter).latency = protocol_end_samples(p);
    EEG.event(event_counter).duration = 1;
    event_counter = event_counter + 1;
end

% Update the EEGLAB data structure
EEG = eeg_checkset(EEG);

% Print protocol timing information
fprintf('Protocol timings:\n');
for p = 1:n_protocols
    fprintf('Protocol %d: Start = %.1f min, End = %.1f min (Duration = %.1f s)\n', ...
        p, protocol_starts(p)/60, protocol_ends(p)/60, protocol_duration);
    
    if p < n_protocols
        gap = protocol_starts(p+1) - protocol_ends(p);
        fprintf('  Gap to next protocol: %.1f s (%.1f min)\n', gap, gap/60);
    end
end

% Print time from last protocol end to end of recording
last_protocol_end_time = protocol_ends(end);
time_after_last_protocol = total_duration - last_protocol_end_time;
fprintf('Time after last protocol: %.1f s (%.1f min)\n', ...
    time_after_last_protocol, time_after_last_protocol/60);

% Save the EEG data
save('nrem_sleep_with_entrainment_40min.mat', 'EEG');

% Also save as EEGLAB .set file
pop_saveset(EEG, 'filename', 'nrem_sleep_with_entrainment_40min.set');

fprintf('Data successfully created and saved!\n');

%% Modified filter_bandpass function to improve numerical stability
function filtered_signal = filter_bandpass(signal, fs, low_cutoff, high_cutoff)
% FILTER_BANDPASS Applies a bandpass filter to a signal with improved stability
%   filtered_signal = filter_bandpass(signal, fs, low_cutoff, high_cutoff)
%   
%   Inputs:
%       signal - Input signal to be filtered
%       fs - Sampling frequency in Hz
%       low_cutoff - Lower cutoff frequency in Hz
%       high_cutoff - Upper cutoff frequency in Hz
%
%   Output:
%       filtered_signal - Filtered signal

% Use a more stable approach with two cascaded filters
% First design a lowpass filter
order_low = 3; % Lower order for better stability
[b_low, a_low] = butter(order_low, high_cutoff/(fs/2), 'low');

% Then design a highpass filter
order_high = 3; % Lower order for better stability
[b_high, a_high] = butter(order_high, low_cutoff/(fs/2), 'high');

% Apply filters in sequence with error handling
try
    % First apply lowpass
    temp_signal = filter(b_low, a_low, signal);
    
    % Then apply highpass
    filtered_signal = filter(b_high, a_high, temp_signal);
catch
    % If filtering fails, use a simpler approach for robustness
    warning('Standard filtering failed. Using alternative approach.');
    
    % Try a simple FFT-based filter as fallback
    L = length(signal);
    NFFT = 2^nextpow2(L);
    f = fs/2*linspace(0,1,NFFT/2+1);
    
    % FFT of the signal
    signal_fft = fft(signal, NFFT);
    
    % Create a bandpass filter in frequency domain
    filt = ones(size(f));
    filt(f < low_cutoff) = 0;
    filt(f > high_cutoff) = 0;
    
    % Apply smooth transitions at cutoffs to reduce ringing
    transition_width = 0.5; % Hz
    for i = 1:length(f)
        if f(i) >= low_cutoff - transition_width && f(i) < low_cutoff
            filt(i) = 0.5 * (1 + sin(pi * (f(i) - low_cutoff + transition_width/2) / transition_width));
        elseif f(i) <= high_cutoff + transition_width && f(i) > high_cutoff
            filt(i) = 0.5 * (1 + sin(pi * (high_cutoff + transition_width/2 - f(i)) / transition_width));
        end
    end
    
    % Create a symmetric filter for the negative frequencies
    full_filt = [filt, fliplr(filt(2:end-1))];
    if mod(NFFT, 2) == 1
        full_filt = [filt, fliplr(filt(2:end))];
    end
    
    % Apply filter and convert back to time domain
    filtered_fft = signal_fft .* full_filt;
    filtered_full = real(ifft(filtered_fft, NFFT));
    
    % Trim to original length
    filtered_signal = filtered_full(1:L);
end
end
