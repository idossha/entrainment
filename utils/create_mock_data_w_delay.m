%% Initialize EEGLAB
clear; close all; clc;
% Make sure EEGLAB is in your MATLAB path, otherwise you'll need to add it:
addpath('/Users/idohaber/Documents/MATLAB/eeglab2024.0/');
eeglab nogui;

%% Parameters for simulated EEG data
fs = 250;                   % Sampling frequency (Hz)
total_duration = 4*60;      % 4 minutes in seconds
segment_duration = 2*60;    % 2 minutes in seconds (for each segment)
n_channels = 64;            % Number of main EEG channels
sine_freq = 1;              % Frequency of the entrainment sine wave (Hz)

% ADDED: Define entrainment lag
lag_ms = 200;               % Lag in milliseconds
lag_samples = round(lag_ms/1000 * fs);  % Convert lag to samples

% Create separate time vectors for each segment
t1 = 0:1/fs:segment_duration-1/fs;   % First 2 minutes
t2 = 0:1/fs:segment_duration-1/fs;   % Second 2 minutes
t_full = [t1, t1(end) + t2];         % Full 4 minutes

n_samples_segment = length(t1);      % Number of samples per 2-min segment
n_samples_total = length(t_full);    % Total number of samples

%% Generate simulated NREM sleep EEG data
% Create the 64 channels of simulated NREM sleep data
eeg_data = zeros(n_channels, n_samples_total);

% Create 65th channel with 1Hz sine wave for entrainment (present throughout)
sine_wave_full = 10 * sin(2*pi*sine_freq*t_full);
sine_wave_full = sine_wave_full + 0.5*randn(size(sine_wave_full)); % Add some noise


%% PART 1: First 2 minutes - Random EEG with MINIMAL chance of entrainment
for ch = 1:n_channels
    % Generate dominant delta waves (0.5-4 Hz) - NREM characteristic
    delta = 15 * filter_bandpass(randn(1, n_samples_segment), fs, 0.5, 4);
    
    % Add some theta waves (4-8 Hz) - less power than delta
    theta = 5 * filter_bandpass(randn(1, n_samples_segment), fs, 4, 8);
    
    % Add minimal alpha waves (8-12 Hz)
    alpha = 2 * filter_bandpass(randn(1, n_samples_segment), fs, 8, 12);
    
    % Add minimal beta waves (13-30 Hz)
    beta = 1 * filter_bandpass(randn(1, n_samples_segment), fs, 13, 30);
    
    % Add minimal gamma waves (>30 Hz)
    gamma = 0.5 * filter_bandpass(randn(1, n_samples_segment), fs, 30, 45);
    
    % Combine all frequency bands - completely random with no entrainment
    eeg_data(ch, 1:n_samples_segment) = delta + theta + alpha + beta + gamma;
    
    % Add some spatial correlation between channels (nearby channels are more correlated)
    if ch > 1
        % Add correlation with neighboring channels
        neighbor_contrib = 0.3 * eeg_data(max(1, ch-1), 1:n_samples_segment) + 0.2 * eeg_data(max(1, ch-2), 1:n_samples_segment);
        eeg_data(ch, 1:n_samples_segment) = 0.5 * eeg_data(ch, 1:n_samples_segment) + neighbor_contrib;
    end
    
    % Add some random channel-specific characteristics
    eeg_data(ch, 1:n_samples_segment) = eeg_data(ch, 1:n_samples_segment) * (0.9 + 0.2*rand());
end

%% PART 2: Second 2 minutes - Clear ENTRAINMENT to the 1Hz stimulus WITH LAG
% Extract sine wave for second segment
sine_wave_segment2 = sine_wave_full(n_samples_segment+1:end);

% Create a lagged version of the sine wave for entrainment
% Circular shift to maintain the same length
lagged_sine_wave = circshift(sine_wave_segment2, lag_samples);

for ch = 1:n_channels
    % Generate baseline NREM sleep EEG components as before
    delta = 15 * filter_bandpass(randn(1, n_samples_segment), fs, 0.5, 4);
    theta = 5 * filter_bandpass(randn(1, n_samples_segment), fs, 4, 8);
    alpha = 2 * filter_bandpass(randn(1, n_samples_segment), fs, 8, 12);
    beta = 1 * filter_bandpass(randn(1, n_samples_segment), fs, 13, 30);
    gamma = 0.5 * filter_bandpass(randn(1, n_samples_segment), fs, 30, 45);
    
    % Combine baseline frequency bands
    baseline = delta + theta + alpha + beta + gamma;
    
    % KEY CHANGE: Add entrainment with LAG by using the lagged sine wave
    entrainment_strength = 0.2 + 0.6 * (1 - min(1, ch/n_channels)); % Higher for lower channel numbers
    
    % Add phase-locked component to create entrainment with the lagged sine
    entrained_signal = baseline + entrainment_strength * lagged_sine_wave;
    
    % Store the entrained signal in the second half of the data
    eeg_data(ch, n_samples_segment+1:end) = entrained_signal;
    
    % Add some spatial correlation between channels
    if ch > 1
        neighbor_contrib = 0.3 * eeg_data(max(1, ch-1), n_samples_segment+1:end) + 0.2 * eeg_data(max(1, ch-2), n_samples_segment+1:end);
        eeg_data(ch, n_samples_segment+1:end) = 0.5 * eeg_data(ch, n_samples_segment+1:end) + neighbor_contrib;
    end
    
    % Add some random channel-specific characteristics
    eeg_data(ch, n_samples_segment+1:end) = eeg_data(ch, n_samples_segment+1:end) * (0.9 + 0.2*rand());
end

% Combine all channels including the 65th sine wave channel (ORIGINAL unlagged sine)
all_data = [eeg_data; sine_wave_full];

%% Create an EEGLAB data structure
EEG = eeg_emptyset();
EEG.data = all_data;
EEG.srate = fs;
EEG.pnts = n_samples_total;
EEG.trials = 1;
EEG.xmin = 0;
EEG.xmax = total_duration;
EEG.times = t_full;
EEG.nbchan = size(all_data, 1);

% Create channel labels
for ch = 1:n_channels
    EEG.chanlocs(ch).labels = sprintf('Chan%d', ch);
end
EEG.chanlocs(n_channels+1).labels = 'SI-ENV';

% Update the EEGLAB data structure
EEG = eeg_checkset(EEG);

% Save the EEG data
save('nrem_sleep_with_lagged_entrainment.mat', 'EEG');

% Also save as EEGLAB .set file
pop_saveset(EEG, 'filename', 'nrem_sleep_with_lagged_entrainment.set');

% Output lag information to console
fprintf('Created dataset with %d ms lag (%d samples) in entrainment\n', lag_ms, lag_samples);

%% 
function filtered_signal = filter_bandpass(signal, fs, low_cutoff, high_cutoff)
% FILTER_BANDPASS Applies a bandpass filter to a signal
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

% Design a Butterworth bandpass filter
order = 4; % Filter order
[b, a] = butter(order, [low_cutoff high_cutoff]/(fs/2), 'bandpass');

% Apply filter to the signal
filtered_signal = filtfilt(b, a, signal);
end
