%% Simple EEG Channel Plot with Event Markers
% Minimal script to plot EEG channels vs time with event markers

% Initialize and load data
clear; close all; clc;
addpath('/Users/idohaber/Documents/MATLAB/eeglab2024.0/');
addpath('/Users/idohaber/Git-Projects/entrainment/src/');
addpath('/Users/idohaber/Git-Projects/entrainment/data/');
eeglab nogui;

% Load EEG data
% EEG = pop_loadset('/Users/idohaber/Git-Projects/entrainment/data/nrem_sleep_with_entrainment_40min.set');
EEG = pop_loadset('/Users/idohaber/test_data/Strength_101_N1_forSW.set');

% Create figure
figure('Position', [50, 50, 1200, 800]);

% Parameters
channel_spacing = 1;      % Spacing between channels
amplitude_scale = 0.05;   % Scale amplitude for visibility

% Downsample to make plotting faster 
downsample_factor = 10;   % Adjust as needed for your display
time_minutes = EEG.times(1:downsample_factor:end) / 60000; % Time in minutes
data = EEG.data(:, 1:downsample_factor:end);

% Plot each channel
hold on;
for ch = 1:EEG.nbchan
    % Scale and offset data
    scaled_data = data(ch, :) * amplitude_scale + (EEG.nbchan - ch);
    
    % Plot channel (use different color for SI-ENV)
    if ch == EEG.nbchan
        plot(time_minutes, scaled_data, 'r', 'LineWidth', 1);
    else
        plot(time_minutes, scaled_data, 'k', 'LineWidth', 0.5);
    end
    
    % Add channel labels on y-axis
    text(min(time_minutes) - 0.5, (EEG.nbchan - ch), EEG.chanlocs(ch).labels, ...
        'FontSize', 8, 'HorizontalAlignment', 'right');
end

% Add event markers
for e = 1:length(EEG.event)
    % Convert event time to minutes
    event_time_min = EEG.times(round(EEG.event(e).latency)) / 60000;
    
    % Plot vertical line for event
    if strcmp(EEG.event(e).type, 'stim start')
        line([event_time_min, event_time_min], [0, EEG.nbchan+1], ...
            'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
        text(event_time_min, EEG.nbchan+0.5, 'Start', ...
            'Color', 'r', 'FontWeight', 'bold', 'Rotation', 90);
    else
        line([event_time_min, event_time_min], [0, EEG.nbchan+1], ...
            'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5);
        text(event_time_min, EEG.nbchan+0.5, 'End', ...
            'Color', 'g', 'FontWeight', 'bold', 'Rotation', 90);
    end
end

% Set plot properties
xlim([min(time_minutes), max(time_minutes)]);
ylim([0, EEG.nbchan+1]);
xlabel('Time (minutes)');
ylabel('Channels');
title('EEG Channels with Event Markers');
grid on;

% Save figure
saveas(gcf, 'simple_eeg_overview.png');
fprintf('Saved figure as simple_eeg_overview.png\n');
