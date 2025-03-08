%% 3. Cross-Frequency Coupling Analysis (Phase-Amplitude Coupling)
% This shows if the phase of 1Hz modulates amplitude of higher frequencies

% Define parameters
frequencies = 2:2:30; % Frequencies to analyze for amplitude
n_freqs = length(frequencies);
window_size = 2*fs; % 2 seconds
step_size = 0.5*fs; % 0.5 seconds
n_windows_seg1 = floor((n_samples_segment - window_size) / step_size) + 1;
n_windows_seg2 = floor((n_samples_segment - window_size) / step_size) + 1;

% Pre-allocate arrays
mi_seg1 = zeros(n_freqs, 1);
mi_seg2 = zeros(n_freqs, 1);

% Extract phase of 1Hz signal
ref_filtered_seg1 = filtfilt(b, a, EEG.data(ref_channel, 1:n_samples_segment));
ref_filtered_seg2 = filtfilt(b, a, EEG.data(ref_channel, n_samples_segment+1:end));
phase_1hz_seg1 = angle(hilbert(ref_filtered_seg1));
phase_1hz_seg2 = angle(hilbert(ref_filtered_seg2));

% Phase bins for MI calculation
n_bins = 18; % Number of phase bins
phase_bins = linspace(-pi, pi, n_bins+1);

% Choose a frontal channel for PAC analysis
ch_to_analyze = 1;

% Calculate Modulation Index (MI) for each frequency
for f = 1:n_freqs
    % Design bandpass filter for current frequency
    freq = frequencies(f);
    [b_freq, a_freq] = butter(4, [freq-1 freq+1]/(fs/2), 'bandpass');
    
    % Filter and get amplitude for segment 1
    filtered_seg1 = filtfilt(b_freq, a_freq, EEG.data(ch_to_analyze, 1:n_samples_segment));
    amp_seg1 = abs(hilbert(filtered_seg1));
    
    % Filter and get amplitude for segment 2
    filtered_seg2 = filtfilt(b_freq, a_freq, EEG.data(ch_to_analyze, n_samples_segment+1:end));
    amp_seg2 = abs(hilbert(filtered_seg2));
    
    % Calculate mean amplitude at each phase bin
    mean_amp_at_phase_seg1 = zeros(n_bins, 1);
    mean_amp_at_phase_seg2 = zeros(n_bins, 1);
    
    for bin = 1:n_bins
        bin_indices_seg1 = phase_1hz_seg1 >= phase_bins(bin) & phase_1hz_seg1 < phase_bins(bin+1);
        bin_indices_seg2 = phase_1hz_seg2 >= phase_bins(bin) & phase_1hz_seg2 < phase_bins(bin+1);
        
        mean_amp_at_phase_seg1(bin) = mean(amp_seg1(bin_indices_seg1));
        mean_amp_at_phase_seg2(bin) = mean(amp_seg2(bin_indices_seg2));
    end
    
    % Normalize amplitudes
    norm_amp_seg1 = mean_amp_at_phase_seg1 / sum(mean_amp_at_phase_seg1);
    norm_amp_seg2 = mean_amp_at_phase_seg2 / sum(mean_amp_at_phase_seg2);
    
    % Calculate MI (deviation from uniform distribution)
    uniform_dist = ones(n_bins, 1) / n_bins;
    mi_seg1(f) = kldiv(uniform_dist, norm_amp_seg1);
    mi_seg2(f) = kldiv(uniform_dist, norm_amp_seg2);
end

% Plot Cross-Frequency Coupling results
figure('Position', [100 100 900 600]);

% Plot Modulation Index across frequencies
subplot(2,2,[1,3]);
bar(frequencies, [mi_seg1, mi_seg2]);
xlabel('Frequency (Hz)');
ylabel('Modulation Index');
title('Phase-Amplitude Coupling Strength');
legend({'Segment 1 (Random)', 'Segment 2 (Entrained)'});
grid on;

% Select a frequency with high MI for detailed phase-amplitude plot
% For example, let's use 10Hz (alpha)
target_freq = 10;
[~, target_idx] = min(abs(frequencies - target_freq));

% Filter and get amplitude at target frequency
[b_freq, a_freq] = butter(4, [target_freq-1 target_freq+1]/(fs/2), 'bandpass');
filtered_seg1 = filtfilt(b_freq, a_freq, EEG.data(ch_to_analyze, 1:n_samples_segment));
amp_seg1 = abs(hilbert(filtered_seg1));
filtered_seg2 = filtfilt(b_freq, a_freq, EEG.data(ch_to_analyze, n_samples_segment+1:end));
amp_seg2 = abs(hilbert(filtered_seg2));

% Calculate mean amplitude at each phase bin for target frequency
mean_amp_at_phase_seg1 = zeros(n_bins, 1);
mean_amp_at_phase_seg2 = zeros(n_bins, 1);
phase_centers = zeros(n_bins, 1);

for bin = 1:n_bins
    bin_indices_seg1 = phase_1hz_seg1 >= phase_bins(bin) & phase_1hz_seg1 < phase_bins(bin+1);
    bin_indices_seg2 = phase_1hz_seg2 >= phase_bins(bin) & phase_1hz_seg2 < phase_bins(bin+1);
    
    mean_amp_at_phase_seg1(bin) = mean(amp_seg1(bin_indices_seg1));
    mean_amp_at_phase_seg2(bin) = mean(amp_seg2(bin_indices_seg2));
    phase_centers(bin) = (phase_bins(bin) + phase_bins(bin+1)) / 2;
end

% Normalize for better comparison
norm_amp_seg1 = mean_amp_at_phase_seg1 / max(mean_amp_at_phase_seg1);
norm_amp_seg2 = mean_amp_at_phase_seg2 / max(mean_amp_at_phase_seg2);

% Plot phase-amplitude relationship for target frequency
subplot(2,2,2);
polar(phase_centers, norm_amp_seg1);
title(sprintf('%d Hz Amp by 1 Hz Phase - Segment 1', target_freq));

subplot(2,2,4);
polar(phase_centers, norm_amp_seg2);
title(sprintf('%d Hz Amp by 1 Hz Phase - Segment 2', target_freq));

% Save the figure
saveas(gcf, 'cross_frequency_coupling.png');


%% Helper function for Kullback-Leibler divergence
function kl = kldiv(P, Q)
    % Calculates KL divergence between distributions P and Q
    % Avoid zeros in Q for log calculation
    Q(Q == 0) = 1e-10;
    P(P == 0) = 1e-10;
    
    % Calculate KL divergence
    kl = sum(P .* log2(P ./ Q));
end
