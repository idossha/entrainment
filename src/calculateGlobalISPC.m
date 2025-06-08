function [globalISPC, channelISPC, timeISPC, globalStd, nonStimGlobalISPC, nonStimChannelISPC, nonStimTimeISPC, nonStimGlobalStd] = calculateGlobalISPC(EEG, stim_channel_idx, config)
    % CALCULATEGLOBALISPC Calculates global and per-channel ISPC across the entire dataset
    %   Also calculates ISPC metrics excluding stimulation periods
    %
    % Inputs:
    %   EEG - EEG data structure
    %   stim_channel_idx - Index of the stimulus channel
    %   config - Configuration structure
    %
    % Outputs:
    %   globalISPC - Global mean ISPC across all channels and time
    %   channelISPC - ISPC value for each channel (averaged across time)
    %   timeISPC - ISPC value for each time point (averaged across channels)
    %   globalStd - Global standard deviation of ISPC
    %   nonStimGlobalISPC - Global mean ISPC excluding stimulation periods
    %   nonStimChannelISPC - Channel ISPC values excluding stimulation periods
    %   nonStimTimeISPC - Time point ISPC values excluding stimulation periods
    %   nonStimGlobalStd - Global standard deviation excluding stimulation periods
    
    fprintf('Calculating global and per-channel ISPC across entire dataset...\n');
    
    % Use the entire dataset
    time_range = [1, EEG.pnts];
    
    % Calculate ISPC for the full dataset
    [ispc_results, ~, ~] = calculateISPC(EEG, stim_channel_idx, time_range, config);
    
    % Calculate per-channel mean ISPC (spatial distribution)
    channelISPC = mean(ispc_results, 2, 'omitnan');
    
    % Calculate per-time mean ISPC (temporal distribution)
    timeISPC = mean(ispc_results, 1, 'omitnan')';
    
    % Calculate global mean and standard deviation
    globalISPC = mean(channelISPC, 'omitnan'); % Should be equal to mean(timeISPC)
    globalStd = std(channelISPC, 'omitnan');
    
    fprintf('Global ISPC mean: %.4f, std: %.4f\n', globalISPC, globalStd);
    
    % Find stimulation events to create the exclusion mask
    try
        % Get marker settings from config
        start_marker = config.protocol.start_marker;
        end_marker = config.protocol.end_marker;
        
        % Find event indices for start and end markers
        start_marker_events = strcmpi({EEG.event.type}, start_marker);
        end_marker_events = strcmpi({EEG.event.type}, end_marker);
        
        if ~any(start_marker_events) || ~any(end_marker_events)
            warning('Start or end markers not found: %s/%s. Cannot calculate non-stim ISPC.', start_marker, end_marker);
            nonStimGlobalISPC = NaN;
            nonStimChannelISPC = NaN(size(channelISPC));
            nonStimTimeISPC = NaN(size(timeISPC));
            nonStimGlobalStd = NaN;
        else
            % Get all occurrences of start and end events
            start_indices = find(start_marker_events);
            end_indices = find(end_marker_events);
            
            % Get the sample indices
            start_samples = round([EEG.event(start_indices).latency]);
            end_samples = round([EEG.event(end_indices).latency]);
            
            % Create a mask of times to include (1 = include, 0 = exclude)
            include_mask = ones(1, EEG.pnts);
            
            % Set up post-stim exclusion period in samples
            post_stim_exclusion_samples = round(config.post_stim_exclusion * EEG.srate);
            
            % Match starts and ends - simple pairing approach
            for i = 1:length(start_indices)
                start_sample = start_samples(i);
                
                % Find the next end marker after this start
                next_end_indices = find(end_samples > start_sample, 1, 'first');
                
                if ~isempty(next_end_indices)
                    end_sample = end_samples(next_end_indices);
                    
                    % Exclude samples from start through end plus post-stim period
                    exclusion_end = min(EEG.pnts, end_sample + post_stim_exclusion_samples);
                    include_mask(start_sample:exclusion_end) = 0;
                end
            end
            
            % Apply the mask to the ISPC results
            non_stim_samples = find(include_mask);
            fprintf('Excluding %d stimulation samples (%.2f%% of data) for non-stim ISPC calculation\n', ...
                EEG.pnts - length(non_stim_samples), (EEG.pnts - length(non_stim_samples))/EEG.pnts*100);
            
            % Calculate non-stim ISPC metrics
            non_stim_ispc = ispc_results(:, non_stim_samples);
            
            % Calculate non-stim channel ISPC (spatial distribution)
            nonStimChannelISPC = mean(non_stim_ispc, 2, 'omitnan');
            
            % Calculate non-stim time ISPC (temporal distribution)
            % Note: this will only include the non-stimulation time points
            nonStimTimeISPC = mean(non_stim_ispc, 1, 'omitnan')';
            
            % Calculate non-stim global metrics
            nonStimGlobalISPC = mean(nonStimChannelISPC, 'omitnan');
            nonStimGlobalStd = std(nonStimChannelISPC, 'omitnan');
            
            fprintf('Non-stim global ISPC mean: %.4f, std: %.4f\n', nonStimGlobalISPC, nonStimGlobalStd);
        end
    catch err
        warning('Error calculating non-stim ISPC: %s', err.message);
        nonStimGlobalISPC = NaN;
        nonStimChannelISPC = NaN(size(channelISPC));
        nonStimTimeISPC = NaN(size(1, length(find(include_mask))));
        nonStimGlobalStd = NaN;
    end
    
    % Get information for filename
    [~, setname, ~] = fileparts(config.eeg_filename);
    
    % Create a unique filename based on the EEG file being processed
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    
    % Generate the output filename for the MAT file
    mat_filename = fullfile(config.results_dir, sprintf('ispc_analysis_%s_%s.mat', setname, timestamp));
    
    % Save ISPC values to a MAT file
    save(mat_filename, 'channelISPC', 'timeISPC', 'globalISPC', 'globalStd', ...
        'nonStimGlobalISPC', 'nonStimChannelISPC', 'nonStimTimeISPC', 'nonStimGlobalStd', 'config');
    fprintf('ISPC values saved to: %s\n', mat_filename);
    
    % 1. Create histogram for spatial distribution (across channels)
    figure('Position', [100, 100, 800, 600]);
    
    % Create histogram
    histogram(channelISPC, 20, 'Normalization', 'count');
    
    % Add title and labels
    title_str = sprintf('ISPC Channel Distribution - %s', setname);
    title(title_str, 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('ISPC Value', 'FontSize', 12);
    ylabel('Count (Number of Channels)', 'FontSize', 12);
    
    % Add descriptive statistics as text
    stats_text = sprintf('Channels: %d\nMean = %.4f\nMedian = %.4f\nStd = %.4f\nMin = %.4f\nMax = %.4f', ...
        length(channelISPC), ...
        mean(channelISPC, 'omitnan'), ...
        median(channelISPC, 'omitnan'), ...
        std(channelISPC, 'omitnan'), ...
        min(channelISPC), ...
        max(channelISPC));
    
    % Add text box with statistics
    annotation('textbox', [0.65, 0.6, 0.3, 0.3], ...
        'String', stats_text, ...
        'EdgeColor', 'none', ...
        'FitBoxToText', 'on', ...
        'BackgroundColor', [0.95, 0.95, 0.95], ...
        'FontSize', 11);
    
    % Add a vertical line at the mean
    hold on;
    ylims = ylim;
    plot([mean(channelISPC, 'omitnan'), mean(channelISPC, 'omitnan')], [0, ylims(2)], 'r--', 'LineWidth', 2);
    legend('ISPC Distribution', 'Mean ISPC', 'Location', 'northwest');
    
    grid on;
    box on;
    
    % Generate the output filename for the spatial distribution figure
    fig_spatial_filename = fullfile(config.results_dir, sprintf('ISPC_channel_distribution_%s_%s.png', setname, timestamp));
    
    % Save the figure
    saveas(gcf, fig_spatial_filename);
    fprintf('Channel distribution figure saved to: %s\n', fig_spatial_filename);
    
    % Save as a MATLAB figure as well for later editing
    fig_spatial_mat_filename = fullfile(config.results_dir, sprintf('ISPC_channel_distribution_%s_%s.fig', setname, timestamp));
    saveas(gcf, fig_spatial_mat_filename);
    
    % 2. Create histogram for temporal distribution (across time)
    figure('Position', [100, 100, 800, 600]);
    
    % Create histogram
    histogram(timeISPC, 20, 'Normalization', 'count');
    
    % Add title and labels
    title_str = sprintf('ISPC Temporal Distribution - %s', setname);
    title(title_str, 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('ISPC Value', 'FontSize', 12);
    ylabel('Count (Number of Time Points)', 'FontSize', 12);
    
    % Add descriptive statistics as text
    stats_text = sprintf('Time Points: %d\nMean = %.4f\nMedian = %.4f\nStd = %.4f\nMin = %.4f\nMax = %.4f', ...
        length(timeISPC), ...
        mean(timeISPC, 'omitnan'), ...
        median(timeISPC, 'omitnan'), ...
        std(timeISPC, 'omitnan'), ...
        min(timeISPC), ...
        max(timeISPC));
    
    % Add text box with statistics
    annotation('textbox', [0.65, 0.6, 0.3, 0.3], ...
        'String', stats_text, ...
        'EdgeColor', 'none', ...
        'FitBoxToText', 'on', ...
        'BackgroundColor', [0.95, 0.95, 0.95], ...
        'FontSize', 11);
    
    % Add a vertical line at the mean
    hold on;
    ylims = ylim;
    plot([mean(timeISPC, 'omitnan'), mean(timeISPC, 'omitnan')], [0, ylims(2)], 'r--', 'LineWidth', 2);
    legend('ISPC Distribution', 'Mean ISPC', 'Location', 'northwest');
    
    grid on;
    box on;
    
    % Generate the output filename for the temporal distribution figure
    fig_temporal_filename = fullfile(config.results_dir, sprintf('ISPC_time_distribution_%s_%s.png', setname, timestamp));
    
    % Save the figure
    saveas(gcf, fig_temporal_filename);
    fprintf('Temporal distribution figure saved to: %s\n', fig_temporal_filename);
    
    % Save as a MATLAB figure as well for later editing
    fig_temporal_mat_filename = fullfile(config.results_dir, sprintf('ISPC_time_distribution_%s_%s.fig', setname, timestamp));
    saveas(gcf, fig_temporal_mat_filename);
    
    % Only proceed with non-stim visualizations if data is available
    if ~isnan(nonStimGlobalISPC)
        % 3. Create histogram for non-stim channel ISPC (spatial distribution)
        figure('Position', [100, 100, 800, 600]);
        
        % Create histogram
        histogram(nonStimChannelISPC, 20, 'Normalization', 'count');
        
        % Add title and labels
        title_str = sprintf('Non-Stim ISPC Channel Distribution - %s', setname);
        title(title_str, 'FontSize', 14, 'FontWeight', 'bold');
        xlabel('ISPC Value', 'FontSize', 12);
        ylabel('Count (Number of Channels)', 'FontSize', 12);
        
        % Add descriptive statistics as text
        stats_text = sprintf('Channels: %d\nMean = %.4f\nMedian = %.4f\nStd = %.4f\nMin = %.4f\nMax = %.4f\nPost-stim exclusion: %d s', ...
            length(nonStimChannelISPC), ...
            nonStimGlobalISPC, ...
            median(nonStimChannelISPC, 'omitnan'), ...
            nonStimGlobalStd, ...
            min(nonStimChannelISPC), ...
            max(nonStimChannelISPC), ...
            config.post_stim_exclusion);
        
        % Add text box with statistics
        annotation('textbox', [0.65, 0.6, 0.3, 0.3], ...
            'String', stats_text, ...
            'EdgeColor', 'none', ...
            'FitBoxToText', 'on', ...
            'BackgroundColor', [0.95, 0.95, 0.95], ...
            'FontSize', 11);
        
        % Add a vertical line at the mean
        hold on;
        ylims = ylim;
        plot([nonStimGlobalISPC, nonStimGlobalISPC], [0, ylims(2)], 'r--', 'LineWidth', 2);
        legend('ISPC Distribution', 'Mean ISPC', 'Location', 'northwest');
        
        grid on;
        box on;
        
        % Generate the output filename for the non-stim spatial distribution figure
        fig_nonstim_spatial_filename = fullfile(config.results_dir, sprintf('NonStim_ISPC_channel_distribution_%s_%s.png', setname, timestamp));
        
        % Save the figure
        saveas(gcf, fig_nonstim_spatial_filename);
        fprintf('Non-stim channel distribution figure saved to: %s\n', fig_nonstim_spatial_filename);
        
        % Save as a MATLAB figure as well for later editing
        fig_nonstim_spatial_mat_filename = fullfile(config.results_dir, sprintf('NonStim_ISPC_channel_distribution_%s_%s.fig', setname, timestamp));
        saveas(gcf, fig_nonstim_spatial_mat_filename);
        
        % 4. Create histogram for non-stim temporal ISPC (temporal distribution)
        figure('Position', [100, 100, 800, 600]);
        
        % Create histogram
        histogram(nonStimTimeISPC, 20, 'Normalization', 'count');
        
        % Add title and labels
        title_str = sprintf('Non-Stim ISPC Temporal Distribution - %s', setname);
        title(title_str, 'FontSize', 14, 'FontWeight', 'bold');
        xlabel('ISPC Value', 'FontSize', 12);
        ylabel('Count (Number of Time Points)', 'FontSize', 12);
        
        % Add descriptive statistics as text
        stats_text = sprintf('Time Points: %d\nMean = %.4f\nMedian = %.4f\nStd = %.4f\nMin = %.4f\nMax = %.4f\nPost-stim exclusion: %d s', ...
            length(nonStimTimeISPC), ...
            mean(nonStimTimeISPC, 'omitnan'), ...
            median(nonStimTimeISPC, 'omitnan'), ...
            std(nonStimTimeISPC, 'omitnan'), ...
            min(nonStimTimeISPC), ...
            max(nonStimTimeISPC), ...
            config.post_stim_exclusion);
        
        % Add text box with statistics
        annotation('textbox', [0.65, 0.6, 0.3, 0.3], ...
            'String', stats_text, ...
            'EdgeColor', 'none', ...
            'FitBoxToText', 'on', ...
            'BackgroundColor', [0.95, 0.95, 0.95], ...
            'FontSize', 11);
        
        % Add a vertical line at the mean
        hold on;
        ylims = ylim;
        plot([mean(nonStimTimeISPC, 'omitnan'), mean(nonStimTimeISPC, 'omitnan')], [0, ylims(2)], 'r--', 'LineWidth', 2);
        legend('ISPC Distribution', 'Mean ISPC', 'Location', 'northwest');
        
        grid on;
        box on;
        
        % Generate the output filename for the non-stim temporal distribution figure
        fig_nonstim_temporal_filename = fullfile(config.results_dir, sprintf('NonStim_ISPC_time_distribution_%s_%s.png', setname, timestamp));
        
        % Save the figure
        saveas(gcf, fig_nonstim_temporal_filename);
        fprintf('Non-stim temporal distribution figure saved to: %s\n', fig_nonstim_temporal_filename);
        
        % Save as a MATLAB figure as well for later editing
        fig_nonstim_temporal_mat_filename = fullfile(config.results_dir, sprintf('NonStim_ISPC_time_distribution_%s_%s.fig', setname, timestamp));
        saveas(gcf, fig_nonstim_temporal_mat_filename);
        
        % 5. Compare histograms for channel distribution (with/without stimulation)
        figure('Position', [100, 100, 800, 600]);
        
        % Create histograms
        hold on;
        histogram(channelISPC, 20, 'Normalization', 'count', 'FaceColor', [0.3 0.5 0.8], 'FaceAlpha', 0.6);
        histogram(nonStimChannelISPC, 20, 'Normalization', 'count', 'FaceColor', [0.8 0.3 0.3], 'FaceAlpha', 0.6);
        
        % Add title and labels
        title_str = sprintf('Channel ISPC Comparison (With vs. Without Stimulation) - %s', setname);
        title(title_str, 'FontSize', 14, 'FontWeight', 'bold');
        xlabel('ISPC Value', 'FontSize', 12);
        ylabel('Count (Number of Channels)', 'FontSize', 12);
        
        % Add descriptive statistics as text
        stats_text = sprintf(['Full Dataset:\n  Mean = %.4f\n  Std = %.4f\n\n' ...
            'Non-Stim Only:\n  Mean = %.4f\n  Std = %.4f\n\n' ...
            'Difference:\n  Mean diff = %.4f\n  %% Change = %.2f%%'], ...
            globalISPC, globalStd, ...
            nonStimGlobalISPC, nonStimGlobalStd, ...
            nonStimGlobalISPC - globalISPC, ...
            100 * (nonStimGlobalISPC - globalISPC) / globalISPC);
        
        % Add text box with statistics
        annotation('textbox', [0.65, 0.6, 0.3, 0.3], ...
            'String', stats_text, ...
            'EdgeColor', 'none', ...
            'FitBoxToText', 'on', ...
            'BackgroundColor', [0.95, 0.95, 0.95], ...
            'FontSize', 11);
        
        % Add vertical lines at the means
        ylims = ylim;
        plot([globalISPC, globalISPC], [0, ylims(2)], 'b--', 'LineWidth', 2);
        plot([nonStimGlobalISPC, nonStimGlobalISPC], [0, ylims(2)], 'r--', 'LineWidth', 2);
        
        legend('Full Dataset', 'Non-Stim Only', 'Full Mean', 'Non-Stim Mean', 'Location', 'northwest');
        
        grid on;
        box on;
        
        % Generate the output filename for the comparison figure
        fig_compare_filename = fullfile(config.results_dir, sprintf('ISPC_channel_comparison_%s_%s.png', setname, timestamp));
        
        % Save the figure
        saveas(gcf, fig_compare_filename);
        fprintf('Channel comparison figure saved to: %s\n', fig_compare_filename);
        
        % Save as a MATLAB figure as well for later editing
        fig_compare_mat_filename = fullfile(config.results_dir, sprintf('ISPC_channel_comparison_%s_%s.fig', setname, timestamp));
        saveas(gcf, fig_compare_mat_filename);
    end
end
