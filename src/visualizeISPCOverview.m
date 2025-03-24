function visualizeISPCOverview(EEG, ispc_results, phases, filtered_data, segments, stim_samples, stim_channel_idx, config)
    fprintf('Creating overview figure...\n');
    
    % Select a central channel for demonstration
    channel_labels = {EEG.chanlocs(1:end-1).labels};
    central_channels = {'Cz', 'CPz', 'FCz', 'Pz', 'Fz'};
    channel_idx = find(ismember(channel_labels, central_channels), 1);
    if isempty(channel_idx), channel_idx = 1; end
    
    % Get segment ranges to determine protocol boundaries
    segment_names = fieldnames(segments);
    
    % Find pre-stim and post-stim segments to determine protocol start/end
    pre_stim_start = [];
    post_stim_end = [];
    
    for i = 1:length(segment_names)
        segment = segment_names{i};
        if contains(segment, 'pre_stim')
            pre_stim_start = [pre_stim_start; segments.(segment)(1)];
        end
        if contains(segment, 'post_stim')
            post_stim_end = [post_stim_end; segments.(segment)(2)];
        end
    end
    
    % Ensure we have valid boundaries
    if isempty(pre_stim_start) || isempty(post_stim_end)
        error('Could not find pre-stim or post-stim segments to determine protocol boundaries');
    end
    
    % Protocol time range
    protocol_start = min(pre_stim_start);
    protocol_end = max(post_stim_end);
    
    fprintf('Protocol time range: [%d-%d] (%d samples)\n', ...
        protocol_start, protocol_end, protocol_end-protocol_start+1);
    
    % Create absolute and relative time vectors
    abs_time = (0:EEG.pnts-1) / EEG.srate;
    rel_time = abs_time - abs_time(protocol_start); % 0 = start of pre-stim
    
    % Limit the data to show only the protocol time range
    plot_range = protocol_start:protocol_end;
    plot_time = rel_time(plot_range);
    
    % Create figure
    figure('Name', 'Analysis Overview', 'Position', [100, 100, 1200, 800], 'visible', 'off');
    
    % Plot raw EEG data
    subplot(5, 1, 1);
    plot(plot_time, EEG.data(channel_idx, plot_range));
    title(['Raw EEG - Channel ' EEG.chanlocs(channel_idx).labels]);
    ylabel('Amplitude (µV)');
    xlim([0 plot_time(end) - plot_time(1)]);
    
    % Plot stimulus
    subplot(5, 1, 2);
    plot(plot_time, EEG.data(stim_channel_idx, plot_range));
    title([num2str(config.stim_freq) ' Hz Stimulus Signal']);
    ylabel('Amplitude');
    xlim([0 plot_time(end) - plot_time(1)]);
    
    % Plot filtered data
    subplot(5, 1, 3);
    plot(plot_time, filtered_data(channel_idx, plot_range));
    title(['Filtered EEG at ' num2str(config.stim_freq) ' Hz']);
    ylabel('Amplitude (µV)');
    xlim([0 plot_time(end) - plot_time(1)]);
    
    % Plot phase angles
    subplot(5, 1, 4);
    plot(plot_time, phases(channel_idx, plot_range));
    hold on;
    plot(plot_time, phases(stim_channel_idx, plot_range));
    title('Phase Angles');
    ylabel('Phase (rad)');
    legend('EEG Channel', 'Stimulus');
    ylim([-2*pi 2*pi]); % Setting y-limits to [-2π, 2π]
    xlim([0 plot_time(end) - plot_time(1)]);
    
    % Plot ISPC over time
    subplot(5, 1, 5);
    plot(plot_time, ispc_results(channel_idx, plot_range));
    title(['ISPC over Time']);
    xlabel('Time (s from protocol start)');
    ylabel('ISPC');
    ylim([0 1]);
    xlim([0 plot_time(end) - plot_time(1)]);
    
    % Add vertical lines for stim start/end
    colors = {'g', 'r', 'b', 'm', 'c', 'y'};
    stim_color_idx = 1;
    
    for stim_idx = 1:size(stim_samples, 1)
        % Convert to relative time
        stim_start_rel = rel_time(stim_samples(stim_idx, 1));
        stim_end_rel = rel_time(stim_samples(stim_idx, 2));
        
        for sp = 1:5
            subplot(5, 1, sp);
            hold on;
            
            % Only add markers if stim is within the plot range
            if stim_start_rel >= plot_time(1) && stim_start_rel <= plot_time(end)
                xline(stim_start_rel, '--', 'Color', colors{stim_color_idx}, 'LineWidth', 1.5);
                text(stim_start_rel, 0.1, ['Stim ' num2str(stim_idx) ' Start'], ...
                    'Color', colors{stim_color_idx}, 'FontSize', 8, ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                    'Rotation', 90);
            end
            
            if stim_end_rel >= plot_time(1) && stim_end_rel <= plot_time(end)
                xline(stim_end_rel, '--', 'Color', colors{stim_color_idx}, 'LineWidth', 1.5);
                text(stim_end_rel, 0.1, ['Stim ' num2str(stim_idx) ' End'], ...
                    'Color', colors{stim_color_idx}, 'FontSize', 8, ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                    'Rotation', 90);
            end
        end
        
        stim_color_idx = mod(stim_color_idx, length(colors)) + 1;
    end
    
    % Add segment boundaries
    seg_color_idx = 1;
    
    for i = 1:length(segment_names)
        segment = segment_names{i};
        segment_range = segments.(segment);
        
        % Convert to relative time
        segment_start_rel = rel_time(segment_range(1));
        segment_end_rel = rel_time(segment_range(2));
        
        if ~contains(segment, 'stim') % Skip actual stim segments to avoid clutter
            for sp = 1:5
                subplot(5, 1, sp);
                hold on;
                
                % Only add markers if segment is within the plot range
                if segment_start_rel >= plot_time(1) && segment_start_rel <= plot_time(end)
                    xline(segment_start_rel, ':', 'Color', colors{seg_color_idx}, 'Alpha', 0.5);
                end
                
                if segment_end_rel >= plot_time(1) && segment_end_rel <= plot_time(end)
                    xline(segment_end_rel, ':', 'Color', colors{seg_color_idx}, 'Alpha', 0.5);
                end
            end
            
            seg_color_idx = mod(seg_color_idx, length(colors)) + 1;
        end
    end
    
    % Create directory if it doesn't exist (safeguard)
    if ~exist(config.results_dir, 'dir')
        [success, msg] = mkdir(config.results_dir);
        if ~success
            warning('Failed to create directory %s: %s', config.results_dir, msg);
            fprintf('Will attempt to save to parent directory instead.\n');
            % Try to save to parent directory as fallback
            [parent_dir, ~, ~] = fileparts(config.results_dir);
            config.results_dir = parent_dir;
        else
            fprintf('Created directory: %s\n', config.results_dir);
        end
    end
    
    % Save the figure
    try
        save_path = fullfile(config.results_dir, 'overview_figure.png');
        saveas(gcf, save_path);
        fprintf('Saved figure to: %s\n', save_path);
    catch err
        warning('Failed to save figure: %s', err.message);
    end
end
