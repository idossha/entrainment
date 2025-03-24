function visualizeISPCTimecourses(EEG, ispc_results, segments, config)
    fprintf('Creating ISPC time course plots...\n');
    
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
    
    % Create figure for selected channels
    figure('Name', 'ISPC Time Courses', 'Position', [100, 100, 1200, 600], 'visible', 'off');
    
    % Select a few channels to display (evenly distributed across the head)
    channel_indices = round(linspace(1, EEG.nbchan-1, 5));
    
    % Plot ISPC time courses
    hold on;
    for i = 1:length(channel_indices)
        chan_idx = channel_indices(i);
        plot(plot_time, ispc_results(chan_idx, plot_range), 'DisplayName', EEG.chanlocs(chan_idx).labels);
    end
    
    % Add segment markers
    colors = {'b', 'g', 'r', 'm', 'c', 'y'};
    
    % Identify unique stimulation IDs from segment names
    stim_ids = [];
    for s = 1:length(segment_names)
        segment = segment_names{s};
        parts = strsplit(segment, '_');
        if length(parts) > 2
            stim_id = str2double(parts{end});
            if ~ismember(stim_id, stim_ids) && ~isnan(stim_id)
                stim_ids = [stim_ids, stim_id];
            end
        end
    end
    
    % If no specific stim IDs found, assume single stimulation
    if isempty(stim_ids)
        stim_ids = 1;
    end
    
    % Add vertical lines and labels for each stimulation
    for stim_id = stim_ids
        % Get segments for this stimulation
        if length(stim_ids) > 1
            suffix = ['_' num2str(stim_id)];
            segments_for_stim = fieldnames(segments);
            segments_for_stim = segments_for_stim(endsWith(segments_for_stim, suffix));
        else
            suffix = '';
            segments_for_stim = fieldnames(segments);
        end
        
        for s = 1:length(segments_for_stim)
            segment = segments_for_stim{s};
            segment_range = segments.(segment);
            
            % Convert to relative time
            segment_start_rel = rel_time(segment_range(1));
            segment_end_rel = rel_time(segment_range(2));
            
            % Only add markers if segment is within the plot range
            if segment_start_rel >= plot_time(1) && segment_start_rel <= plot_time(end)
                % Add vertical lines
                color_idx = mod(s-1, length(colors))+1;
                xline(segment_start_rel, '--', 'Color', colors{color_idx});
                
                % Add label for segment start
                segment_base = strrep(segment, suffix, '');
                segment_label = strrep(segment_base, '_', ' ');
                if length(stim_ids) > 1
                    segment_label = [segment_label ' ' num2str(stim_id)];
                end
                
                text(segment_start_rel, 0.95, segment_label, ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
                    'Color', colors{color_idx}, 'FontSize', 8);
            end
            
            % Add end marker if within plot range
            if segment_end_rel >= plot_time(1) && segment_end_rel <= plot_time(end)
                xline(segment_end_rel, '--', 'Color', colors{color_idx});
            end
        end
    end
    
    % Add labels and legend
    title('ISPC Time Courses for Selected Channels');
    xlabel('Time (s from protocol start)');
    ylabel('ISPC');
    legend('Location', 'eastoutside');
    ylim([0 1]);
    grid on;
    
    % Ensure x-axis shows only protocol time
    xlim([0 plot_time(end) - plot_time(1)]);
    
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
        save_path = fullfile(config.results_dir, 'ispc_timecourses.png');
        saveas(gcf, save_path);
        fprintf('Saved figure to: %s\n', save_path);
    catch err
        warning('Failed to save figure: %s', err.message);
    end
    
    % Create a figure for average ISPC across regions
    figure('Name', 'Regional Average ISPC', 'Position', [100, 100, 1200, 600], 'visible', 'off');
    
    % Calculate and plot overall average across all channels
    avg_ispc = mean(ispc_results(1:EEG.nbchan-1, plot_range), 1);
    h1 = plot(plot_time, avg_ispc, 'k-', 'LineWidth', 2, 'DisplayName', 'All Channels');
    hold on;
    
    % Calculate and plot regional averages if defined
    if isfield(config, 'regions') && ~isempty(fieldnames(config.regions))
        reg_colors = {'b', 'r', 'g', 'm', 'c'};
        color_idx = 1;
        
        % Process all regions found in the config
        region_names = fieldnames(config.regions);
        for r = 1:length(region_names)
            region_name = region_names{r};
            
            % Find channel indices for this region
            region_channels = config.regions.(region_name);
            region_indices = findChannelIndices(EEG, region_channels);
            
            % Format the region name for display
            display_name = strrep(region_name, '_', ' ');
            display_name = [upper(display_name(1)) display_name(2:end)];
            
            fprintf('Found %d channels for %s region\n', length(region_indices), display_name);
            
            if ~isempty(region_indices)
                % Calculate average ISPC for this region
                region_avg = mean(ispc_results(region_indices, plot_range), 1);
                
                % Plot with cycling through colors
                current_color = reg_colors{mod(color_idx-1, length(reg_colors))+1};
                plot(plot_time, region_avg, '-', 'Color', current_color, ...
                    'LineWidth', 1.5, 'DisplayName', display_name);
                color_idx = color_idx + 1;
            end
        end
    else
        fprintf('No regions defined or empty regions structure. Skipping regional averages.\n');
    end
    
    % Add segment markers again for regional average plot
    for stim_id = stim_ids
        % Get segments for this stimulation
        if length(stim_ids) > 1
            suffix = ['_' num2str(stim_id)];
            segments_for_stim = fieldnames(segments);
            segments_for_stim = segments_for_stim(endsWith(segments_for_stim, suffix));
        else
            suffix = '';
            segments_for_stim = fieldnames(segments);
        end
        
        for s = 1:length(segments_for_stim)
            segment = segments_for_stim{s};
            segment_range = segments.(segment);
            
            % Convert to relative time
            segment_start_rel = rel_time(segment_range(1));
            segment_end_rel = rel_time(segment_range(2));
            
            % Only add markers if segment is within the plot range
            if segment_start_rel >= plot_time(1) && segment_start_rel <= plot_time(end)
                % Add vertical lines
                color_idx = mod(s-1, length(colors))+1;
                xline(segment_start_rel, '--', 'Color', colors{color_idx});
                
                % Add label for segment start
                segment_base = strrep(segment, suffix, '');
                segment_label = strrep(segment_base, '_', ' ');
                if length(stim_ids) > 1
                    segment_label = [segment_label ' ' num2str(stim_id)];
                end
                
                text(segment_start_rel, 0.95, segment_label, ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
                    'Color', colors{color_idx}, 'FontSize', 8);
            end
            
            % Add end marker if within plot range
            if segment_end_rel >= plot_time(1) && segment_end_rel <= plot_time(end)
                xline(segment_end_rel, '--', 'Color', colors{color_idx});
            end
        end
    end
    
    title('Regional Average ISPC');
    xlabel('Time (s from protocol start)');
    ylabel('ISPC');
    legend('Location', 'eastoutside');
    ylim([0 1]);
    xlim([0 plot_time(end) - plot_time(1)]);
    grid on;
    
    % Save the regional average figure
    try
        save_path = fullfile(config.results_dir, 'regional_average_ispc.png');
        saveas(gcf, save_path);
        fprintf('Saved figure to: %s\n', save_path);
    catch err
        warning('Failed to save figure: %s', err.message);
    end
    
    %% NEW ADDITION: Top 5% Channels Analysis
    fprintf('Analyzing top 5%% channels with greatest ISPC increase...\n');
    
    % Find pre-stim and late-stim segments for the first stimulation (or without suffix)
    pre_stim_seg = '';
    late_stim_seg = '';
    
    % First try to find segments without numerical suffix
    for i = 1:length(segment_names)
        segment = segment_names{i};
        if strcmp(segment, 'pre_stim')
            pre_stim_seg = segment;
        elseif strcmp(segment, 'late_stim')
            late_stim_seg = segment;
        end
    end
    
    % If not found, look for segments with suffix _1
    if isempty(pre_stim_seg) || isempty(late_stim_seg)
        for i = 1:length(segment_names)
            segment = segment_names{i};
            if strcmp(segment, 'pre_stim_1')
                pre_stim_seg = segment;
            elseif strcmp(segment, 'late_stim_1')
                late_stim_seg = segment;
            end
        end
    end
    
    % If still not found, use the first segment that contains pre_stim or late_stim
    if isempty(pre_stim_seg)
        for i = 1:length(segment_names)
            if contains(segment_names{i}, 'pre_stim')
                pre_stim_seg = segment_names{i};
                break;
            end
        end
    end
    
    if isempty(late_stim_seg)
        for i = 1:length(segment_names)
            if contains(segment_names{i}, 'late_stim')
                late_stim_seg = segment_names{i};
                break;
            end
        end
    end
    
    % Ensure we found the segments
    if isempty(pre_stim_seg) || isempty(late_stim_seg)
        warning('Could not find pre_stim or late_stim segments. Skipping top 5%% analysis.');
    else
        fprintf('Using segments %s and %s for top 5%% analysis\n', pre_stim_seg, late_stim_seg);
        
        % Calculate average ISPC for each segment
        pre_stim_range = segments.(pre_stim_seg);
        late_stim_range = segments.(late_stim_seg);
        
        pre_stim_ispc = mean(ispc_results(:, pre_stim_range(1):pre_stim_range(2)), 2);
        late_stim_ispc = mean(ispc_results(:, late_stim_range(1):late_stim_range(2)), 2);
        
        % Calculate ISPC difference for each channel
        ispc_diff = late_stim_ispc(1:EEG.nbchan-1) - pre_stim_ispc(1:EEG.nbchan-1);
        
        % Find top 5% of channels with highest increase
        num_channels = EEG.nbchan - 1; % Exclude stimulus channel
        num_top_channels = ceil(0.05 * num_channels); % Top 5%
        
        [~, sorted_indices] = sort(ispc_diff, 'descend');
        top_channel_indices = sorted_indices(1:num_top_channels);
        
        % Get channel labels for top channels
        top_channel_labels = cell(num_top_channels, 1);
        for i = 1:num_top_channels
            chan_idx = top_channel_indices(i);
            top_channel_labels{i} = EEG.chanlocs(chan_idx).labels;
        end
        
        % List the top channels
        fprintf('\nTop %d channels with highest ISPC increase (late-stim vs pre-stim):\n', num_top_channels);
        for i = 1:num_top_channels
            fprintf('  %d. %s (increase: %.3f)\n', i, top_channel_labels{i}, ispc_diff(top_channel_indices(i)));
        end
        
        % Create a figure for top channels time course
        figure('Name', 'Top 5% Channels ISPC Time Course', 'Position', [100, 100, 1200, 600], 'visible', 'off');
        
        % Calculate average ISPC for top channels
        top_channels_avg = mean(ispc_results(top_channel_indices, plot_range), 1);
        
        % Plot the average
        plot(plot_time, top_channels_avg, 'r-', 'LineWidth', 2, 'DisplayName', 'Top 5% Channels');
        hold on;
        
        % Add overall average for comparison
        plot(plot_time, avg_ispc, 'k-', 'LineWidth', 1.5, 'DisplayName', 'All Channels');
        
        % Add segment markers
        for stim_id = stim_ids
            % Get segments for this stimulation
            if length(stim_ids) > 1
                suffix = ['_' num2str(stim_id)];
                segments_for_stim = fieldnames(segments);
                segments_for_stim = segments_for_stim(endsWith(segments_for_stim, suffix));
            else
                suffix = '';
                segments_for_stim = fieldnames(segments);
            end
            
            for s = 1:length(segments_for_stim)
                segment = segments_for_stim{s};
                segment_range = segments.(segment);
                
                % Convert to relative time
                segment_start_rel = rel_time(segment_range(1));
                segment_end_rel = rel_time(segment_range(2));
                
                % Only add markers if segment is within the plot range
                if segment_start_rel >= plot_time(1) && segment_start_rel <= plot_time(end)
                    % Add vertical lines
                    color_idx = mod(s-1, length(colors))+1;
                    xline(segment_start_rel, '--', 'Color', colors{color_idx});
                    
                    % Add label for segment start
                    segment_base = strrep(segment, suffix, '');
                    segment_label = strrep(segment_base, '_', ' ');
                    if length(stim_ids) > 1
                        segment_label = [segment_label ' ' num2str(stim_id)];
                    end
                    
                    text(segment_start_rel, 0.95, segment_label, ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
                        'Color', colors{color_idx}, 'FontSize', 8);
                end
                
                % Add end marker if within plot range
                if segment_end_rel >= plot_time(1) && segment_end_rel <= plot_time(end)
                    xline(segment_end_rel, '--', 'Color', colors{color_idx});
                end
            end
        end
        
        % Add labels and legend
        title('ISPC Time Course for Top 5% Responsive Channels');
        subtitle(['Channels with greatest ISPC increase between ' strrep(pre_stim_seg, '_', ' ') ' and ' strrep(late_stim_seg, '_', ' ')]);
        xlabel('Time (s from protocol start)');
        ylabel('ISPC');
        legend('Location', 'eastoutside');
        ylim([0 1]);
        xlim([0 plot_time(end) - plot_time(1)]);
        grid on;
        
        % Save the top channels figure
        try
            save_path = fullfile(config.results_dir, 'top_channels_ispc.png');
            saveas(gcf, save_path);
            fprintf('Saved figure to: %s\n', save_path);
            
            % Also save the list of top channels to a text file
            channels_file = fullfile(config.results_dir, 'top_channels_list.txt');
            fid = fopen(channels_file, 'w');
            if fid > 0
                fprintf(fid, 'Top %d channels with highest ISPC increase (late-stim vs pre-stim):\n\n', num_top_channels);
                for i = 1:num_top_channels
                    fprintf(fid, '%d. %s (increase: %.3f)\n', i, top_channel_labels{i}, ispc_diff(top_channel_indices(i)));
                end
                fclose(fid);
                fprintf('Saved top channels list to: %s\n', channels_file);
            else
                warning('Failed to save top channels list.');
            end
        catch err
            warning('Failed to save top channels figure: %s', err.message);
        end
    end
    
    %% NEW ADDITION: ISPC Time Course with Transition Regions
    % Calculate half window size in samples for transition region marking
    half_window_samples = floor(config.window_size * EEG.srate / 2);

    % Create a figure for average ISPC across regions with transition shading
    figure('Name', 'ISPC Time Course with Transition Regions', 'Position', [100, 100, 1200, 600], 'visible', 'off');

    % Plot the same data as in the regional average figure
    % Calculate and plot overall average across all channels
    avg_ispc = mean(ispc_results(1:EEG.nbchan-1, plot_range), 1);
    plot(plot_time, avg_ispc, 'k-', 'LineWidth', 2, 'DisplayName', 'All Channels');
    hold on;

    % Add regional plots if you have them (reused from above)
    if isfield(config, 'regions') && ~isempty(fieldnames(config.regions))
        reg_colors = {'b', 'r', 'g', 'm', 'c'};
        color_idx = 1;
        
        % Process all regions found in the config
        region_names = fieldnames(config.regions);
        for r = 1:length(region_names)
            region_name = region_names{r};
            
            % Find channel indices for this region
            region_channels = config.regions.(region_name);
            region_indices = findChannelIndices(EEG, region_channels);
            
            % Format the region name for display
            display_name = strrep(region_name, '_', ' ');
            display_name = [upper(display_name(1)) display_name(2:end)];
            
            if ~isempty(region_indices)
                % Calculate average ISPC for this region
                region_avg = mean(ispc_results(region_indices, plot_range), 1);
                
                % Plot with cycling through colors
                current_color = reg_colors{mod(color_idx-1, length(reg_colors))+1};
                plot(plot_time, region_avg, '-', 'Color', current_color, ...
                    'LineWidth', 1.5, 'DisplayName', display_name);
                color_idx = color_idx + 1;
            end
        end
    end

    % Find all segment boundaries
    boundaries = [];
    for stim_id = stim_ids
        % Get segments for this stimulation
        if length(stim_ids) > 1
            suffix = ['_' num2str(stim_id)];
            segments_for_stim = fieldnames(segments);
            segments_for_stim = segments_for_stim(endsWith(segments_for_stim, suffix));
        else
            suffix = '';
            segments_for_stim = fieldnames(segments);
        end
        
        % Extract all start and end points
        for s = 1:length(segments_for_stim)
            segment = segments_for_stim{s};
            segment_range = segments.(segment);
            
            % Add start and end to boundaries if not already there
            if ~ismember(segment_range(1), boundaries)
                boundaries = [boundaries, segment_range(1)];
            end
            if ~ismember(segment_range(2), boundaries)
                boundaries = [boundaries, segment_range(2)];
            end
        end
    end

    % Sort boundaries
    boundaries = sort(boundaries);

    % Add protocol start/end to boundaries if they aren't already included
    if ~ismember(protocol_start, boundaries)
        boundaries = [protocol_start, boundaries];
    end
    if ~ismember(protocol_end, boundaries)
        boundaries = [boundaries, protocol_end];
    end
    boundaries = sort(boundaries);

    % Determine shading regions for transitions
    fprintf('Identifying transition regions for ISPC visualization...\n');
    for i = 1:length(boundaries)
        boundary = boundaries(i);
        
        % Skip if out of plot range
        if boundary < plot_range(1) || boundary > plot_range(end)
            continue;
        end
        
        % Calculate transition region around this boundary
        trans_start = max(plot_range(1), boundary - half_window_samples);
        trans_end = min(plot_range(end), boundary + half_window_samples);
        
        % Calculate transition region time values
        trans_time_start = rel_time(trans_start);
        trans_time_end = rel_time(trans_end);
        
        % Create shaded region
        x = [trans_time_start, trans_time_end, trans_time_end, trans_time_start];
        y = [0, 0, 1, 1];
        patch(x, y, [0.8, 0.8, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        
        % Add text label for transition
        text(rel_time(boundary), 0.05, 'Transition', 'FontSize', 8, 'Color', 'r', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        
        fprintf('Added transition region shading around sample %d (%.2f - %.2f seconds)\n', ...
            boundary, trans_time_start, trans_time_end);
    end

    % Special case: Check for beginning of recording
    if half_window_samples > plot_range(1)
        % There's a partial window at the beginning
        trans_end = min(plot_range(1) + half_window_samples, plot_range(end));
        
        % Calculate transition region time values
        trans_time_start = rel_time(plot_range(1));
        trans_time_end = rel_time(trans_end);
        
        % Create shaded region
        x = [trans_time_start, trans_time_end, trans_time_end, trans_time_start];
        y = [0, 0, 1, 1];
        patch(x, y, [0.8, 0.8, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        
        % Add text label for partial window
        text(trans_time_start + (trans_time_end - trans_time_start)/2, 0.05, 'Partial Window', ...
            'FontSize', 8, 'Color', 'r', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        
        fprintf('Added partial window shading at beginning (%.2f - %.2f seconds)\n', ...
            trans_time_start, trans_time_end);
    end

    % Special case: Check for end of recording
    if plot_range(end) + half_window_samples > EEG.pnts
        % There's a partial window at the end
        trans_start = max(plot_range(end) - half_window_samples, plot_range(1));
        
        % Calculate transition region time values
        trans_time_start = rel_time(trans_start);
        trans_time_end = rel_time(plot_range(end));
        
        % Create shaded region
        x = [trans_time_start, trans_time_end, trans_time_end, trans_time_start];
        y = [0, 0, 1, 1];
        patch(x, y, [0.8, 0.8, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        
        % Add text label for partial window
        text(trans_time_start + (trans_time_end - trans_time_start)/2, 0.05, 'Partial Window', ...
            'FontSize', 8, 'Color', 'r', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        
        fprintf('Added partial window shading at end (%.2f - %.2f seconds)\n', ...
            trans_time_start, trans_time_end);
    end

    % Add segment markers (as in the original plots)
    for stim_id = stim_ids
        % Get segments for this stimulation
        if length(stim_ids) > 1
            suffix = ['_' num2str(stim_id)];
            segments_for_stim = fieldnames(segments);
            segments_for_stim = segments_for_stim(endsWith(segments_for_stim, suffix));
        else
            suffix = '';
            segments_for_stim = fieldnames(segments);
        end
        
        for s = 1:length(segments_for_stim)
            segment = segments_for_stim{s};
            segment_range = segments.(segment);
            
            % Convert to relative time
            segment_start_rel = rel_time(segment_range(1));
            segment_end_rel = rel_time(segment_range(2));
            
            % Only add markers if segment is within the plot range
            if segment_start_rel >= plot_time(1) && segment_start_rel <= plot_time(end)
                % Add vertical lines
                color_idx = mod(s-1, length(colors))+1;
                xline(segment_start_rel, '--', 'Color', colors{color_idx});
                
                % Add label for segment start
                segment_base = strrep(segment, suffix, '');
                segment_label = strrep(segment_base, '_', ' ');
                if length(stim_ids) > 1
                    segment_label = [segment_label ' ' num2str(stim_id)];
                end
                
                text(segment_start_rel, 0.95, segment_label, ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
                    'Color', colors{color_idx}, 'FontSize', 8);
            end
            
            % Add end marker if within plot range
            if segment_end_rel >= plot_time(1) && segment_end_rel <= plot_time(end)
                xline(segment_end_rel, '--', 'Color', colors{color_idx});
            end
        end
    end

    % Add a legend explaining the shaded regions
    legend_entry = findobj(gca, 'DisplayName', 'All Channels');
    if ~isempty(legend_entry)
        % Add a fake patch for the legend
        h_patch = patch([0 0 0 0], [0 0 0 0], [0.8, 0.8, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
            'DisplayName', 'Transition Regions');
        
        % Show legend
        legend('Location', 'eastoutside');
    end

    % Add explanatory note
    dim = [0.15 0.01 0.7 0.05];
    str = ['Note: Shaded regions indicate transition areas where the ISPC calculation window ',...
           'includes data from multiple segments. The window size is ' ...
           num2str(config.window_size) ' seconds.'];
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'BackgroundColor', [1 1 1 0.7]);

    % Final figure setup
    title('ISPC Time Course with Transition Regions Highlighted');
    xlabel('Time (s from protocol start)');
    ylabel('ISPC');
    ylim([0 1]);
    xlim([0 plot_time(end) - plot_time(1)]);
    grid on;

    % Save the new figure
    try
        save_path = fullfile(config.results_dir, 'ispc_timecourse_with_transitions.png');
        saveas(gcf, save_path);
        fprintf('Saved figure with transition regions to: %s\n', save_path);
    catch err
        warning('Failed to save figure: %s', err.message);
    end
end

% Helper function to find channel indices from labels
function channel_indices = findChannelIndices(EEG, channel_labels)
    channel_indices = [];
    channel_count = 0;
    for c = 1:length(channel_labels)
        for e = 1:EEG.nbchan-1 % Excluding the stimulus channel
            if strcmpi(EEG.chanlocs(e).labels, channel_labels{c})
                channel_indices = [channel_indices, e];
                channel_count = channel_count + 1;
                break;
            end
        end
    end
    
    % Debug info if few channels are found
    if channel_count < length(channel_labels) / 2
        fprintf('Warning: Found only %d out of %d channel labels. Some might be missing in the dataset.\n', ...
            channel_count, length(channel_labels));
        % Show first 5 channel labels that were requested for debugging
        if length(channel_labels) > 5
            fprintf('First 5 requested channels: %s %s %s %s %s\n', ...
                channel_labels{1}, channel_labels{2}, channel_labels{3}, channel_labels{4}, channel_labels{5});
        end
    end
end
