function visualizeNormalizedISPC(EEG, ispc_results, norm_ispc_results, z_scored_ispc, segments, config)
    % VISUALIZENORMALIZEDISPC Creates visualizations for normalized ISPC values
    %
    % Inputs:
    %   EEG - EEG data structure
    %   ispc_results - Original ISPC values
    %   norm_ispc_results - Normalized ISPC values (as % of global mean)
    %   z_scored_ispc - Z-scored ISPC values for each channel
    %   segments - Structure with segment definitions
    %   config - Configuration structure
    
    fprintf('Creating visualizations for normalized ISPC...\n');
    
    % Get segment names and find pre-stim and late-stim
    segment_names = fieldnames(segments);
    pre_stim_seg = '';
    late_stim_seg = '';
    
    for i = 1:length(segment_names)
        if strcmp(segment_names{i}, 'pre_stim')
            pre_stim_seg = segment_names{i};
        elseif strcmp(segment_names{i}, 'late_stim')
            late_stim_seg = segment_names{i};
        end
    end
    
    if isempty(pre_stim_seg) || isempty(late_stim_seg)
        warning('Could not find pre_stim or late_stim segments. Some visualizations may be skipped.');
    end
    
    %% 1. Create time course visualization
    % Get segment ranges to determine protocol boundaries
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
    
    % Protocol time range
    protocol_start = min(pre_stim_start);
    protocol_end = max(post_stim_end);
    
    % Create time vectors
    abs_time = (0:EEG.pnts-1) / EEG.srate;
    rel_time = abs_time - abs_time(protocol_start);
    
    % Limit the data to show only the protocol time range
    plot_range = protocol_start:protocol_end;
    plot_time = rel_time(plot_range);
    
    % Create figure for normalized time courses
    figure('Name', 'Normalized ISPC Time Courses', 'Position', [100, 100, 1200, 600], 'visible', 'off');
    
    % Calculate average normalized ISPC across all channels
    avg_norm_ispc = mean(norm_ispc_results(1:EEG.nbchan-1, plot_range), 1, 'omitnan');
    
    % Calculate standard deviation across channels
    std_norm_ispc = std(norm_ispc_results(1:EEG.nbchan-1, plot_range), 0, 1, 'omitnan');
    
    % Plot mean line
    h1 = plot(plot_time, avg_norm_ispc, 'k-', 'LineWidth', 2, 'DisplayName', 'All Channels (Mean)');
    hold on;
    
    % Find min and max values for y-axis
    y_max = max(avg_norm_ispc + std_norm_ispc) * 1.1;  % Add 10% margin
    y_min = min(avg_norm_ispc - std_norm_ispc) * 0.9;  % Add 10% margin
    y_min = max(0, y_min);  % Ensure y_min is not negative
    
    % Plot regional averages if available
    if isfield(config, 'regions') && ~isempty(fieldnames(config.regions))
        reg_colors = {'b', 'r', 'g', 'm', 'c'};
        color_idx = 1;
        
        region_names = fieldnames(config.regions);
        for r = 1:length(region_names)
            region_name = region_names{r};
            region_channels = config.regions.(region_name);
            
            % Find channel indices
            region_indices = [];
            for c = 1:length(region_channels)
                for e = 1:EEG.nbchan-1
                    if strcmpi(EEG.chanlocs(e).labels, region_channels{c})
                        region_indices = [region_indices, e];
                        break;
                    end
                end
            end
            
            if ~isempty(region_indices)
                % Calculate average normalized ISPC for this region
                region_avg = mean(norm_ispc_results(region_indices, plot_range), 1, 'omitnan');
                
                % Plot with cycling through colors
                display_name = strrep(region_name, '_', ' ');
                display_name = [upper(display_name(1)) display_name(2:end) ' (Mean)'];
                
                current_color = reg_colors{mod(color_idx-1, length(reg_colors))+1};
                plot(plot_time, region_avg, '-', 'Color', current_color, ...
                    'LineWidth', 1.5, 'DisplayName', display_name);
                
                % Update y-axis limits if needed
                y_max = max(y_max, max(region_avg) * 1.1);
                y_min = min(y_min, max(0, min(region_avg) * 0.9));
                
                color_idx = color_idx + 1;
            end
        end
    end
    
    % Add segment markers
    colors = {'b', 'g', 'r', 'm', 'c', 'y'};
    
    for s = 1:length(segment_names)
        segment = segment_names{s};
        segment_range = segments.(segment);
        
        % Convert to relative time
        segment_start_rel = rel_time(segment_range(1));
        segment_end_rel = rel_time(segment_range(2));
        
        if segment_start_rel >= plot_time(1) && segment_start_rel <= plot_time(end)
            color_idx = mod(s-1, length(colors))+1;
            xline(segment_start_rel, '--', 'Color', colors{color_idx});
            
            % Add label
            segment_label = strrep(segment, '_', ' ');
            text(segment_start_rel, 0.95*y_max, segment_label, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
                'Color', colors{color_idx}, 'FontSize', 8);
        end
        
        if segment_end_rel >= plot_time(1) && segment_end_rel <= plot_time(end)
            color_idx = mod(s-1, length(colors))+1;
            xline(segment_end_rel, '--', 'Color', colors{color_idx});
        end
    end
    
    % Calculate global mean and std values (these should be constants)
    global_mean = 1.0; % This is the normalized reference value
    
    % If we have global standard deviation in the config, use it
    if isfield(config, 'globalStd')
        global_std = config.globalStd;
    else
        % Otherwise estimate it from the data
        global_std = mean(std_norm_ispc);
    end
    
    % Add reference line at 1.0 (global mean)
    h_global_mean = yline(global_mean, '-', 'Global Mean', 'LineWidth', 1.5, 'Alpha', 0.6, 'Color', [0.5 0.5 0.5]);
    
    % Add standard deviation lines
    h_global_std_plus = yline(global_mean + global_std, ':', '+1 StdDev', 'LineWidth', 1, 'Alpha', 0.6, 'Color', [0.5 0.5 0.5]);
    h_global_std_minus = yline(global_mean - global_std, ':', '-1 StdDev', 'LineWidth', 1, 'Alpha', 0.6, 'Color', [0.5 0.5 0.5]);
    
    % Create proper display names for legend
    set(h_global_mean, 'DisplayName', 'Global Mean');
    set(h_global_std_plus, 'DisplayName', '+1 StdDev Global');
    set(h_global_std_minus, 'DisplayName', '-1 StdDev Global');
    
    title('Normalized ISPC Time Courses (as % of Global Mean)');
    xlabel('Time (s from protocol start)');
    ylabel('Normalized ISPC');
    legend('Location', 'eastoutside');
    ylim([y_min y_max]);  % Set calculated y-limits to ensure all data is visible
    xlim([0 plot_time(end) - plot_time(1)]);
    grid on;
    
    % Save the figure
    save_path = fullfile(config.results_dir, 'normalized_ispc_timecourses.png');
    saveas(gcf, save_path);
    fprintf('Saved normalized time courses to: %s\n', save_path);
    
    %% 2. Create topographic maps of normalized and percent change data
    if ~isempty(pre_stim_seg) && ~isempty(late_stim_seg)
        % Get segment ranges
        pre_range = segments.(pre_stim_seg);
        late_range = segments.(late_stim_seg);
        
        % Calculate mean values for segments
        norm_pre_mean = mean(norm_ispc_results(:, pre_range(1):pre_range(2)), 2, 'omitnan');
        norm_late_mean = mean(norm_ispc_results(:, late_range(1):late_range(2)), 2, 'omitnan');
        
        % Calculate z-scored means
        z_pre_mean = mean(z_scored_ispc(:, pre_range(1):pre_range(2)), 2, 'omitnan');
        z_late_mean = mean(z_scored_ispc(:, late_range(1):late_range(2)), 2, 'omitnan');
        
        % Calculate percent change
        pct_change = ((norm_late_mean - norm_pre_mean) ./ norm_pre_mean) * 100;
        
        % Calculate z-score difference
        z_diff = z_late_mean - z_pre_mean;
        
        % Create figure for topographic comparison
        figure('Name', 'Normalized ISPC Topography', 'Position', [100, 100, 1200, 400], 'visible', 'off');
        
        % Load excluded channels
        excluded_indices = [];
        if isfield(config, 'exclude_channels') && ~isempty(config.exclude_channels)
            for i = 1:length(config.exclude_channels)
                for c = 1:EEG.nbchan-1
                    if strcmpi(EEG.chanlocs(c).labels, config.exclude_channels{i})
                        excluded_indices = [excluded_indices, c];
                        break;
                    end
                end
            end
        end
        
        % Create a list of channels to include (all minus excluded)
        included_indices = setdiff(1:EEG.nbchan-1, excluded_indices);
        
        % Plot normalized pre-stim
        subplot(1, 3, 1);
        topoplot(norm_pre_mean(1:EEG.nbchan-1), EEG.chanlocs(1:end-1), 'maplimits', [0.5 1.5], 'electrodes', 'on', 'plotchans', included_indices);
        title('Normalized Pre-Stim ISPC');
        colorbar;
        
        % Plot normalized late-stim
        subplot(1, 3, 2);
        topoplot(norm_late_mean(1:EEG.nbchan-1), EEG.chanlocs(1:end-1), 'maplimits', [0.5 1.5], 'electrodes', 'on', 'plotchans', included_indices);
        title('Normalized Late-Stim ISPC');
        colorbar;
        
        % Plot percent change
        subplot(1, 3, 3);
        topoplot(pct_change(1:EEG.nbchan-1), EEG.chanlocs(1:end-1), 'maplimits', [-50 50], 'electrodes', 'on', 'plotchans', included_indices);
        title('Percent Change (Late vs Pre)');
        colorbar;
        
        sgtitle('Normalized ISPC Topography (% of Global Mean)');
        
        % Save the figure
        save_path = fullfile(config.results_dir, 'normalized_ispc_topography.png');
        saveas(gcf, save_path);
        fprintf('Saved normalized topography to: %s\n', save_path);
        
        % Create figure for z-scored comparison
        figure('Name', 'Z-Scored ISPC Topography', 'Position', [100, 100, 1200, 400], 'visible', 'off');
        
        % Plot z-scored pre-stim
        subplot(1, 3, 1);
        topoplot(z_pre_mean(1:EEG.nbchan-1), EEG.chanlocs(1:end-1), 'maplimits', [-2 2], 'electrodes', 'on', 'plotchans', included_indices);
        title('Z-Scored Pre-Stim ISPC');
        colorbar;
        
        % Plot z-scored late-stim
        subplot(1, 3, 2);
        topoplot(z_late_mean(1:EEG.nbchan-1), EEG.chanlocs(1:end-1), 'maplimits', [-2 2], 'electrodes', 'on', 'plotchans', included_indices);
        title('Z-Scored Late-Stim ISPC');
        colorbar;
        
        % Plot z-score difference
        subplot(1, 3, 3);
        topoplot(z_diff(1:EEG.nbchan-1), EEG.chanlocs(1:end-1), 'maplimits', [-2 2], 'electrodes', 'on', 'plotchans', included_indices);
        title('Z-Score Difference (Late - Pre)');
        colorbar;
        
        sgtitle('Z-Scored ISPC Topography');
        
        % Save the figure
        save_path = fullfile(config.results_dir, 'zscored_ispc_topography.png');
        saveas(gcf, save_path);
        fprintf('Saved z-scored topography to: %s\n', save_path);
    end
    
    fprintf('Normalized ISPC visualizations complete\n');
end
