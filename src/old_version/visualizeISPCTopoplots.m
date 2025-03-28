function visualizeISPCTopoplots(EEG, ispc_results, segments, config)
    fprintf('Creating topoplots for each segment and difference maps...\n');
    
    % Get segment names
    segment_names = fieldnames(segments);
    
    % Define the standard segments we want to find
    standard_segments = {'pre_stim', 'early_stim', 'late_stim', 'post_stim'};
    segment_display_names = {'Pre-Stim', 'Early-Stim', 'Late-Stim', 'Post-Stim'};
    
    % Structure to hold the found segments
    found_segments = struct();
    
    % Loop through standard segments and find them
    for i = 1:length(standard_segments)
        seg_name = standard_segments{i};
        
        % Look for segments matching this pattern, preferring the first stimulation
        for j = 1:length(segment_names)
            current_seg = segment_names{j};
            if contains(current_seg, seg_name) && (endsWith(current_seg, '_1') || ~contains(current_seg, '_'))
                found_segments.(seg_name) = current_seg;
                fprintf('Found %s segment: %s\n', seg_name, current_seg);
                break;
            end
        end
        
        % If not found in stim_1 or without suffix, look for any matching segment
        if ~isfield(found_segments, seg_name)
            for j = 1:length(segment_names)
                current_seg = segment_names{j};
                if contains(current_seg, seg_name)
                    found_segments.(seg_name) = current_seg;
                    fprintf('Found %s segment (alternative): %s\n', seg_name, current_seg);
                    break;
                end
            end
        end
    end
    
    % Check if we found the segments we need
    available_segments = fieldnames(found_segments);
    fprintf('Found %d out of %d standard segments: %s\n', length(available_segments), length(standard_segments), ...
        strjoin(available_segments, ', '));
    
    % Need at least 2 segments to proceed
    if length(available_segments) < 2
        warning('Not enough segments found for topoplots. Need at least 2 segments. Skipping.');
        return;
    end
    
    % Filter out excluded channels
    included_channels = true(1, EEG.nbchan-1); % Last channel is the stimulus
    
    % Check for excluded channels
    if isfield(config, 'exclude_channels') && ~isempty(config.exclude_channels)
        fprintf('Checking %d excluded channels...\n', length(config.exclude_channels));
        
        excluded_count = 0;
        for c = 1:EEG.nbchan-1
            if isempty(EEG.chanlocs(c).labels)
                fprintf('Warning: Channel %d has no label\n', c);
                continue;
            end
            
            channel_label = EEG.chanlocs(c).labels;
            if ismember(channel_label, config.exclude_channels)
                included_channels(c) = false;
                excluded_count = excluded_count + 1;
            end
        end
        
        fprintf('Marked %d out of %d channels for exclusion\n', excluded_count, EEG.nbchan-1);
    else
        fprintf('No channels marked for exclusion\n');
    end
    
    % Find excluded channels by index
    excluded_indices = find(~included_channels);
    
    % Get channel locations for topoplot
    included_chanlocs = EEG.chanlocs(1:end-1);
    
    % Calculate average ISPC for each segment
    segment_data = struct();
    
    for i = 1:length(available_segments)
        seg_name = available_segments{i};
        segment = found_segments.(seg_name);
        segment_range = segments.(segment);
        
        % Calculate average ISPC for this segment
        avg_ispc = mean(ispc_results(:, segment_range(1):segment_range(2)), 2);
        
        % Store in both ways for easier access
        segment_data.(seg_name) = avg_ispc;
    end
    
    % Calculate min/max ISPC across segments for consistent color scaling
    all_avg_ispc = [];
    for i = 1:length(available_segments)
        seg_name = available_segments{i};
        all_avg_ispc = [all_avg_ispc; segment_data.(seg_name)(included_channels)];
    end
    
    % Remove any NaN or Inf values for min/max calculation
    all_avg_ispc = all_avg_ispc(~isnan(all_avg_ispc) & ~isinf(all_avg_ispc));
    
    ispc_min = min(all_avg_ispc);
    ispc_max = max(all_avg_ispc);
    
    fprintf('ISPC range for colormap: [%.4f, %.4f]\n', ispc_min, ispc_max);
    
    % Create directory if it doesn't exist (safeguard)
    if ~exist(config.results_dir, 'dir')
        [success, msg] = mkdir(config.results_dir);
        if ~success
            warning('Failed to create directory %s: %s', config.results_dir, msg);
            fprintf('Will attempt to save to parent directory instead.\n');
            [parent_dir, ~, ~] = fileparts(config.results_dir);
            config.results_dir = parent_dir;
        else
            fprintf('Created directory: %s\n', config.results_dir);
        end
    end
    
    % Create a table of segments so we can order them correctly
    segment_table = struct('name', {}, 'display_name', {}, 'data', {});
    
    % Map standard segments to their display names
    for i = 1:length(standard_segments)
        std_seg = standard_segments{i};
        if isfield(found_segments, std_seg)
            idx = length(segment_table) + 1;
            segment_table(idx).name = std_seg;
            segment_table(idx).display_name = segment_display_names{i};
            segment_table(idx).data = segment_data.(std_seg);
        end
    end
    
    % Define differences we want to calculate
    diff_pairs = {};
    
    % Check which segments we have available to make differences
    if isfield(found_segments, 'pre_stim')
        % Early-stim vs pre-stim
        if isfield(found_segments, 'early_stim')
            diff_pairs{end+1} = struct('name', 'early_minus_pre', ...
                'display_name', 'Early-Stim - Pre-Stim', ...
                'minuend', 'early_stim', 'subtrahend', 'pre_stim');
        end
        
        % Late-stim vs pre-stim
        if isfield(found_segments, 'late_stim')
            diff_pairs{end+1} = struct('name', 'late_minus_pre', ...
                'display_name', 'Late-Stim - Pre-Stim', ...
                'minuend', 'late_stim', 'subtrahend', 'pre_stim');
        end
        
        % Post-stim vs pre-stim
        if isfield(found_segments, 'post_stim')
            diff_pairs{end+1} = struct('name', 'post_minus_pre', ...
                'display_name', 'Post-Stim - Pre-Stim', ...
                'minuend', 'post_stim', 'subtrahend', 'pre_stim');
        end
    end
    
    % Late-stim vs early-stim
    if isfield(found_segments, 'late_stim') && isfield(found_segments, 'early_stim')
        diff_pairs{end+1} = struct('name', 'late_minus_early', ...
            'display_name', 'Late-Stim - Early-Stim', ...
            'minuend', 'late_stim', 'subtrahend', 'early_stim');
    end
    
    % Calculate difference data
    diff_data = struct();
    
    for i = 1:length(diff_pairs)
        diff_pair = diff_pairs{i};
        diff_data.(diff_pair.name) = segment_data.(diff_pair.minuend) - segment_data.(diff_pair.subtrahend);
    end
    
    % Find min/max values across all difference data for symmetric colormap
    all_diffs = [];
    diff_names = fieldnames(diff_data);
    
    for i = 1:length(diff_names)
        valid_diff = diff_data.(diff_names{i})(included_channels);
        valid_diff = valid_diff(~isnan(valid_diff) & ~isinf(valid_diff));
        all_diffs = [all_diffs; valid_diff];
    end
    
    diff_abs_max = max(abs(all_diffs));
    diff_limit = [-diff_abs_max diff_abs_max]; % Symmetric colormap
    
    fprintf('Difference ISPC range for colormap: [%.4f, %.4f]\n', -diff_abs_max, diff_abs_max);
    
    % Get number of plots and calculate layout
    num_segments = length(segment_table);
    num_diffs = length(diff_pairs);
    total_plots = num_segments + num_diffs;
    
    % Create comprehensive figure with segments and differences
    figure('Name', 'ISPC Segment and Difference Topoplots', 'Position', [50, 50, 1600, 900], 'visible', 'off');
    
    % Calculate a suitable grid layout
    grid_cols = ceil(sqrt(total_plots));
    grid_rows = ceil(total_plots / grid_cols);
    
    % First plot the segments
    for i = 1:num_segments
        subplot(grid_rows, grid_cols, i);
        
        % Get segment data
        topodata = segment_table(i).data;
        
        % Create a temporary data vector where we set excluded channels to NaN
        temp_data = topodata;
        if ~isempty(excluded_indices)
            temp_data(excluded_indices) = NaN;
        end
        
        % Plot with NaN values for excluded channels
        try
            topoplot(temp_data, included_chanlocs, 'maplimits', [ispc_min ispc_max], ...
                'electrodes', 'on');
            title(segment_table(i).display_name);
            colorbar;
        catch err
            warning('Error in topoplot for %s: %s', segment_table(i).display_name, err.message);
            % Try simpler approach
            try
                topoplot(temp_data, included_chanlocs);
                title(segment_table(i).display_name);
                colorbar;
            catch
                text(0.5, 0.5, 'Topoplot Error', 'HorizontalAlignment', 'center');
                axis off;
            end
        end
    end
    
    % Then plot the differences
    for i = 1:num_diffs
        subplot(grid_rows, grid_cols, num_segments + i);
        
        diff_pair = diff_pairs{i};
        diff_name = diff_pair.name;
        
        % Create a temporary data vector where we set excluded channels to NaN
        temp_diff_data = diff_data.(diff_name);
        if ~isempty(excluded_indices)
            temp_diff_data(excluded_indices) = NaN;
        end
        
        % Plot with NaN values for excluded channels
        try
            topoplot(temp_diff_data, included_chanlocs, 'maplimits', diff_limit, ...
                'electrodes', 'on');
            title(diff_pair.display_name);
            colorbar;
        catch err
            warning('Error in difference topoplot for %s: %s', diff_pair.display_name, err.message);
            % Try simpler approach
            try
                topoplot(temp_diff_data, included_chanlocs);
                title(diff_pair.display_name);
                colorbar;
            catch
                text(0.5, 0.5, 'Topoplot Error', 'HorizontalAlignment', 'center');
                axis off;
            end
        end
    end
    
    sgtitle('ISPC by Segment and Segment Differences', 'FontSize', 16);
    
    % Save the comprehensive figure
    try
        save_path = fullfile(config.results_dir, 'ispc_comprehensive_topoplots.png');
        saveas(gcf, save_path);
        fprintf('Saved comprehensive topoplots to: %s\n', save_path);
    catch err
        warning('Failed to save comprehensive topoplots: %s', err.message);
    end
end
