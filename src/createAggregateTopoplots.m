function [final_topos, diff_topos] = createAggregateTopoplots(EEG, all_ispc_results, all_stim_segments, config, save_data)
    % Creates topoplots showing pre, early, late, post conditions and their differences
    % collapsed across all protocols
    % 
    % Inputs:
    %   EEG - EEG structure
    %   all_ispc_results - Cell array of ISPC results for all stimulations
    %   all_stim_segments - Cell array of segment definitions for all stimulations
    %   config - Configuration structure
    %   save_data - (Optional) Whether to save aggregate data (default: false)
    %
    % Outputs:
    %   final_topos - Structure containing averaged ISPC values for each segment
    %   diff_topos - Structure containing difference maps between segments
    
    fprintf('Creating aggregate topoplots across all stimulations...\n');
    
    % Check optional input
    if nargin < 5
        save_data = false;
    end
    
    % Number of stimulations
    num_stims = length(all_ispc_results);
    
    if num_stims < 1
        fprintf('No stimulation data found. Skipping aggregate topoplots.\n');
        final_topos = [];
        diff_topos = [];
        return;
    end
    
    % Define standard segment types we want to include
    standard_segments = {'pre_stim', 'early_stim', 'late_stim', 'post_stim'};
    display_names = {'Pre-Stim', 'Early-Stim', 'Late-Stim', 'Post-Stim'};
    
    % Find all segments across all stimulations
    all_segment_types = {};
    for s = 1:num_stims
        stim_segments = all_stim_segments{s};
        segment_names = fieldnames(stim_segments);
        
        for i = 1:length(segment_names)
            segment = segment_names{i};
            if ~ismember(segment, all_segment_types)
                all_segment_types{end+1} = segment;
            end
        end
    end
    
    % Filter for segments that match our standard names
    segment_indices = struct();
    for i = 1:length(standard_segments)
        segment_indices.(standard_segments{i}) = [];
        % Find all segments matching this standard type
        matches = find(contains(all_segment_types, standard_segments{i}));
        if ~isempty(matches)
            segment_indices.(standard_segments{i}) = matches;
        end
    end
    
    % Initialize aggregated data for each segment type
    aggregated_data = struct();
    for i = 1:length(standard_segments)
        segment_type = standard_segments{i};
        if isfield(segment_indices, segment_type) && ~isempty(segment_indices.(segment_type))
            aggregated_data.(segment_type) = [];
        end
    end
    
    % Collect and aggregate ISPC values for each segment type across all stimulations
    for s = 1:num_stims
        ispc_results = all_ispc_results{s};
        stim_segments = all_stim_segments{s};
        
        for i = 1:length(standard_segments)
            segment_type = standard_segments{i};
            
            % Find matching segment in this stimulation
            found_segment = false;
            segment_names = fieldnames(stim_segments);
            
            for j = 1:length(segment_names)
                if contains(segment_names{j}, segment_type)
                    segment_range = stim_segments.(segment_names{j});
                    
                    % Calculate mean ISPC over time for this segment
                    segment_mean = mean(ispc_results(:, segment_range(1):segment_range(2)), 2);
                    
                    % Add to aggregated data
                    if isempty(aggregated_data.(segment_type))
                        aggregated_data.(segment_type) = segment_mean;
                    else
                        aggregated_data.(segment_type) = [aggregated_data.(segment_type), segment_mean];
                    end
                    
                    found_segment = true;
                    break;
                end
            end
            
            if ~found_segment
                fprintf('Segment %s not found in stimulation %d\n', segment_type, s);
            end
        end
    end
    
    % Calculate mean across stimulations for each segment type
    final_topos = struct();
    for i = 1:length(standard_segments)
        segment_type = standard_segments{i};
        if isfield(aggregated_data, segment_type) && ~isempty(aggregated_data.(segment_type))
            % Mean across stimulations
            final_topos.(segment_type) = mean(aggregated_data.(segment_type), 2);
        end
    end
    
    % Check if we have all required segments
    if ~all(isfield(final_topos, standard_segments))
        missing = standard_segments(~isfield(final_topos, standard_segments));
        fprintf('Warning: Missing segments for aggregate topoplots: %s\n', strjoin(missing, ', '));
        % Continue with what we have
    end
    
    % Calculate differences between segments
    diff_topos = struct();
    
    % Define differences we want to calculate (same as in visualizeISPCTopoplots.m)
    if isfield(final_topos, 'pre_stim')
        if isfield(final_topos, 'early_stim')
            diff_topos.early_minus_pre = final_topos.early_stim - final_topos.pre_stim;
        end
        
        if isfield(final_topos, 'late_stim')
            diff_topos.late_minus_pre = final_topos.late_stim - final_topos.pre_stim;
        end
        
        if isfield(final_topos, 'post_stim')
            diff_topos.post_minus_pre = final_topos.post_stim - final_topos.pre_stim;
        end
    end
    
    if isfield(final_topos, 'early_stim') && isfield(final_topos, 'late_stim')
        diff_topos.late_minus_early = final_topos.late_stim - final_topos.early_stim;
    end
    
    % Find min/max values for consistent color scales
    segment_values = [];
    for i = 1:length(standard_segments)
        segment_type = standard_segments{i};
        if isfield(final_topos, segment_type)
            segment_values = [segment_values; final_topos.(segment_type)];
        end
    end
    
    segment_min = min(segment_values);
    segment_max = max(segment_values);
    
    % Find min/max values for difference maps
    diff_values = [];
    diff_fields = fieldnames(diff_topos);
    for i = 1:length(diff_fields)
        diff_values = [diff_values; diff_topos.(diff_fields{i})];
    end
    
    diff_abs_max = max(abs(diff_values));
    diff_limit = [-diff_abs_max diff_abs_max]; % Symmetric colormap for differences
    
    % Create comprehensive figure with segments and differences
    figure('Name', 'Aggregate Topoplots', 'Position', [50, 50, 1200, 800], 'visible', 'off');
    
    % Calculate total number of plots
    num_segments = sum(isfield(final_topos, standard_segments));
    num_diffs = length(fieldnames(diff_topos));
    total_plots = num_segments + num_diffs;
    
    % Calculate a suitable grid layout
    grid_cols = ceil(sqrt(total_plots));
    grid_rows = ceil(total_plots / grid_cols);
    
    % Keep track of current subplot index
    plot_idx = 1;
    
    % First plot the segments
    for i = 1:length(standard_segments)
        segment_type = standard_segments{i};
        if isfield(final_topos, segment_type)
            subplot(grid_rows, grid_cols, plot_idx);
            
            % Create topoplot
            topoplot(final_topos.(segment_type), EEG.chanlocs(1:end-1), ...
                'maplimits', [segment_min segment_max], 'electrodes', 'on');
            
            title(display_names{i});
            colorbar;
            
            plot_idx = plot_idx + 1;
        end
    end
    
    % Then plot the differences
    diff_fields = fieldnames(diff_topos);
    diff_display_names = {
        'early_minus_pre', 'Early-Stim - Pre-Stim';
        'late_minus_pre', 'Late-Stim - Pre-Stim';
        'post_minus_pre', 'Post-Stim - Pre-Stim';
        'late_minus_early', 'Late-Stim - Early-Stim'
    };
    
    for i = 1:length(diff_fields)
        diff_type = diff_fields{i};
        
        subplot(grid_rows, grid_cols, plot_idx);
        
        % Find display name for this difference
        display_idx = find(strcmp(diff_display_names(:,1), diff_type));
        if ~isempty(display_idx)
            display_name = diff_display_names{display_idx, 2};
        else
            display_name = strrep(diff_type, '_', ' ');
        end
        
        % Create topoplot
        topoplot(diff_topos.(diff_type), EEG.chanlocs(1:end-1), ...
            'maplimits', diff_limit, 'electrodes', 'on');
        
        title(display_name);
        colorbar;
        
        plot_idx = plot_idx + 1;
    end
    
    sgtitle('Aggregate Topoplots Across All Stimulations', 'FontSize', 16);
    
    % Save the comprehensive figure
    saveas(gcf, fullfile(config.results_dir, 'aggregate_topoplots.png'));
    fprintf('Saved aggregate topoplots to: %s\n', fullfile(config.results_dir, 'aggregate_topoplots.png'));
    
    % Save the subject-level average data if requested
    if save_data
        subject_id = config.subject_id; % Assuming subject_id is in config
        save_path = fullfile(config.results_dir, sprintf('subject_%s_average.mat', subject_id));
        save(save_path, 'final_topos', 'diff_topos');
        fprintf('Saved subject-level average data to: %s\n', save_path);
    end
end