function [norm_final_topos, pct_change_topos, z_score_topos] = createNormalizedAggregateTopoplots(EEG, all_ispc_results, all_norm_ispc_results, all_stim_segments, config)
    % CREATENORMALIZEDAGGREGATETOPOPLOTS - Creates normalized topoplots across stimulations
    %
    % Inputs:
    %   EEG - EEG data structure
    %   all_ispc_results - Cell array with original ISPC results for each stimulation
    %   all_norm_ispc_results - Cell array with normalized ISPC results
    %   all_stim_segments - Cell array with segment definitions
    %   config - Configuration structure
    
    fprintf('Creating normalized aggregate topoplots across all stimulations...\n');
    
    % Number of stimulations
    num_stims = length(all_ispc_results);
    
    if num_stims < 1
        fprintf('No stimulation data found. Skipping aggregate topoplots.\n');
        norm_final_topos = [];
        pct_change_topos = [];
        z_score_topos = [];
        return;
    end
    
    % Define standard segment types
    standard_segments = {'pre_stim', 'early_stim', 'late_stim', 'post_stim'};
    display_names = {'Pre-Stim', 'Early-Stim', 'Late-Stim', 'Post-Stim'};
    
    % Initialize aggregated data for each segment type
    aggregated_norm_data = struct();
    aggregated_z_score_data = struct();
    
    for i = 1:length(standard_segments)
        segment_type = standard_segments{i};
        aggregated_norm_data.(segment_type) = [];
        aggregated_z_score_data.(segment_type) = [];
    end
    
    % Collect normalized ISPC values for each segment type across all stimulations
    for s = 1:num_stims
        norm_ispc_results = all_norm_ispc_results{s};
        stim_segments = all_stim_segments{s};
        
        % Calculate z-scored ISPC if not directly provided
        if exist('all_z_scored_ispc', 'var') && iscell(all_z_scored_ispc) && length(all_z_scored_ispc) >= s
            z_scored_ispc = all_z_scored_ispc{s};
        else
            % If all_z_scored_ispc is not provided or not a cell array,
            % we'll calculate z-scores based on the normalized data
            % This is a crude approximation, but better than crashing
            z_scored_ispc = (norm_ispc_results - 1) / 0.25;  % Assuming mean=1, std=0.25
            fprintf('Warning: Using approximated z-scores for stimulation %d\n', s);
        end
        
        for i = 1:length(standard_segments)
            segment_type = standard_segments{i};
            
            % Find matching segment in this stimulation
            found_segment = false;
            segment_names = fieldnames(stim_segments);
            
            for j = 1:length(segment_names)
                if contains(segment_names{j}, segment_type)
                    segment_range = stim_segments.(segment_names{j});
                    
                    % Calculate mean normalized ISPC over time for this segment
                    norm_segment_mean = mean(norm_ispc_results(:, segment_range(1):segment_range(2)), 2, 'omitnan');
                    
                    % Calculate mean z-scored ISPC over time for this segment
                    z_segment_mean = mean(z_scored_ispc(:, segment_range(1):segment_range(2)), 2, 'omitnan');
                    
                    % Add to aggregated data
                    if isempty(aggregated_norm_data.(segment_type))
                        aggregated_norm_data.(segment_type) = norm_segment_mean;
                        aggregated_z_score_data.(segment_type) = z_segment_mean;
                    else
                        aggregated_norm_data.(segment_type) = [aggregated_norm_data.(segment_type), norm_segment_mean];
                        aggregated_z_score_data.(segment_type) = [aggregated_z_score_data.(segment_type), z_segment_mean];
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
    norm_final_topos = struct();
    z_score_topos = struct();
    
    for i = 1:length(standard_segments)
        segment_type = standard_segments{i};
        if ~isempty(aggregated_norm_data.(segment_type))
            % Mean normalized ISPC across stimulations
            norm_final_topos.(segment_type) = mean(aggregated_norm_data.(segment_type), 2, 'omitnan');
            
            % Mean z-scored ISPC across stimulations
            z_score_topos.(segment_type) = mean(aggregated_z_score_data.(segment_type), 2, 'omitnan');
        end
    end
    
    % Calculate percent changes between segments
    pct_change_topos = struct();
    
    if isfield(norm_final_topos, 'pre_stim')
        if isfield(norm_final_topos, 'early_stim')
            pct_change_topos.early_vs_pre = ((norm_final_topos.early_stim - norm_final_topos.pre_stim) ./ norm_final_topos.pre_stim) * 100;
        end
        
        if isfield(norm_final_topos, 'late_stim')
            pct_change_topos.late_vs_pre = ((norm_final_topos.late_stim - norm_final_topos.pre_stim) ./ norm_final_topos.pre_stim) * 100;
        end
        
        if isfield(norm_final_topos, 'post_stim')
            pct_change_topos.post_vs_pre = ((norm_final_topos.post_stim - norm_final_topos.pre_stim) ./ norm_final_topos.pre_stim) * 100;
        end
    end
    
    if isfield(norm_final_topos, 'early_stim') && isfield(norm_final_topos, 'late_stim')
        pct_change_topos.late_vs_early = ((norm_final_topos.late_stim - norm_final_topos.early_stim) ./ norm_final_topos.early_stim) * 100;
    end
    
    % Calculate z-score differences between segments
    z_diff_topos = struct();
    
    if isfield(z_score_topos, 'pre_stim')
        if isfield(z_score_topos, 'early_stim')
            z_diff_topos.early_vs_pre = z_score_topos.early_stim - z_score_topos.pre_stim;
        end
        
        if isfield(z_score_topos, 'late_stim')
            z_diff_topos.late_vs_pre = z_score_topos.late_stim - z_score_topos.pre_stim;
        end
        
        if isfield(z_score_topos, 'post_stim')
            z_diff_topos.post_vs_pre = z_score_topos.post_stim - z_score_topos.pre_stim;
        end
    end
    
    if isfield(z_score_topos, 'early_stim') && isfield(z_score_topos, 'late_stim')
        z_diff_topos.late_vs_early = z_score_topos.late_stim - z_score_topos.early_stim;
    end
    
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
    fprintf('Found %d excluded channels, using %d channels for topoplots\n', ...
        length(excluded_indices), length(included_indices));
    
    % Create comprehensive figure with normalized segments and percent changes
    figure('Name', 'Normalized Aggregate Topoplots', 'Position', [50, 50, 1200, 800], 'visible', 'off');
    
    % Calculate total number of plots
    num_segments = sum(isfield(norm_final_topos, standard_segments));
    num_changes = length(fieldnames(pct_change_topos));
    total_plots = num_segments + num_changes;
    
    % Calculate a suitable grid layout
    grid_cols = ceil(sqrt(total_plots));
    grid_rows = ceil(total_plots / grid_cols);
    
    % Keep track of current subplot index
    plot_idx = 1;
    
    % First plot the normalized segments
    for i = 1:length(standard_segments)
        segment_type = standard_segments{i};
        if isfield(norm_final_topos, segment_type)
            subplot(grid_rows, grid_cols, plot_idx);
            
            % Create topoplot
            topoplot(norm_final_topos.(segment_type)(1:EEG.nbchan-1), EEG.chanlocs(1:end-1), ...
                'maplimits', [0.5 1.5], 'electrodes', 'on', 'plotchans', included_indices);
            
            title([display_names{i} ' (Normalized)']);
            colorbar;
            
            plot_idx = plot_idx + 1;
        end
    end
    
    % Then plot the percent changes
    change_fields = fieldnames(pct_change_topos);
    change_display_names = {
        'early_vs_pre', 'Early-Stim vs Pre-Stim';
        'late_vs_pre', 'Late-Stim vs Pre-Stim';
        'post_vs_pre', 'Post-Stim vs Pre-Stim';
        'late_vs_early', 'Late-Stim vs Early-Stim'
    };
    
    for i = 1:length(change_fields)
        change_type = change_fields{i};
        
        subplot(grid_rows, grid_cols, plot_idx);
        
        % Find display name for this change
        display_idx = find(strcmp(change_display_names(:,1), change_type));
        if ~isempty(display_idx)
            display_name = change_display_names{display_idx, 2};
        else
            display_name = strrep(change_type, '_', ' ');
        end
        
        % Create topoplot
        topoplot(pct_change_topos.(change_type)(1:EEG.nbchan-1), EEG.chanlocs(1:end-1), ...
            'maplimits', [-50 50], 'electrodes', 'on', 'plotchans', included_indices);
        
        title(['% Change: ' display_name]);
        colorbar;
        
        plot_idx = plot_idx + 1;
    end
    
    sgtitle('Normalized Aggregate Topoplots Across All Stimulations', 'FontSize', 16);
    
    % Save the comprehensive figure
    saveas(gcf, fullfile(config.results_dir, 'normalized_aggregate_topoplots.png'));
    fprintf('Saved normalized aggregate topoplots to: %s\n', fullfile(config.results_dir, 'normalized_aggregate_topoplots.png'));
    
    % Create a second figure for Z-score topoplots
    figure('Name', 'Z-Score Aggregate Topoplots', 'Position', [50, 50, 1200, 800], 'visible', 'off');
    
    % Calculate total number of plots for Z-score figure (segments and differences)
    num_z_segments = sum(isfield(z_score_topos, standard_segments));
    num_z_diffs = length(fieldnames(z_diff_topos));
    total_z_plots = num_z_segments + num_z_diffs;
    
    % Calculate a suitable grid layout
    z_grid_cols = ceil(sqrt(total_z_plots));
    z_grid_rows = ceil(total_z_plots / z_grid_cols);
    
    % Keep track of current subplot index
    z_plot_idx = 1;
    
    % First plot the Z-scored segments
    for i = 1:length(standard_segments)
        segment_type = standard_segments{i};
        if isfield(z_score_topos, segment_type)
            subplot(z_grid_rows, z_grid_cols, z_plot_idx);
            
            % Create topoplot
            topoplot(z_score_topos.(segment_type)(1:EEG.nbchan-1), EEG.chanlocs(1:end-1), ...
                'maplimits', [-2 2], 'electrodes', 'on', 'plotchans', included_indices);
            
            title([display_names{i} ' (Z-Score)']);
            colorbar;
            
            z_plot_idx = z_plot_idx + 1;
        end
    end
    
    % Then plot the Z-score differences
    z_diff_fields = fieldnames(z_diff_topos);
    
    for i = 1:length(z_diff_fields)
        z_diff_type = z_diff_fields{i};
        
        subplot(z_grid_rows, z_grid_cols, z_plot_idx);
        
        % Find display name for this difference
        display_idx = find(strcmp(change_display_names(:,1), z_diff_type));
        if ~isempty(display_idx)
            display_name = change_display_names{display_idx, 2};
        else
            display_name = strrep(z_diff_type, '_', ' ');
        end
        
        % Create topoplot
        topoplot(z_diff_topos.(z_diff_type)(1:EEG.nbchan-1), EEG.chanlocs(1:end-1), ...
            'maplimits', [-2 2], 'electrodes', 'on', 'plotchans', included_indices);
        
        title(['Z-Score Diff: ' display_name]);
        colorbar;
        
        z_plot_idx = z_plot_idx + 1;
    end
    
    sgtitle('Z-Score Aggregate Topoplots Across All Stimulations', 'FontSize', 16);
    
    % Save the Z-score figure
    saveas(gcf, fullfile(config.results_dir, 'zscored_aggregate_topoplots.png'));
    fprintf('Saved Z-score aggregate topoplots to: %s\n', fullfile(config.results_dir, 'zscored_aggregate_topoplots.png'));
    
    % Save the subject-level normalized average data
    subject_id = config.subject_id; % Assuming subject_id is in config
    save_path = fullfile(config.results_dir, sprintf('subject_%s_normalized_average.mat', subject_id));
    save(save_path, 'norm_final_topos', 'pct_change_topos', 'z_score_topos', 'z_diff_topos');
    fprintf('Saved subject-level normalized average data to: %s\n', save_path);
end
