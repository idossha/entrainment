function createStimulationComparison(EEG, all_ispc_results, all_stim_segments, config)
    % Creates comparative visualizations across all stimulation periods
    fprintf('Creating comparison visualizations across all stimulations...\n');
    
    % Number of stimulations
    num_stims = length(all_ispc_results);
    
    if num_stims < 2
        fprintf('Only one stimulation found. Skipping comparison.\n');
        return;
    end
    
    % Find common segment types across stimulations (pre_stim, early_stim, etc.)
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
    
    % Initialize data structures for comparison
    segment_ispc = cell(length(all_segment_types), num_stims);
    segment_topos = cell(length(all_segment_types), num_stims);
    
    % Collect ISPC values for each segment type and stimulation
    for s = 1:length(all_segment_types)
        segment_type = all_segment_types{s};
        
        for stim_idx = 1:num_stims
            ispc_results = all_ispc_results{stim_idx};
            stim_segments = all_stim_segments{stim_idx};
            
            if isfield(stim_segments, segment_type)
                segment_range = stim_segments.(segment_type);
                
                % Calculate mean ISPC over time for this segment
                segment_mean_over_time = mean(ispc_results(:, segment_range(1):segment_range(2)), 2);
                
                % Calculate overall mean (across channels and time)
                segment_ispc{s, stim_idx} = mean(segment_mean_over_time);
                
                % Store topographic data
                segment_topos{s, stim_idx} = segment_mean_over_time;
            end
        end
    end
    
    % Create bar plot comparison
    figure('Name', 'Stimulation Comparison', 'Position', [100, 100, 1200, 800], 'visible', 'off');
    
    % Plot each segment type
    for s = 1:length(all_segment_types)
        subplot(2, ceil(length(all_segment_types)/2), s);
        
        % Extract data for this segment type
        segment_data = cell2mat(segment_ispc(s, :));
        
        % Create bar plot
        bar(segment_data);
        
        % Add labels
        title(strrep(all_segment_types{s}, '_', ' '));
        xlabel('Stimulation Instance');
        ylabel('Mean ISPC');
        ylim([0 1]);
        grid on;
    end
    
    sgtitle('Mean ISPC by Segment Across Stimulation Instances');
    
    % Save figure
    saveas(gcf, fullfile(config.results_dir, 'stimulation_comparison_bar.png'));
    
    % Create topographic comparison for each segment type
    for s = 1:length(all_segment_types)
        segment_type = all_segment_types{s};
        
        % Skip if any stimulation is missing this segment
        topo_data = segment_topos(s, :);
        if any(cellfun(@isempty, topo_data))
            continue;
        end
        
        figure('Name', ['Topographic Comparison - ' segment_type], 'Position', [100, 100, 250*num_stims, 250], 'visible', 'off');
        
        % Find common color scale for consistent visualization
        all_values = [];
        for stim_idx = 1:num_stims
            all_values = [all_values; segment_topos{s, stim_idx}];
        end
        clim_min = min(all_values);
        clim_max = max(all_values);
        
        % Create topoplots for each stimulation
        for stim_idx = 1:num_stims
            subplot(1, num_stims, stim_idx);
            
            % Create topoplot
            topoplot(segment_topos{s, stim_idx}, EEG.chanlocs(1:end-1), ...
                'maplimits', [clim_min clim_max], 'electrodes', 'on');
            
            title(['Stim ' num2str(stim_idx)]);
            
            % Add colorbar to last plot
            if stim_idx == num_stims
                colorbar;
            end
        end
        
        sgtitle(['Topographic Comparison - ' strrep(segment_type, '_', ' ')]);
        
        % Save figure
        saveas(gcf, fullfile(config.results_dir, ['topo_compare_' segment_type '.png']));
    end
    
    % NEW CODE FOR DIFFERENCE TOPOPLOTS
    % Define the segments we want to compare (similar to what's in visualizeISPCTopoplots.m)
    if ismember('pre_stim', all_segment_types) && ismember('late_stim', all_segment_types)
        createDifferenceTopoplot(EEG, segment_topos, all_segment_types, 'late_stim', 'pre_stim', 'Late-Stim minus Pre-Stim', config);
    end
    
    if ismember('pre_stim', all_segment_types) && ismember('early_stim', all_segment_types)
        createDifferenceTopoplot(EEG, segment_topos, all_segment_types, 'early_stim', 'pre_stim', 'Early-Stim minus Pre-Stim', config);
    end
    
    if ismember('pre_stim', all_segment_types) && ismember('post_stim', all_segment_types)
        createDifferenceTopoplot(EEG, segment_topos, all_segment_types, 'post_stim', 'pre_stim', 'Post-Stim minus Pre-Stim', config);
    end
    
    if ismember('early_stim', all_segment_types) && ismember('late_stim', all_segment_types)
        createDifferenceTopoplot(EEG, segment_topos, all_segment_types, 'late_stim', 'early_stim', 'Late-Stim minus Early-Stim', config);
    end
    
    % Create detailed comparison table
    fprintf('\n==== Stimulation Comparison ====\n');
    fprintf('%-15s', 'Segment');
    for stim_idx = 1:num_stims
        fprintf('Stim %-9d', stim_idx);
    end
    fprintf('\n');
    fprintf('%s\n', repmat('-', 1, 15 + 15*num_stims));
    
    for s = 1:length(all_segment_types)
        fprintf('%-15s', strrep(all_segment_types{s}, '_', ' '));
        
        for stim_idx = 1:num_stims
            if ~isempty(segment_ispc{s, stim_idx})
                fprintf('%-15.3f', segment_ispc{s, stim_idx});
            else
                fprintf('%-15s', 'N/A');
            end
        end
        fprintf('\n');
    end
    
    % Calculate average across all stimulations for each segment type
    fprintf('\n==== Average Across All Stimulations ====\n');
    for s = 1:length(all_segment_types)
        segment_data = segment_ispc(s, :);
        valid_data = segment_data(~cellfun(@isempty, segment_data));
        
        if ~isempty(valid_data)
            fprintf('%-15s: %.3f\n', strrep(all_segment_types{s}, '_', ' '), mean(cell2mat(valid_data)));
        end
    end
    
    fprintf('\nStimulation comparison complete.\n');
end

% New helper function for creating difference topoplots
function createDifferenceTopoplot(EEG, segment_topos, all_segment_types, minuend_name, subtrahend_name, display_title, config)
    % Find indices for both segments
    minuend_idx = find(strcmp(all_segment_types, minuend_name));
    subtrahend_idx = find(strcmp(all_segment_types, subtrahend_name));
    
    if isempty(minuend_idx) || isempty(subtrahend_idx)
        fprintf('Could not find segments for comparison: %s and %s\n', minuend_name, subtrahend_name);
        return;
    end
    
    % Number of stimulations
    num_stims = size(segment_topos, 2);
    
    figure('Name', ['Difference Topoplot - ' display_title], 'Position', [100, 100, 250*num_stims, 250], 'visible', 'off');
    
    % Calculate all differences to find proper color scale
    all_diffs = [];
    for stim_idx = 1:num_stims
        % Skip any stimulations where data is missing
        if isempty(segment_topos{minuend_idx, stim_idx}) || isempty(segment_topos{subtrahend_idx, stim_idx})
            continue;
        end
        
        diff_data = segment_topos{minuend_idx, stim_idx} - segment_topos{subtrahend_idx, stim_idx};
        all_diffs = [all_diffs; diff_data];
    end
    
    % Compute symmetric color limits for difference maps
    diff_max = max(abs(all_diffs));
    diff_limit = [-diff_max diff_max];
    
    % Create a topoplot for each stimulation
    for stim_idx = 1:num_stims
        % Skip any stimulations where data is missing
        if isempty(segment_topos{minuend_idx, stim_idx}) || isempty(segment_topos{subtrahend_idx, stim_idx})
            continue;
        end
        
        % Calculate difference
        diff_data = segment_topos{minuend_idx, stim_idx} - segment_topos{subtrahend_idx, stim_idx};
        
        % Create subplot
        subplot(1, num_stims, stim_idx);
        
        % Create topoplot with symmetric color scale
        topoplot(diff_data, EEG.chanlocs(1:end-1), 'maplimits', diff_limit, 'electrodes', 'on');
        
        title(['Stim ' num2str(stim_idx)]);
        
        % Add colorbar to last plot
        if stim_idx == num_stims
            colorbar;
        end
    end
    
    sgtitle(['Difference Map: ' display_title]);
    
    % Create sanitized filename by replacing spaces with underscores
    filename = strrep(display_title, ' ', '_');
    filename = strrep(filename, '-', '_');
    
    % Save figure
    saveas(gcf, fullfile(config.results_dir, ['topo_compare_' lower(filename) '.png']));
    fprintf('Created difference topoplot: %s\n', display_title);
end
