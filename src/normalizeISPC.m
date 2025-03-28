function [norm_ispc_results, z_scored_ispc] = normalizeISPC(ispc_results, stim_segments, globalISPC, globalStd)
    % NORMALIZEISPC Normalizes ISPC values using the global ISPC
    %
    % Inputs:
    %   ispc_results - Matrix of ISPC values [channels x time]
    %   stim_segments - Structure with segment definitions
    %   globalISPC - Global mean ISPC (precalculated across entire dataset)
    %   globalStd - Global standard deviation of ISPC
    %
    % Outputs:
    %   norm_ispc_results - Normalized ISPC values (as % of global mean)
    %   z_scored_ispc - Z-scored ISPC values for each channel
    
    fprintf('Normalizing ISPC values using global ISPC: %.4f\n', globalISPC);
    
    % Normalize ISPC values as percentage of global mean
    norm_ispc_results = ispc_results / globalISPC;
    
    % Calculate z-score based on global mean and std
    z_scored_ispc = (ispc_results - globalISPC) / globalStd;
    
    % Calculate normalized segment means and percent changes
    segment_names = fieldnames(stim_segments);
    
    % Find pre-stim and late-stim segments (for percent change calculation)
    pre_stim_seg = '';
    late_stim_seg = '';
    
    % Look for pre_stim and late_stim segments
    for i = 1:length(segment_names)
        if strcmp(segment_names{i}, 'pre_stim')
            pre_stim_seg = segment_names{i};
        elseif strcmp(segment_names{i}, 'late_stim')
            late_stim_seg = segment_names{i};
        end
    end
    
    % If standard segments found, calculate percent change
    if ~isempty(pre_stim_seg) && ~isempty(late_stim_seg)
        pre_range = stim_segments.(pre_stim_seg);
        late_range = stim_segments.(late_stim_seg);
        
        % Calculate normalized segment means
        norm_pre_mean = mean(norm_ispc_results(:, pre_range(1):pre_range(2)), 2, 'omitnan');
        norm_late_mean = mean(norm_ispc_results(:, late_range(1):late_range(2)), 2, 'omitnan');
        
        % Calculate percent change
        pct_change = ((norm_late_mean - norm_pre_mean) ./ norm_pre_mean) * 100;
        
        fprintf('Mean normalized pre-stim ISPC: %.4f\n', mean(norm_pre_mean, 'omitnan'));
        fprintf('Mean normalized late-stim ISPC: %.4f\n', mean(norm_late_mean, 'omitnan'));
        fprintf('Mean percent change: %.2f%%\n', mean(pct_change, 'omitnan'));
    end
    
    fprintf('ISPC normalization complete\n');
end
