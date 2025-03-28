function time_range = calculateTimeRange(stim_segments)
    segment_names = fieldnames(stim_segments);
    all_segment_bounds = zeros(length(segment_names), 2);
    
    for seg_idx = 1:length(segment_names)
        all_segment_bounds(seg_idx, :) = stim_segments.(segment_names{seg_idx});
    end
    
    time_range = [min(all_segment_bounds(:,1)), max(all_segment_bounds(:,2))];
end
