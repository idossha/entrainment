function stim_segments = extractSegmentsForStimulation(segments, all_segment_names, stim_idx)
    stim_segments = struct();
    for i = 1:length(all_segment_names)
        segment_name = all_segment_names{i};
        if endsWith(segment_name, ['_' num2str(stim_idx)]) || (~contains(segment_name, '_') && stim_idx == 1)
            % Extract base segment name
            if contains(segment_name, '_')
                parts = strsplit(segment_name, '_');
                base_name = strjoin(parts(1:end-1), '_');
            else
                base_name = segment_name;
            end
            
            % Store with base name
            stim_segments.(base_name) = segments.(segment_name);
        end
    end
end
