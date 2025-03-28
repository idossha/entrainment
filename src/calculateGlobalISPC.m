function [globalISPC, globalStd] = calculateGlobalISPC(EEG, stim_channel_idx, config)
    % CALCULATEGLOBALISPC Calculates global ISPC across the entire dataset
    %
    % Inputs:
    %   EEG - EEG data structure
    %   stim_channel_idx - Index of the stimulus channel
    %   config - Configuration structure
    %
    % Outputs:
    %   globalISPC - Global mean ISPC across all channels and time
    %   globalStd - Global standard deviation of ISPC
    
    fprintf('Calculating global ISPC across entire dataset...\n');
    
    % Use the entire dataset
    time_range = [1, EEG.pnts];
    
    % Calculate ISPC for the full dataset
    [ispc_results, ~, ~] = calculateISPC(EEG, stim_channel_idx, time_range, config);
    
    % Calculate global mean and standard deviation
    globalISPC = mean(mean(ispc_results, 'omitnan'), 'omitnan');
    globalStd = std(reshape(ispc_results, [], 1), 'omitnan');
    
    fprintf('Global ISPC mean: %.4f, std: %.4f\n', globalISPC, globalStd);
end
