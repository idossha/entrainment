function extract_with_reject()
% EXTRACT_WITH_REJECT - Extract data using the custom reject.m function
% This function extracts the first protocol by rejecting everything else

% Add EEGLAB and project paths
addpath('/Users/idohaber/Documents/MATLAB/eeglab2024.0/');
addpath('/Users/idohaber/Git-Projects/entrainment/src/main/data_processing/'); % For reject.m

% Start EEGLAB
eeglab nogui;

% Create test_data directory if it doesn't exist
test_data_dir = '/Users/idohaber/Git-Projects/entrainment/main/tests/data/';
if ~exist(test_data_dir, 'dir')
    mkdir(test_data_dir);
    fprintf('Created test data directory: %s\n', test_data_dir);
end

try
    % Load the dataset
    data_file = '/Users/idohaber/test_data/Strength_101_N1_forSW.set';
    fprintf('Loading dataset...\n');
    EEG = pop_loadset(data_file);
    
    % Get the protocol boundaries (known from inspection)
    % Protocol 1: Samples 285419-393454
    % Add 60 seconds before and after
    pre_samples = 60 * EEG.srate;
    protocol_start = max(1, 285419 - pre_samples);
    protocol_end = min(EEG.pnts, 393454 + pre_samples);
    
    fprintf('Protocol region: %d to %d (%.1f to %.1f seconds)\n', protocol_start, protocol_end, protocol_start/EEG.srate, protocol_end/EEG.srate);
    
    % Define regions to reject (everything except the protocol)
    reject_regions = [];
    
    % Reject from beginning to protocol start (if applicable)
    if protocol_start > 1
        reject_regions = [1, protocol_start-1];
    end
    
    % Reject from protocol end to end of data (if applicable)
    if protocol_end < EEG.pnts
        reject_regions = [reject_regions; protocol_end+1, EEG.pnts];
    end
    
    fprintf('Regions to reject:\n');
    disp(reject_regions);
    
    % Use the custom reject function
    if ~isempty(reject_regions)
        fprintf('Calling custom reject function...\n');
        EEG = reject(EEG, reject_regions);
        fprintf('Rejection completed.\n');
    else
        fprintf('No regions to reject.\n');
    end
    
    % Verify the new data length
    expected_length = protocol_end - protocol_start + 1;
    fprintf('Expected data length: %d samples\n', expected_length);
    fprintf('Actual data length: %d samples\n', EEG.pnts);
    
    % Rename the dataset
    EEG.setname = 'Protocol 1 Test Data';
    EEG.filename = 'protocol1_test_data.set';
    
    % Save the dataset
    save_file = fullfile(test_data_dir, 'protocol1_test_data.set');
    fprintf('Saving dataset to: %s\n', save_file);
    pop_saveset(EEG, 'filename', 'protocol1_test_data.set', 'filepath', test_data_dir);
    
    fprintf('Protocol 1 test dataset created successfully!\n');
    fprintf('Dataset duration: %.1f minutes\n', EEG.pnts/EEG.srate/60);
    
catch e
    fprintf('Error: %s\n', e.message);
    disp(e.stack);
end

end
