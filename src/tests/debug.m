%% Simple script to check results.mat structure
% Run this script to see what's in a results file

% Choose a subject
subject_id = '101';
session_id = 'N1';
experiment_path = '/Volumes/Ido/analyze';

% Construct path to results file
results_path = fullfile(experiment_path, subject_id, session_id, 'output', ...
    ['entrainment_' subject_id '_' session_id], 'stim_1', 'results.mat');

% Display path
disp(['Checking file: ' results_path]);

% Check if file exists
if ~exist(results_path, 'file')
    disp('File does not exist!');
    return;
end

% Load the file
disp('Loading file...');
data = load(results_path);

% Display fields
disp('Fields in results.mat:');
disp(fieldnames(data));

% Check for specific fields
if isfield(data, 'ispc_results')
    disp(['ispc_results size: ' mat2str(size(data.ispc_results))]);
else
    disp('ispc_results not found!');
end

if isfield(data, 'stim_segments')
    disp('stim_segments found, contains:');
    disp(fieldnames(data.stim_segments));
elseif isfield(data, 'segments')
    disp('segments found, contains:');
    disp(fieldnames(data.segments));
else
    disp('No segments or stim_segments found!');
    
    % Try to find any fields that might contain segment information
    fields = fieldnames(data);
    for i = 1:length(fields)
        if isstruct(data.(fields{i}))
            subfields = fieldnames(data.(fields{i}));
            if any(contains(lower(subfields), 'stim'))
                disp(['Found potential segment info in: ' fields{i}]);
                disp(subfields);
            end
        end
    end
end
