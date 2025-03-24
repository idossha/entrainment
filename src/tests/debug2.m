%% Script to debug the structure of a results.mat file
% This will print the fields and their sizes in the file

% Pick a subject that exists (from the previous output)
subject_id = '101';
session_id = 'N1';

% Define the path to the results file
experiment_path = '/Volumes/Ido/analyze';
subject_path = fullfile(experiment_path, subject_id, session_id);
output_path = fullfile(subject_path, 'output');
entrainment_path = fullfile(output_path, sprintf('entrainment_%s_%s', subject_id, session_id));
stim_dir = fullfile(entrainment_path, 'stim_1');
results_file = fullfile(stim_dir, 'results.mat');

fprintf('Looking for file: %s\n', results_file);

% Check if file exists
if ~exist(results_file, 'file')
    error('Results file not found: %s', results_file);
end

% Load the file and display its structure
fprintf('Loading results file...\n');
results = load(results_file);

% Display the structure of the loaded data
fprintf('\n=== Structure of results.mat file ===\n');
fields = fieldnames(results);
for i = 1:length(fields)
    field = fields{i};
    value = results.(field);
    
    % Display information about the field
    fprintf('Field: %s\n', field);
    
    % Show type and size
    if isstruct(value)
        fprintf('  Type: struct\n');
        if isempty(fieldnames(value))
            fprintf('  (empty struct)\n');
        else
            sub_fields = fieldnames(value);
            for j = 1:length(sub_fields)
                sub_field = sub_fields{j};
                sub_value = value.(sub_field);
                
                % Display information about the subfield
                fprintf('  Subfield: %s\n', sub_field);
                
                % Show type and size of subfield
                if isnumeric(sub_value)
                    fprintf('    Type: numeric, Size: %s\n', mat2str(size(sub_value)));
                elseif ischar(sub_value)
                    fprintf('    Type: char, Size: %s\n', mat2str(size(sub_value)));
                elseif iscell(sub_value)
                    fprintf('    Type: cell, Size: %s\n', mat2str(size(sub_value)));
                elseif isstruct(sub_value)
                    fprintf('    Type: struct, Size: %s\n', mat2str(size(sub_value)));
                else
                    fprintf('    Type: other, Size: %s\n', mat2str(size(sub_value)));
                end
            end
        end
    elseif isnumeric(value)
        fprintf('  Type: numeric, Size: %s\n', mat2str(size(value)));
    elseif ischar(value)
        fprintf('  Type: char, Size: %s\n', mat2str(size(value)));
    elseif iscell(value)
        fprintf('  Type: cell, Size: %s\n', mat2str(size(value)));
    else
        fprintf('  Type: other, Size: %s\n', mat2str(size(value)));
    end
    
    fprintf('\n');
end

fprintf('=== End of structure ===\n');
