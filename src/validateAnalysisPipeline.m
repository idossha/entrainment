function validateAnalysisPipeline()
    % VALIDATEANALYSISPIPELINE - Checks if all required paths and dependencies exist
    %
    % This function validates if all the paths, files, and dependencies required
    % for the EEG entrainment analysis pipeline exist.
    
    fprintf('Validating EEG Entrainment Analysis Pipeline...\n\n');
    
    % Define the base paths
    code_path = '/Users/idohaber/Git-Projects/entrainment';
    src_path = fullfile(code_path, 'src');
    assets_path = fullfile(code_path, 'src/assets');
    eeglab_path = '/Users/idohaber/Documents/MATLAB/eeglab2024.0/';
    
    % Validate experiment path (this might be different for each user)
    experiment_path = '/Volumes/Ido/analyze';
    if ~exist(experiment_path, 'dir')
        warning('Experiment path does not exist: %s\nYou will need to specify a correct path when running the analysis.', experiment_path);
    else
        fprintf('✓ Experiment path exists: %s\n', experiment_path);
    end
    
    % Check required source directories
    required_dirs = {
        src_path, 'Source code directory';
        assets_path, 'Assets directory';
        fullfile(src_path, 'eeglab'), 'EEGLAB custom functions directory';
        fullfile(src_path, 'utils'), 'Utilities directory'
    };
    
    fprintf('\nChecking required directories:\n');
    for i = 1:size(required_dirs, 1)
        dir_path = required_dirs{i, 1};
        description = required_dirs{i, 2};
        
        if exist(dir_path, 'dir')
            fprintf('✓ %s exists: %s\n', description, dir_path);
        else
            warning('%s does not exist: %s', description, dir_path);
        end
    end
    
    % Check required files
    required_files = {
        fullfile(assets_path, 'exclude.json'), 'Excluded channels JSON file';
        fullfile(assets_path, 'subject_condition.json'), 'Subject condition JSON file';
        fullfile(src_path, 'utils/egi256_GSN_HydroCel.sfp'), 'EGI electrode positions file';
        fullfile(assets_path, 'GSN-HydroCel-256.sfp'), 'Alternative EGI electrode positions file';
        fullfile(src_path, 'eeglab/pop_interp.m'), 'EEGLAB interpolation function';
        fullfile(src_path, 'eeglab/eeg_topoplot.m'), 'EEGLAB topoplot function';
        fullfile(src_path, 'eeglab/pop_topoplot.m'), 'EEGLAB pop_topoplot function'
    };
    
    fprintf('\nChecking required files:\n');
    for i = 1:size(required_files, 1)
        file_path = required_files{i, 1};
        description = required_files{i, 2};
        
        if exist(file_path, 'file')
            fprintf('✓ %s exists: %s\n', description, file_path);
        else
            warning('%s does not exist: %s', description, file_path);
        end
    end
    
    % Check for external EEGLAB
    if exist(eeglab_path, 'dir')
        fprintf('✓ External EEGLAB directory exists: %s\n', eeglab_path);
    else
        warning('External EEGLAB directory does not exist: %s', eeglab_path);
    end
    
    % Check for required MATLAB functions
    if isempty(which('readlocs'))
        warning('EEGLAB function "readlocs" not found in MATLAB path. Make sure EEGLAB is initialized.');
    else
        fprintf('✓ EEGLAB function "readlocs" is available.\n');
    end
    
    if isempty(which('topoplot'))
        warning('EEGLAB function "topoplot" not found in MATLAB path. Make sure EEGLAB is initialized.');
    else
        fprintf('✓ EEGLAB function "topoplot" is available.\n');
    end
    
    % Check for scripts in the pipeline
    pipeline_scripts = {
        fullfile(src_path, 'interpolateSubjectData.m'), 'Channel interpolation script';
        fullfile(src_path, 'runGroupAnalysisWithInterp.m'), 'Group analysis script';
        fullfile(src_path, 'runProtocolGroupAnalysis.m'), 'Protocol group analysis script';
        fullfile(src_path, 'setConfiguration.m'), 'Configuration script';
        fullfile(src_path, 'main.m'), 'Main analysis script'
    };
    
    fprintf('\nChecking pipeline scripts:\n');
    for i = 1:size(pipeline_scripts, 1)
        script_path = pipeline_scripts{i, 1};
        description = pipeline_scripts{i, 2};
        
        if exist(script_path, 'file')
            fprintf('✓ %s exists: %s\n', description, script_path);
        else
            warning('%s does not exist: %s', description, script_path);
        end
    end
    
    % Validate paths in scripts
    fprintf('\nValidating paths in scripts:\n');
    validateScriptPaths(fullfile(src_path, 'interpolateSubjectData.m'), code_path, eeglab_path);
    validateScriptPaths(fullfile(src_path, 'runGroupAnalysisWithInterp.m'), code_path, eeglab_path);
    validateScriptPaths(fullfile(src_path, 'runProtocolGroupAnalysis.m'), code_path, eeglab_path);
    
    fprintf('\nValidation complete. Address any warnings before running the analysis pipeline.\n');
    
    % Display the recommended pipeline execution order
    fprintf('\n=== Recommended Pipeline Execution Order ===\n');
    fprintf('1. Run interpolateSubjectData.m to create interpolated datasets.\n');
    fprintf('2. Run runGroupAnalysisWithInterp.m to generate global group averages.\n');
    fprintf('3. Run runProtocolGroupAnalysis.m to create protocol-specific averages.\n\n');
    
    % Provide example commands
    fprintf('Example MATLAB commands:\n');
    fprintf('>> interpolateSubjectData(''%s'', ''N1'');\n', experiment_path);
    fprintf('>> runGroupAnalysisWithInterp(''%s'', ''N1'');\n', experiment_path);
    fprintf('>> runProtocolGroupAnalysis(''%s'', ''N1'', 5);\n', experiment_path);
    fprintf('(The last parameter of runProtocolGroupAnalysis specifies the maximum number of stimulation protocols to process)\n');
end

function validateScriptPaths(script_path, code_path, eeglab_path)
    if ~exist(script_path, 'file')
        return;  % Skip if script doesn't exist
    end
    
    % Read the script content
    fid = fopen(script_path, 'r');
    if fid == -1
        warning('Could not open script for validation: %s', script_path);
        return;
    end
    
    script_content = fread(fid, '*char')';
    fclose(fid);
    
    [~, script_name, ~] = fileparts(script_path);
    
    % Check for hardcoded paths
    if contains(script_content, code_path)
        fprintf('✓ %s: Code path is correctly set to %s\n', script_name, code_path);
    else
        warning('%s: Code path might be incorrect. Expected: %s', script_name, code_path);
    end
    
    if contains(script_content, eeglab_path)
        fprintf('✓ %s: EEGLAB path is correctly set to %s\n', script_name, eeglab_path);
    else
        warning('%s: EEGLAB path might be incorrect. Expected: %s', script_name, eeglab_path);
    end
end