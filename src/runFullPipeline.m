function runFullPipeline(experiment_path, session_id, max_stim)
    % RUNFULLPIPELINE - Runs the complete EEG entrainment analysis pipeline
    %
    % This function runs the full analysis pipeline from subject-level interpolation
    % to group-level analysis for both global and protocol-specific averages.
    %
    % Inputs:
    %   experiment_path - Base path for the experiment data
    %   session_id - Session ID (e.g., 'N1')
    %   max_stim - Maximum number of stimulation protocols to process (default: 5)
    
    % Set default session if not provided
    if nargin < 2
        session_id = 'N1';
    end
    
    % Set default max stimulation number
    if nargin < 3
        max_stim = 5;
    end
    
    % Configure execution
    diary_file = fullfile(experiment_path, sprintf('pipeline_log_%s.txt', datestr(now, 'yyyymmdd_HHMMSS')));
    diary(diary_file);
    
    % Start tracking execution time
    start_time = tic;
    fprintf('Starting EEG Entrainment Analysis Pipeline at %s\n', datestr(now));
    fprintf('Experiment Path: %s\n', experiment_path);
    fprintf('Session ID: %s\n', session_id);
    fprintf('Maximum Stimulation Protocols: %d\n\n', max_stim);
    
    try
        % Step 1: Validate the pipeline
        fprintf('\n===== STEP 1: VALIDATING PIPELINE =====\n\n');
        validateAnalysisPipeline();
        
        % Step 2: Run subject-level interpolation
        fprintf('\n===== STEP 2: INTERPOLATING SUBJECT DATA =====\n\n');
        interpolation_start = tic;
        interpolateSubjectData(experiment_path, session_id);
        fprintf('\nInterpolation completed in %.2f minutes\n', toc(interpolation_start)/60);
        
        % Step 3: Run group analysis with interpolated data
        fprintf('\n===== STEP 3: RUNNING GROUP ANALYSIS =====\n\n');
        group_start = tic;
        runGroupAnalysisWithInterp(experiment_path, session_id);
        fprintf('\nGroup analysis completed in %.2f minutes\n', toc(group_start)/60);
        
        % Step 4: Run protocol-specific analysis
        fprintf('\n===== STEP 4: RUNNING PROTOCOL-SPECIFIC ANALYSIS =====\n\n');
        protocol_start = tic;
        runProtocolGroupAnalysis(experiment_path, session_id, max_stim);
        fprintf('\nProtocol analysis completed in %.2f minutes\n', toc(protocol_start)/60);
        
        % Final status
        fprintf('\n===== PIPELINE COMPLETED SUCCESSFULLY =====\n\n');
        fprintf('Total execution time: %.2f minutes\n', toc(start_time)/60);
        fprintf('Log saved to: %s\n', diary_file);
        
    catch err
        % Handle errors
        fprintf('\n===== ERROR IN PIPELINE EXECUTION =====\n\n');
        fprintf('Error message: %s\n', err.message);
        fprintf('Error in: %s (line %d)\n', err.stack(1).name, err.stack(1).line);
        fprintf('Stack trace:\n');
        disp(err.stack);
        fprintf('\nPipeline execution failed after %.2f minutes\n', toc(start_time)/60);
    end
    
    % End logging
    diary off;
    
    % Display final message
    fprintf('\nAnalysis complete. Output files are stored in the following locations:\n');
    fprintf('1. Subject-level interpolated data: <subject-dir>/output/entrainment_<subject>_<session>/subject_<subject>_interp.mat\n');
    fprintf('2. Group analysis results: %s/group_results/%s/\n', experiment_path, session_id);
    fprintf('3. Protocol-specific results: %s/group_results/%s/stim_*/\n', experiment_path, session_id);
    fprintf('4. Execution log: %s\n', diary_file);
end
