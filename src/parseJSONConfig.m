function parsedConfig = parseJSONConfig(config)
    % PARSEJSONCONFIG Parses JSON configuration files for EEG analysis
    %
    % This function reads and parses JSON files for excluded channels and region
    % definitions, then adds them to the configuration structure.
    %
    % Input:
    %   config - Configuration structure that contains paths to JSON files
    %
    % Output:
    %   parsedConfig - Updated configuration with excluded channels and regions
    
    % Initialize output as a copy of input
    parsedConfig = config;
    
    % Read excluded channels from JSON file if it exists
    if isfield(config, 'exclude_json') && exist(config.exclude_json, 'file')
        try
            excluded_channels = readExcludeChannelsFromJSON(config.exclude_json);
            parsedConfig.exclude_channels = excluded_channels;
            fprintf('Loaded %d excluded channels from %s\n', length(excluded_channels), config.exclude_json);
        catch err
            warning('Error reading exclude.json: %s', err.message);
            fprintf('Using default empty exclude list\n');
            parsedConfig.exclude_channels = {};
        end
    else
        if isfield(config, 'exclude_json')
            fprintf('No exclude.json found at %s. Using empty exclude list.\n', config.exclude_json);
        else
            fprintf('No exclude_json path specified. Using empty exclude list.\n');
        end
        parsedConfig.exclude_channels = {};
    end
    
    % Read region definitions from JSON file if it exists
    if isfield(config, 'regions_json') && exist(config.regions_json, 'file')
        try
            regions = readRegionsFromJSON(config.regions_json);
            parsedConfig.regions = regions;
            fprintf('Loaded %d region definitions from %s\n', length(fieldnames(regions)), config.regions_json);
        catch err
            warning('Error reading regions.json: %s', err.message);
            fprintf('No region-specific analysis will be performed\n');
            parsedConfig.regions = struct();
        end
    else
        if isfield(config, 'regions_json')
            fprintf('No regions.json found at %s. No region-specific analysis will be performed.\n', config.regions_json);
        else
            fprintf('No regions_json path specified. No region-specific analysis will be performed.\n');
        end
        parsedConfig.regions = struct();
    end
end

% Helper function to read excluded channels from JSON
function excluded_channels = readExcludeChannelsFromJSON(json_file)
    % Read the entire file content
    fid = fopen(json_file, 'r');
    if fid == -1
        error('Could not open file: %s', json_file);
    end
    raw_json = fread(fid, '*char')';
    fclose(fid);
    
    % Parse JSON content
    try
        % If the file contains a simple array of channel names
        if startsWith(strtrim(raw_json), '[')
            % Parse as a simple array
            channels_cell = jsondecode(raw_json);
            excluded_channels = cell(size(channels_cell));
            for i = 1:length(channels_cell)
                excluded_channels{i} = channels_cell{i};
            end
        else
            % Parse as a JSON object and extract channels
            json_data = jsondecode(raw_json);
            
            % If it's a structure, assume channels are in a field or need to be concatenated
            % This is a simplified approach - modify as needed based on your JSON structure
            excluded_channels = {};
            
            % Check if it's a direct array of channels
            if isfield(json_data, 'channels')
                excluded_channels = json_data.channels;
            else
                % If it's a list of channels, concatenate them
                excluded_channels = {};
            end
        end
    catch err
        error('Error parsing JSON: %s', err.message);
    end
end

% Helper function to read region definitions from JSON
function regions = readRegionsFromJSON(json_file)
    % Read the entire file content
    fid = fopen(json_file, 'r');
    if fid == -1
        error('Could not open file: %s', json_file);
    end
    raw_json = fread(fid, '*char')';
    fclose(fid);
    
    % Parse JSON content
    try
        json_data = jsondecode(raw_json);
        
        % Convert JSON structure to MATLAB structure
        regions = struct();
        
        % Get field names from JSON object
        region_names = fieldnames(json_data);
        
        % Process each region
        for i = 1:length(region_names)
            region_name = region_names{i};
            region_channels = json_data.(region_name);
            
            % Handle different possible formats of the channel list
            if iscell(region_channels) && length(region_channels) == 1 && iscell(region_channels{1})
                % Nested array format: {"region": [[ch1, ch2, ...]]}
                channels = region_channels{1};
                fprintf('Region %s: Extracted %d channels from nested array\n', ...
                    region_name, length(channels));
            elseif iscell(region_channels)
                % Direct array format: {"region": [ch1, ch2, ...]}
                channels = region_channels;
                fprintf('Region %s: Using %d direct channel entries\n', ...
                    region_name, length(channels));
            else
                % Some other unexpected format
                warning('Unexpected format for region %s. Skipping.', region_name);
                continue;
            end
            
            % Store channels in regions structure using the region_name directly
            % as it should already be in snake_case format in the JSON
            regions.(region_name) = channels;
            
            % Print first few channels for verification
            if length(channels) > 3
                fprintf('First 3 channels: %s, %s, %s\n', ...
                    channels{1}, channels{2}, channels{3});
            else
                first_channels = strjoin(channels, ', ');
                fprintf('Channels: %s\n', first_channels);
            end
        end
    catch err
        error('Error parsing JSON: %s', err.message);
    end
end
