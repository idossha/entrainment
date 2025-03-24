

clear all; close all;
experiment_path = '/Volumes/Ido/analyze';
eeglab_path = '/Users/idohaber/Documents/MATLAB/eeglab2024.0/';
addpath(eeglab_path);
eeglab nogui;
electrode_file = '../assets/GSN-HydroCel-256.sfp';
exclude_file = '../assets/exclude.json';

% Define the function
function visualizeChannel(electrode_file, exclude_file)
  % Load electrode locations
  chanlocs = readlocs(electrode_file);
  
  % Read excluded channels
  exclude_data = jsondecode(fileread(exclude_file));
  
  % Fix: Access the exclude_channels field correctly
  % The error was accessing it with dot notation
  exclude_channels = exclude_data;
  if isstruct(exclude_data) && isfield(exclude_data, 'exclude_channels')
      exclude_channels = exclude_data.exclude_channels;
  end
  
  % Fix: Typo in "lengr" should be "length"
  data = ones(length(chanlocs), 1);
  
  for i = 1:length(chanlocs)
    if ismember(chanlocs(i).labels, exclude_channels)
      data(i) = 0;
    end
  end
  
  figure;
  topoplot(data, chanlocs, 'electrodes', 'labels', 'style', 'map');
  title('Channels (0 = excluded, 1 = included)');
  caxis([-1 1]);
end

% Call the function
visualizeChannel(electrode_file, exclude_file);
