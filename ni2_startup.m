function ni2_startup

% NI2_STARTUP helps to get the path-settings correct for doing the
% MATLAB-exercises for the Neuroimaging2 course.

fprintf('clearing all variables from the workspace\n');
evalin('base', 'clear')

fprintf('restoring the default MATLAB-path\n');
restoredefaultpath;

ni2_dir = fileparts(mfilename('fullpath'));
ft_dir  = fullfile(fileparts(ni2_dir), 'fieldtrip');

% add the ni2 and fieldtrip directories to the path
fprintf('adding % s to the MATLAB-path\n', ni2_dir);
addpath(ni2_dir);

fprintf('adding % s to the MATLAB-path\n', ft_dir);
addpath(ft_dir);

% add the necessary subdirectories of FieldTrip to the path
ft_defaults;
