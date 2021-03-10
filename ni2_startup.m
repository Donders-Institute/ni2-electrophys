function ni2_startup

% NI2_STARTUP ensures that the path-settings are correct when doing the
% MATLAB-exercises for the Neuroimaging2 course.

fprintf('clearing all variables from the workspace\n');
clear all

fprintf('restoring the default MATLAB-path\n');
restoredefaultpath;

curr_dir = pwd;
ft_dir   = [curr_dir(1:(strfind(curr_dir,'ni2')-1)) filesep 'fieldtrip'];

% add the ni2 and fieldtrip directories to the path
fprintf('adding % s to the MATLAB-path\n', curr_dir);
addpath(curr_dir);

fprintf('adding % s to the MATLAB-path\n', ft_dir);
addpath(ft_dir);

% add the necessary subdirectories of FieldTrip to the path
ft_defaults;
