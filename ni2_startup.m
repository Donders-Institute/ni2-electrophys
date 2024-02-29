function ni2_startup

% NI2_STARTUP helps to get the path-settings correct for doing the MATLAB-exercises
% for the Source Reconstruction part of the Neuroimaging2 electrophysiology course.
%
% This expects the code to be organized as
%  somewhere/ni2
%  somewhere/fieldtrip

fprintf('clearing all variables from the workspace\n');
evalin('base', 'clear')

ni2_dir = fileparts(mfilename('fullpath')); % this file itself is contained in the ni2 directory
ft_dir  = fileparts(which('ft_defaults'));

if ~exist(ft_dir, 'dir')
  % see whether we can find it next to the ni2 directory
  d = dir(fileparts(ni2_dir));
  select = find(startsWith({d.name}, 'fieldtrip'), 1);
  ft_dir = fullfile(fileparts(ni2_dir), d(select).name);
  if ~exist(fullfile(ft_dir, 'ft_defaults.m'), 'file')
    % no luck
    ft_dir = [];
  end
end

if ~exist(ft_dir, 'dir')
  warning('cannot locate fieldtrip directory')
end

fprintf('restoring the default MATLAB-path\n');
restoredefaultpath;

% add the ni2 and fieldtrip directories to the path
fprintf('adding % s to the MATLAB-path\n', ni2_dir);
addpath(ni2_dir);

fprintf('adding % s to the MATLAB-path\n', ft_dir);
addpath(ft_dir);

try
  % add the necessary subdirectories of FieldTrip to the path
  ft_defaults;
catch
  warning('cannot add the necessary subdirectories of FieldTrip to the path')
end