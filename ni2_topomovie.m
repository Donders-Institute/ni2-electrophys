function ni2_topomovie(sens, data, time, cfg)

% NI2_TOPOMOVIE makes a movie of the 2D-projection of the spatial EEG or MEG
% topography as it changes over time.
%
% Use as
%   ni2_topoplot(sens, data, time)
%
% Input arguments:
%   sens = a sensor-array as obtained with NI2_SENSORS
%   data = a NxT matrix representing a spatial topography over time. The number of rows (N) 
%          should be equal to the number of sensors in the sens structure, the number of 
%          columns should be equal to the number of time points.

if nargin<4
  cfg = [];
end

% ensure the sensor description to be according to what FieldTrip expects
sens = ft_datatype_sens(sens);

% construct a timelocked ERP data structure that FieldTrip can deal with
timelock = [];
timelock.time = time;
timelock.avg = data;
timelock.label = sens.label;

if ~isfield(cfg, 'layout')
  % create layout
  if strncmp(sens.chantype{1}, 'meg', 3)
    tmpcfg.grad = sens;
  else
    tmpcfg.elec = sens;
  end
  cfg.layout = ft_prepare_layout(tmpcfg);
end

cfg.speed = 20;
cfg.zlim = 'maxabs';
cfg.colormap = '*RdBu';
ft_movieplotER(cfg, timelock);

% update the figure immediately
drawnow