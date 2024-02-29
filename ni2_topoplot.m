function ni2_topoplot(sens, topo, cfg)

% NI2_TOPOPLOT plots a 2D-projection of the spatial topography of the EEG or MEG
% distribution over the head.
%
% Use as
%   ni2_topoplot(sens, topo)
%
% Input arguments:
%   sens = a sensor-array as obtained with NI2_SENSORS
%   topo = a Nx1 vector representing a spatial topography. The number of elements 
%          or rows (N) should be equal to the number of sensors in the sens structure.

if nargin<3
  cfg = [];
  cfg.interactive = 'no';
end

% ensure the sensor description to be according to what FieldTrip expects
sens = ft_datatype_sens(sens);

% do a sanity check on the topo
if all(size(topo)>1)
  error('topo should be either a column or row vector');
end

% make column vector
topo = topo(:);

% do another sanity check on the topo
if numel(topo)~=size(sens.chanpos,1)
  % error('the number of data points in the topography is different from the number of channels');
end

% construct a timelocked ERP data structure that FieldTrip can deal with
data        = [];
data.avg    = topo;
data.label  = sens.label;
data.time   = 1;
data.dimord = 'chan_time';

if ~isfield(cfg, 'layout')
  % create layout
  if strncmp(sens.chantype{1}, 'meg', 3)
    tmpcfg.grad = sens;
  else
    tmpcfg.elec = sens;
  end
  cfg.layout = ft_prepare_layout(tmpcfg);
end

cfg.colormap = '*RdBu';
cfg.zlim = 'maxabs';
ft_topoplotER(cfg, data);

% update the figure immediately
drawnow
