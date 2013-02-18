function ni2_topoplot(sens, topo, cfg)

% NI2_TOPOPLOT plots a 2D-projection of a spatial topography.
%
% Use as:
%   ni2_topoplot(sens, topo, cfg)
%
% Input arguments:
%   - sens = a sensor-array obtained with ni2_sensors
%   - topo = a Nx1 vector representing a spatial topography. The number of
%   elements should be equal to the number of sensors in the sens structure
%   - cfg  = a FieldTrip style configuration that specifies optional
%   parameters that influence the plotting (see FT_TOPOPLOTER for more
%   info).

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
  %error('the number of data points in the topography is different from the number of channels');
end

% create a data structure that FieldTrip can deal with
data        = [];
data.avg    = topo;
data.label  = sens.label;
data.time   = 1;
data.dimord = 'chan_time';

if nargin<3
  cfg = [];
end

if ~isfield(cfg, 'layout')
  % create layout
  if strncmp(sens.chantype{1}, 'meg', 3);
    tmpcfg.grad = sens;
  else
    tmpcfg.elec = sens;
  end
  cfg.layout = ft_prepare_layout(tmpcfg);
end

figure;ft_topoplotER(cfg, data);
