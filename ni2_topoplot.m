function ni2_topoplot(sens, topo, cfg)

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
  error('the number of data points in the topography is different from the number of channels');
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
