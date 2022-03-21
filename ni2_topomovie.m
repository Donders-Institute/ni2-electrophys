function ni2_topomovie(sens, dat, time, cfg)

% ensure the sensor description to be according to what FieldTrip expects
sens = ft_datatype_sens(sens);

% construct a timelocked ERP data structure
timelock = [];
timelock.time = time;
timelock.avg = dat;
timelock.label = sens.label;

if nargin<4
  cfg = [];
end

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