function leadfield = ni2_leadfield(sens, headmodel, dippar)

% NI2_LEADFIELD generates a leadfield for a given source (or set of
% sources) using a specified sensor array and volume conductor model.
%
% Use as:
%   leadfield = ni2_leadfield(sens, headmodel, dippar)
%
% Input arguments:
%   - sens      = a sensor array, obtained with ni2_sensors
%   - headmodel = a volume conductor model, obtained with ni2_headmodel
%   - dippar    = Nx6, or Nx3 matrix with dipole parameters.
%
% Each row in the dippar-matrix represents a source, the first 3 columns
% are the position parameters (x,y,z coordinates in Cartesian space), and
% the options 4-6 columns are the dipole moment parameters (x,y,z
% orientation and amplitude).
% 
% Output argument:
%   - leadfield = MxN or Mx(Nx3) matrix with the leadfield, where M is the
%                 number of sensors, and N the number of sources

ndip = size(dippar, 1);

cfg = [];
if strcmp(sens.chantype{1}, 'eeg')
  cfg.elec = sens;
else
  cfg.grad = sens;
end
cfg.headmodel   = headmodel;
cfg.grid.pos    = dippar(:,1:3);
cfg.grid.inside = 1:ndip;
cfg.grid.unit   = 'cm';
cfg.reducerank  = 'no';
leadf           = ft_prepare_leadfield(cfg);

if size(dippar, 2)==6
  leadfield = zeros(size(leadf.leadfield{1},1), ndip);
  for k = 1:ndip
    leadfield(:, k) = leadf.leadfield{k}*dippar(k,4:6)';
  end
else
  leadfield = zeros(size(leadf.leadfield{1},1), ndip*3);
  for k = 1:ndip
    indx = (k-1)*3 + (1:3);
    leadfield(:, indx) = leadf.leadfield{k};
  end
end
