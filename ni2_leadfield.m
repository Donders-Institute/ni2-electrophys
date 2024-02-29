function leadfield = ni2_leadfield(sens, headmodel, dippar)

% NI2_LEADFIELD generates a leadfield for a given source (or set of
% sources) using a specified sensor array and volume conductor model.
%
% Use as
%   leadfield = ni2_leadfield(sens, headmodel, dippar)
%
% Input arguments:
%   sens      = a sensor array, obtained with NI2_SENSORS
%   headmodel = a volume conductor model, obtained with NI2_HEADMODEL
%   dippar    = Nx3 or Nx6 matrix with dipole parameters, where N is the number of dipoles
%
% Each row in the dippar-matrix represents a source, the first columns 1-3 are the
% dipole position parameters (x,y,z coordinates in Cartesian space), and the optional
% columns 4-6 are the dipole moment parameters (x,y,z orientation, these can include
% the amplitude).
%
% Output argument:
%   leadfield = MxN or Mx(Nx3) matrix with the leadfield, where M is the number of sensors
%               and N is the number of dipoles. There are three columns for dipoles with
%               an unspecified orientation.

ndip = size(dippar, 1);

cfg = [];

if strcmp(sens.chantype{1}, 'eeg')
  cfg.elec = sens;
else
  cfg.grad = sens;
end

cfg.headmodel          = headmodel;
cfg.sourcemodel.pos    = dippar(:,1:3);
cfg.sourcemodel.inside = 1:ndip;
cfg.sourcemodel.unit   = headmodel.unit;
cfg.reducerank         = 'no';

if strcmp(headmodel.type, 'singleshell')
  cfg.singleshell.batchsize = 2500;
end

leadf = ft_prepare_leadfield(cfg);

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
