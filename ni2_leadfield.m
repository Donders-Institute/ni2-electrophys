function leadfield = ni2_leadfield(sens, headmodel, dippar)

ndip = size(dippar, 1);

cfg = [];
if strcmp(sens.chantype{1}, 'eeg')
  cfg.elec = sens;
else
  cfg.grad = sens;
end
cfg.vol         = headmodel;
cfg.grid.pos    = dippar(:,1:3);
cfg.grid.inside = 1:ndip;
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
