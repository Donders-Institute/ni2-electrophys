% create some data
[S,time] = ni2_activation;
sens        = ni2_sensors('type','eeg','n',64);

% in order to fool around with the reference of the leadfield, we need to
% add an explicit 'tra', otherwise ft_compute_leadfield will subtract the
% mean by default. By adding a tra as Identity matrix, the leadfields will
% be returned 'without' reference. Here we add a reference on the left, 'eeg74'
ref            = 25;
sens.tra       = eye(numel(sens.label));
sens.tra(:,ref) = sens.tra(:,ref)-1;

% create a volume conductor model
headmodel   = ni2_headmodel('type', 'spherical', 'nshell', 3);

% create a 3D grid for the initial gridsearch
sourcemodel = ni2_sourcemodel('type','grid','resolution',1);

% simulate some sensor level data
dippar      = [-8./sqrt(2) 0.5 8./sqrt(2) [1 1 1].*(3^(-1/3))];
leadfield   = ni2_leadfield(sens, headmodel, dippar);
noise       = sens.tra*randn(numel(sens.label),1000)*1e-3;
sensordata  = leadfield*S+noise;
assert(all(sensordata(ref,:)==0));

% visualize it
topo_observed = sensordata(:,500);
ni2_topoplot(sens, topo_observed);

data        = [];
data.avg    = sensordata;
data.time   = time;
data.label  = sens.label;
data.elec   = sens;
data.dimord ='chan_time';

% this is what I think Sarang means, i.e. data with an implicitref present,
% and then svd'ed in order to remove the last component
[u,s,v]  = svd(data.avg, 'econ');
montage2.tra = u(:,1:end-1)';
montage2.labelold = data.label;
for k = 1:size(montage2.tra,1)
  montage2.labelnew{k,1} = sprintf('comp% 02d',k);
end
data_svd      = ft_apply_montage(data, montage2);
data_svd.elec = ft_apply_montage(data.elec, montage2);

% this is the data with the reference channel removed
cfg         = [];
cfg.channel = data.label;
cfg.channel(ref) = [];
data        = ft_selectdata(cfg, data);

% do some referencing by hand to obtain the average referenced data
montage1.tra  = eye(numel(data.label))- ones(numel(data.label))./numel(data.label);
montage1.labelold = data.label;
montage1.labelnew = data.label;
data_avg      = ft_apply_montage(data, montage1);
data_avg.elec = ft_apply_montage(data.elec, montage1);

cfg            = [];
cfg.gridsearch = 'yes';
cfg.model      = 'regional';
cfg.headmodel  = headmodel;
cfg.sourcemodel = sourcemodel;
cfg.latency    = [0.49 0.51];
cfg.nonlinear  = 'yes';
cfg.numdipoles = 1;
dip            = ft_dipolefitting(cfg, data);
dip_avg        = ft_dipolefitting(cfg, data_avg);
dip_svd        = ft_dipolefitting(cfg, data_svd);

% ni2_topoplot(sens,dip.Vmodel(:,6));
% ni2_topoplot(sens,dip.Vdata(:,6));

sens_jittered = ni2_sensors('type', 'eeg', 'jitter', 0.1, 'n', 64);

dataj = data;
dataj_avg = data_avg;
dataj_svd = data_svd;

dataj.elec.elecpos = sens_jittered.chanpos;
dipj = ft_dipolefitting(cfg, dataj);
dataj_avg.elec.elecpos = sens_jittered.chanpos;
dipj_avg = ft_dipolefitting(cfg, dataj_avg);
dataj_svd.elec.elecpos = sens_jittered.chanpos;
dipj_svd = ft_dipolefitting(cfg, dataj_svd);
