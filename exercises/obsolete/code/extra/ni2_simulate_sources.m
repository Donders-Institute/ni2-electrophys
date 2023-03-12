%% noiseless case

[data1,time1] = ni2_activation;
[data2,time2] = ni2_activation('frequency',11,'latency',0.48);
sens = ni2_sensors('type', 'meg');
headmodel = ni2_headmodel('type','spherical','nshell',1);
leadfield = ni2_leadfield(sens,headmodel,[4.9 0 6.2 0 1 0]); % position 2351 in grid
leadfield2 = ni2_leadfield(sens,headmodel,[-5.3 0 5.9 1 0 0]); % position 2342 in grid
sensordata = leadfield*data1+leadfield2*data2;%+randn(301,1000)*.7e-10;

data        = [];
data.avg    = sensordata;
data.time   = time1;
data.label  = sens.label;
data.grad   = sens;
data.cov    = eye(numel(sens.label))*((.7e-9).^2);
data.dimord = 'chan_time';

sourcemodel = ni2_sourcemodel('type','grid','resolution',1);

cfg                    = [];
cfg.sourcemodel        = sourcemodel;
cfg.headmodel          = headmodel;
cfg.method             = 'mne';
cfg.mne.prewhiten      = 'yes';
cfg.mne.scalesourcecov = 'yes';
cfg.mne.lambda         = 0;
cfg.keepleadfield      = 'yes';
source = ft_sourceanalysis(cfg, data);

%% with noise

[data1,time1] = ni2_activation;
[data2,time2] = ni2_activation('frequency',11,'latency',0.48);
sens = ni2_sensors('type', 'meg');
headmodel = ni2_headmodel('type','spherical','nshell',1);
leadfield = ni2_leadfield(sens,headmodel,[4.9 0 6.2 0 1 0]); % position 2351 in grid
leadfield2 = ni2_leadfield(sens,headmodel,[-5.3 0 5.9 1 0 0]); % position 2342 in grid
sensordata = leadfield*data1+leadfield2*data2+randn(301,1000)*.7e-10;

data        = [];
data.avg    = sensordata;
data.time   = time1;
data.label  = sens.label;
data.grad   = sens;
data.cov    = cov(randn(301,1000)'*.7e-10);
data.dimord = 'chan_time';

sourcemodel = ni2_sourcemodel('type','grid','resolution',1);

cfg                    = [];
cfg.sourcemodel        = sourcemodel;
cfg.headmodel          = headmodel;
cfg.method             = 'mne';
cfg.mne.prewhiten      = 'yes';
cfg.mne.scalesourcecov = 'yes';
cfg.mne.lambda         = 0;
cfg.keepleadfield      = 'yes';
source_noise = ft_sourceanalysis(cfg, data);

%% with noise regularisation
cfg.mne.lambda = 0.5;
source_noise_reg = ft_sourceanalysis(cfg, data);


%% with depth weighting
[data1,time1] = ni2_activation;
[data2,time2] = ni2_activation('frequency',11,'latency',0.48);
sens = ni2_sensors('type', 'meg');
headmodel = ni2_headmodel('type','spherical','nshell',1);
leadfield = ni2_leadfield(sens,headmodel,[4.9 0 6.2 0 1 0]); % position 2351 in grid
leadfield2 = ni2_leadfield(sens,headmodel,[-5.3 0 5.9 1 0 0]); % position 2342 in grid
sensordata = leadfield*data1+leadfield2*data2+randn(301,1000)*2e-10;

data        = [];
data.avg    = sensordata;
data.time   = time1;
data.label  = sens.label;
data.grad   = sens;
data.cov    = cov(randn(301,1000)'*2e-10);
data.dimord = 'chan_time';

sourcemodel = ni2_sourcemodel('type','grid','resolution',1);

cfg                    = [];
cfg.sourcemodel        = sourcemodel;
cfg.headmodel          = headmodel;
cfg.method             = 'mne';
cfg.mne.prewhiten      = 'yes';
cfg.mne.scalesourcecov = 'yes';
cfg.keepleadfield      = 'yes';
cfg.mne.lambda         = 0.5;
source_noise = ft_sourceanalysis(cfg, data);

cfg.normalize = 'yes';
cfg.normalizeparam = 1;
source_noise_lfnorm = ft_sourceanalysis(cfg, data);
