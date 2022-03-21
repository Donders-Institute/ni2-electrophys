L=[1 0 1; 0 1 1];
y=[2 1]';

Lpinv = L'*(L*L')^-1;

Lpinv = pinv(L);

% inv(L)

L*Lpinv;
Lpinv*L;


[data1,time1]=ni2_activation;
[data2,time2]=ni2_activation('frequency',11,'latency',0.48);
% sens = ni2_sensors('type','meg');
load('ni2_megsensors.mat'); sens = sensmeg;
headmodel = ni2_headmodel('type','spherical','nshell',1);
leadfield1 = ni2_leadfield(sens,headmodel,[4.9 0 6.2 0 1 0]); % position 2352 in grid
leadfield2 = ni2_leadfield(sens,headmodel,[-5.3 0 5.9 1 0 0]); % position 2342 in grid
sensordata = leadfield1*data1+leadfield2*data2;

data        = [];
data.avg    = sensordata;
data.time   = time1;
data.label  = sens.label;
data.grad   = sens;
data.cov    = eye(numel(sens.label));
data.dimord = 'chan_time';

sourcemodel = ni2_sourcemodel('type','grid','resolution',1);

cfg                    = [];
cfg.grid               = sourcemodel;
cfg.headmodel          = headmodel;
cfg.method             = 'mne';
cfg.mne.prewhiten      = 'yes';
cfg.mne.scalesourcecov = 'yes';
cfg.mne.lambda         = 0;
cfg.keepleadfield      = 'yes';
source = ft_sourceanalysis(cfg, data);

figure; plot(source.time,source.avg.mom{2352}); legend({'x' 'y' 'z'});
figure; plot(source.time,source.avg.mom{2342}); legend({'x' 'y' 'z'});
figure; plot(source.time,source.avg.mom{2347}); legend({'x' 'y' 'z'});
figure; plot(source.time,source.avg.mom{2713}); legend({'x' 'y' 'z'});


[data1,time1]=ni2_activation;
[data2,time2]=ni2_activation('frequency',11,'latency',0.48);
% sens = ni2_sensors('type','meg');
load('ni2_megsensors.mat'); sens = sensmeg;
headmodel = ni2_headmodel('type','spherical','nshell',1);
leadfield1 = ni2_leadfield(sens,headmodel,[4.9 0 6.2 0 1 0]); % position 2352 in grid
leadfield2 = ni2_leadfield(sens,headmodel,[-5.3 0 5.9 1 0 0]); % position 2342 in grid
sensordata = leadfield1*data1+leadfield2*data2+randn(301,1000)*.7e-10;

data        = [];
data.avg    = sensordata;
data.time   = time1;
data.label  = sens.label;
data.grad   = sens;
data.cov    = cov(randn(301,1000)'*.7e-10);
data.dimord = 'chan_time';

sourcemodel = ni2_sourcemodel('type','grid','resolution',1);

cfg                    = [];
cfg.grid               = sourcemodel;
cfg.headmodel          = headmodel;
cfg.method             = 'mne';
cfg.mne.prewhiten      = 'yes';
cfg.mne.scalesourcecov = 'yes';
cfg.mne.lambda         = 0;
cfg.keepleadfield      = 'yes';
source_noise = ft_sourceanalysis(cfg, data);

figure; plot(source_noise.time,source_noise.avg.mom{2352},'linewidth', 2); legend({'x' 'y' 'z'});
figure; plot(source_noise.time,source_noise.avg.mom{2342},'linewidth', 2); legend({'x' 'y' 'z'});
figure; plot(source_noise.time,source_noise.avg.mom{2347},'linewidth', 2); legend({'x' 'y' 'z'});
figure; plot(source_noise.time,source_noise.avg.mom{2713},'linewidth', 2); legend({'x' 'y' 'z'});

cfg                    = [];
cfg.grid               = sourcemodel;
cfg.headmodel          = headmodel;
cfg.method             = 'mne';
cfg.mne.prewhiten      = 'yes';
cfg.mne.scalesourcecov = 'yes';
cfg.mne.lambda         = 2;
cfg.keepleadfield      = 'yes';
source_noise_reg = ft_sourceanalysis(cfg, data);

figure; plot(source_noise_reg.time,source_noise_reg.avg.mom{2352},'linewidth', 2); legend({'x' 'y' 'z'});
figure; plot(source_noise_reg.time,source_noise_reg.avg.mom{2342},'linewidth', 2); legend({'x' 'y' 'z'});
figure; plot(source_noise_reg.time,source_noise_reg.avg.mom{2347},'linewidth', 2); legend({'x' 'y' 'z'});
figure; plot(source_noise_reg.time,source_noise_reg.avg.mom{2713},'linewidth', 2); legend({'x' 'y' 'z'});

L = cat(2,source_noise_reg.leadfield{source_noise_reg.inside});
S = cat(1,source_noise_reg.avg.mom{source_noise_reg.inside});
model = L*S;
residual = sensordata-model;
figure; plot(residual');

cfg                    = [];
cfg.grid               = sourcemodel;
cfg.headmodel          = headmodel;
cfg.method             = 'mne';
cfg.mne.prewhiten      = 'yes';
cfg.mne.scalesourcecov = 'yes';
cfg.mne.lambda         = 0.01;
cfg.keepleadfield      = 'yes';
source_noise_reg = ft_sourceanalysis(cfg, data);

figure; plot(source_noise_reg.time,source_noise_reg.avg.mom{2352}); legend({'x' 'y' 'z'});
figure; plot(source_noise_reg.time,source_noise_reg.avg.mom{2342}); legend({'x' 'y' 'z'});
figure; plot(source_noise_reg.time,source_noise_reg.avg.mom{2347}); legend({'x' 'y' 'z'});
figure; plot(source_noise_reg.time,source_noise_reg.avg.mom{2713}); legend({'x' 'y' 'z'});

L = cat(2,source_noise_reg.leadfield{source_noise_reg.inside});
S = cat(1,source_noise_reg.avg.mom{source_noise_reg.inside});
model = L*S;
residual = sensordata-model;
figure; plot(residual');

cfg                    = [];
cfg.grid               = sourcemodel;
cfg.headmodel          = headmodel;
cfg.method             = 'mne';
cfg.mne.prewhiten      = 'yes';
cfg.mne.scalesourcecov = 'yes';
cfg.mne.lambda         = 0.5;
cfg.keepleadfield      = 'yes';
cfg.normalize          = 'yes';
cfg.normalizeparam     = 1;
source_noise_reg    = ft_sourceanalysis(cfg, data);

figure; plot(source_noise_reg.time,source_noise_reg.avg.mom{2352}); legend({'x' 'y' 'z'});
% figure; plot(source_noise_reg.time,source_noise_reg.avg.mom{2342}); legend({'x' 'y' 'z'});
% figure; plot(source_noise_reg.time,source_noise_reg.avg.mom{2347}); legend({'x' 'y' 'z'});
figure; plot(source_noise_reg.time,source_noise_reg.avg.mom{2713}); legend({'x' 'y' 'z'});

figure; plot(source_noise_lfnorm.time,source_noise_lfnorm.avg.mom{2352}); legend({'x' 'y' 'z'});
% figure; plot(source_noise_lfnorm.time,source_noise_lfnorm.avg.mom{2342}); legend({'x' 'y' 'z'});
% figure; plot(source_noise_lfnorm.time,source_noise_lfnorm.avg.mom{2347}); legend({'x' 'y' 'z'});
figure; plot(source_noise_lfnorm.time,source_noise_lfnorm.avg.mom{2713}); legend({'x' 'y' 'z'});
