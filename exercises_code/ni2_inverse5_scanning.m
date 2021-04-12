% dealing with the concept of covariance

% create some random vectors
x = randn(1,100);
x = x-mean(x);
y = randn(1,100);
y = y-mean(y);

covxy = (x*y')/99;
figure; plot(x); hold on; plot(y,'r');
figure; plot(x,y,'.');

c = randn(1,100);
x = randn(1,100)+c;
x = x-mean(x);
y = randn(1,100)+c;
y = y-mean(y);

covxy = (x*y')/99
figure; plot(x,y,'.');
figure; plot(x); hold on; plot(y,'r');

c = randn(1,100);
x = randn(1,100)+2.*c;
x = x-mean(x);
y = randn(1,100)+2.*c;
y = y-mean(y);

covxy = (x*y')/99
figure; plot(x,y,'.');
figure; plot(x); hold on; plot(y,'r');

c = randn(1,100);
x = randn(1,100)+2.*c;
x = x-mean(x);
y = randn(1,100)-2.*c;
y = y-mean(y);

covxy = (x*y')/99
figure; plot(x,y,'.');
figure; plot(x); hold on; plot(y,'r');

% compute the covariance matrix from some simulated data in 2 different
% ways
[data,time]=ni2_activation;
[data2,time]=ni2_activation('frequency',11,'latency',0.48);
sens = ni2_sensors('type','eeg');
headmodel = ni2_headmodel('type','spherical','nshell',3);
leadfield = ni2_leadfield(sens,headmodel,[4.9 0 6.2 0 1 0]);
leadfield2 = ni2_leadfield(sens,headmodel,[-5.3 0 5.9 1 0 0]);

% way 1
sensordata = leadfield*data+leadfield2*data2;
C1 = sensordata*sensordata';

% way 2
sourcecov=[data; data2]*[data; data2]';
C2 = leadfield*leadfield'   * sourcecov(1,1)   + ...
   leadfield*leadfield2'  * sourcecov(1,2)   + ...
   leadfield2*leadfield'  * sourcecov(2,1)   + ...
   leadfield2*leadfield2' * sourcecov(2,2);
 
 

%% do the beamforming on some simulated data

% create a leadfield for a few locations
sens = ni2_sensors('type','meg');

headmodel = ni2_headmodel('type','spherical','nshell',1);
leadfield1 = ni2_leadfield(sens,headmodel,[4.9 0 6.2 0 1 0]); % position 2352 in grid
leadfield2 = ni2_leadfield(sens,headmodel,[-5.3 0 5.9 1 0 0]); % position 2342 in grid
leadfield3 = ni2_leadfield(sens,headmodel,[4.2 -2 7 0 0.2 0.7]); % 2674
leadfield4 = ni2_leadfield(sens,headmodel,[-2 -7.8 3 -1 0 0]); % 1110

% create the time course of activation
[s1,t1]=ni2_activation('latency',.45,'frequency',3);
[s2,t2]=ni2_activation('latency',.5);
[s3,t3]=ni2_activation('latency',.55,'frequency',30);
[s4,t4]=ni2_activation('latency',.5, 'frequency',15);

% create the sensor data
sensordata = leadfield1*s1+leadfield2*s2+leadfield3*s3+leadfield4*s4+randn(301,1000)*0.04e-8;

% create a 3D grid
sourcemodel = ni2_sourcemodel('type','grid','resolution',1);

% use FieldTrip to quickly compute the forward model
cfg = [];
cfg.grid = sourcemodel;
cfg.grad = sens;
cfg.vol  = headmodel;
% cfg.normalize = 'yes';
sourcemodel = ft_prepare_leadfield(cfg);
L = cat(2,sourcemodel.leadfield{sourcemodel.inside});

% compute the covariance
C  = cov(sensordata');
iC = pinv(C); % so that we compute it only once
iCr = inv(C+eye(301)*5e-19);

% compute the beamformer spatial filter
for ii = 1:size(L,2)/3
  indx=(ii-1)*3+(1:3);
  Lr = L(:,indx);  % Lr is the leadfield for source r
  wbfr(indx,:)=pinv(Lr'*iCr*Lr)*Lr'*iCr;
  wbf(indx,:)=pinv(Lr'*iC*Lr)*Lr'*iC;
end
sbf = wbf*sensordata;
sbfr = wbfr*sensordata;

wmn = L'*inv(L*L'+eye(301)*1e-15);
smn = wmn*sensordata;

sel = find(ismember(find(sourcemodel.inside), [1110 2342 2352 2674]));
sel = repmat((sel-1)*3,1,3)+repmat(1:3,numel(sel),1);

figure;
% subplot(1,2,1); plot(t1,sbf(sel(1,:),:));
subplot(1,2,1); plot(t1,sbfr(sel(1,:),:));
subplot(1,2,2); plot(t1,s4);
% subplot(1,2,2); plot(t1,smn(sel(1,:),:));

figure;
subplot(1,2,1); plot(t1,sbf(sel(2,:),:));
subplot(1,2,2); plot(t1,s2);
% subplot(1,2,2); plot(t1,smn(sel(2,:),:));

figure;
subplot(1,2,1); plot(t1,sbfr(sel(3,:),:));
subplot(1,2,2); plot(t1,s1);
% subplot(1,2,2); plot(t1,smn(sel(3,:),:));

figure;
subplot(1,2,1); plot(t1,sbf(sel(4,:),:));
subplot(1,2,2); plot(t1,s3);
% subplot(1,2,2); plot(t1,smn(sel(4,:),:));


% compute the beamformer spatial filter with limited data available
R10  = cov(sensordata(:,491:500)');
iR10 = inv(R10+eye(301)*0.1e-25);
% compute the beamformer spatial filter
for ii = 1:size(L,2)/3
  indx=(ii-1)*3+(1:3);
  Lr = L(:,indx); % Lr is the leadfield for source r
  wbf10(indx,:)=pinv(Lr'*iR10*Lr)*Lr'*iR10;
end
sbf10 = wbf10*sensordata;

figure;
subplot(2,2,1); plot(t1,sbf(sel(1,:),:)); xlabel('sbf');
subplot(2,2,2); plot(t1,sbfr(sel(1,:),:)); xlabel('sbfr');
subplot(2,2,3); plot(t1,sbf10(sel(1,:),:)); xlabel('sbf10');
subplot(2,2,4); plot(t1,s4); xlabel('simulated source 4');

figure;
subplot(2,2,1); plot(t1,sbf(sel(2,:),:)); xlabel('sbf');
subplot(2,2,2); plot(t1,sbfr(sel(2,:),:)); xlabel('sbfr');
subplot(2,2,3); plot(t1,sbf10(sel(2,:),:)); xlabel('sbf10');
subplot(2,2,4); plot(t1,s2); xlabel('simulated source 2');

figure;
subplot(2,2,1); plot(t1,sbf(sel(3,:),:)); xlabel('sbf');
subplot(2,2,2); plot(t1,sbfr(sel(3,:),:)); xlabel('sbfr');
subplot(2,2,3); plot(t1,sbf10(sel(3,:),:)); xlabel('sbf10');
subplot(2,2,4); plot(t1,s1); xlabel('simulated source 1');

figure;
subplot(2,2,1); plot(t1,sbf(sel(4,:),:)); xlabel('sbf');
subplot(2,2,2); plot(t1,sbfr(sel(4,:),:)); xlabel('sbfr');
subplot(2,2,3); plot(t1,sbf10(sel(4,:),:)); xlabel('sbf10');
subplot(2,2,4); plot(t1,s3); xlabel('simulated source 3');

%% depth bias

pbf = var(sbf,[],2);
pbf = sum(reshape(pbf,3,[]));

source.pos = sourcemodel.pos;
source.dim = sourcemodel.dim;
source.inside = sourcemodel.inside;
source.avg.pow = zeros(size(source.pos,1),1);
source.avg.pow(source.inside)=pbf;

cfg = [];
cfg.funparameter='avg.pow';
cfg.method='slice';
cfg.nslices = 10;
cfg.funcolorlim=[0 0.2];
ft_sourceplot(cfg,source);

%% contrast between 2 conditions
% create the sensor data for the second condition
sensordata2 = 1.25.*leadfield1*s1+0.8.*leadfield2*s2+0.8.*leadfield3*s3+1.25.*leadfield4*s4+randn(301,1000)*0.04e-8;

% compute the covariance
C2  = cov([sensordata sensordata2]');
iC2 = pinv(C2); % so that we compute it only once
iCr2 = inv(C2+eye(301)*1e-19);

% compute the beamformer spatial filter
for ii = 1:size(L,2)/3
  indx=(ii-1)*3+(1:3);
  Lr = L(:,indx);  % Lr is the leadfield for source r
  wbfr2(indx,:)=pinv(Lr'*iCr2*Lr)*Lr'*iCr2;
end
sbfr2 = wbfr2*sensordata2;
pbfr2 = var(sbfr2,[],2);
pbfr2 = sum(reshape(pbfr2,3,[]));
sbfr1 = wbfr2*sensordata;
pbfr1 = var(sbfr1,[],2);
pbfr1 = sum(reshape(pbfr1,3,[]));
source.avg.pow(source.inside)=(pbfr1-pbfr2)./(pbfr1+pbfr2);

cfg = [];
cfg.funparameter='avg.pow';
cfg.method='slice';
cfg.nslices = 10;
cfg.funcolorlim=[-.2 .2];
ft_sourceplot(cfg,source);

%% correlated sources
sens = ni2_sensors('type','meg');
headmodel = ni2_headmodel('type','spherical','nshell',1);
leadfield1 = ni2_leadfield(sens,headmodel,[-5.3 0 5.9 1 0 0]); % position 2342 in grid
leadfield2 = ni2_leadfield(sens,headmodel,[4.9 0 6.2 0 1 0]); % position 2352 in grid

% create the time course of activation
[s1,t1]=ni2_activation('latency',.5,'frequency',10);
[s2,t2]=ni2_activation('latency',.478,'frequency',10);
corr(s1',s2')

% create the sensor data
sensordata = leadfield1*s1+leadfield2*s2+randn(301,1000)*0.04e-8;

% create a 3D grid
sourcemodel = ni2_sourcemodel('type','grid','resolution',1);

% use FieldTrip to quickly compute the forward model
cfg = [];
cfg.grid = sourcemodel;
cfg.grad = sens;
cfg.vol  = headmodel;
% cfg.normalize = 'yes';
sourcemodel = ft_prepare_leadfield(cfg);
L = cat(2,sourcemodel.leadfield{sourcemodel.inside});

% compute the covariance
C  = cov(sensordata');
iCr = inv(C+eye(301)*5e-19);

% compute the beamformer spatial filter
for ii = 1:size(L,2)/3
  indx=(ii-1)*3+(1:3);
  Lr = L(:,indx);  % Lr is the leadfield for source r
  wbfr(indx,:)=pinv(Lr'*iCr*Lr)*Lr'*iCr;
end
sbfr = wbfr*sensordata;

sel = find(ismember(find(sourcemodel.inside), [2342 2352]));
sel = repmat((sel-1)*3,1,3)+repmat(1:3,numel(sel),1);

figure;
subplot(1,2,1); plot(t1,sbfr(sel(1,:),:));
subplot(1,2,2); plot(t1,s1);

figure;
subplot(1,2,1); plot(t1,sbfr(sel(2,:),:));
subplot(1,2,2); plot(t1,s2);

pbf = var(sbfr,[],2);
pbf = sum(reshape(pbf,3,numel(pbf)/3));

source.pos = sourcemodel.pos;
source.dim = sourcemodel.dim;
source.inside = sourcemodel.inside;
source.avg.pow = zeros(size(source.pos,1),1);
source.avg.pow(source.inside)=pbf;
 
cfg=[];
cfg.funparameter='avg.pow';
cfg.method='slice';
cfg.nslices = 4;
cfg.slicerange=[5 8];
cfg.funcolorlim=[0 0.08];
ft_sourceplot(cfg,source);
