%% Section 2

% page 2
[data,time] = ni2_activation;
sens = ni2_sensors('type','eeg');
headmodel = ni2_headmodel('type', 'spherical', 'nshell', 3);
leadfield = ni2_leadfield(sens, headmodel, [4.9 0 6.2 0 1 0]);
sensordata = leadfield*data+randn(91,1000)*1e-3;

% visualize it
topo_observed = sensordata(:,500);
figure;ni2_topoplot(sens, topo_observed);

% make a model of the data at a 'arbitrary' location
leadfield = ni2_leadfield(sens,headmodel,[5 5 4]);

% fit the moment with a linear equation
dipmom = leadfield\topo_observed;

% create the modelled data
topo_modelled = leadfield*dipmom;
figure;ni2_topoplot(sens, topo_modelled);

% compute the 'error'
sumsq=sum((topo_observed-topo_modelled).^2)./sum(topo_observed.^2);

% make a model of the data at a another location
leadfield = ni2_leadfield(sens,headmodel,[-5 5 4]);

% fit the moment with a linear equation
dipmom = leadfield\topo_observed;

% create the modelled data
topo_modelled = leadfield*dipmom;
figure;ni2_topoplot(sens, topo_modelled);

% compute the 'error'
sumsq2=sum((topo_observed-topo_modelled).^2)./sum(topo_observed.^2);

%% Section 3

% page 4

% create a 3D grid
sourcemodel=ni2_sourcemodel('type','grid','resolution',1);

% take all positions in a plane with z-coordinate = 6;
pos=sourcemodel.pos(sourcemodel.inside,:); 
sel=find(pos(:,3)==6); 
pos=pos(sel,:);

% brute force approach in a for-loop
sumsq=zeros(size(pos,1),1); 
for k=1:size(pos,1)
  leadfield=ni2_leadfield(sens,headmodel,pos(k,:));
  dipmom=leadfield\topo_observed;
  topo_modelled=leadfield*dipmom; 
  sumsq(k)=sum((topo_observed-topo_modelled).^2)./sum(topo_observed.^2);
end

figure;plot(sumsq);

[m,ix]=min(sumsq);
pos(ix,:) 	

% use the pre-computed leadfields
pos=sourcemodel.pos(sourcemodel.inside,:); 
load('leadfields');
sumsq=ones(size(pos,1),1); 
for k=1:size(pos,1)
  ik=(k-1)*3+(1:3); 
  dipmom=leadfield(:,ik)\topo_observed; 
  topo_modelled=leadfield(:,ik)*dipmom; 
  sumsq(k,1)=sum((topo_observed(:)- topo_modelled(:)).^2)./sum(topo_observed(:).^2); 
end

[m,ix]=min(sumsq);
sourcemodel.pos(sourcemodel.inside(ix),:)

%% Section 4

% use a single time point
topo_observed=sensordata(:,500);
pos=sourcemodel.pos(sourcemodel.inside,:);
load('leadfields');
sumsq=ones(size(pos,1),1); 
for k=1:size(pos,1)
  ik=(k-1)*3+(1:3); 
  dipmom=leadfield(:,ik)\topo_observed;
  topo_modelled=leadfield(:,ik)*dipmom;
  sumsq(k,1)=sum((topo_observed(:)-topo_modelled(:)).^2)./sum(topo_observed(:).^2);
end
[m,ix]=min(sumsq);
pos(ix,:) 	

% use more time points
topo_observed=sensordata(:,490:510);
pos=sourcemodel.pos(sourcemodel.inside,:);
load('leadfields');
sumsq=ones(size(pos,1),1); 
for k=1:size(pos,1)
  ik=(k-1)*3+(1:3); 
  dipmom=leadfield(:,ik)\topo_observed;
  topo_modelled=leadfield(:,ik)*dipmom;
  sumsq(k,1)=sum((topo_observed(:)-topo_modelled(:)).^2)./sum(topo_observed(:).^2);
end
[m,ix]=min(sumsq);
pos(ix,:) 	

%% Section 5

data=[]; 
data.avg=sensordata; 
data.time=time; 
data.label=sens.label; 
data.elec=sens; 
data.dimord='chan_time';

cfg=[]; 
cfg.gridsearch='no'; 
cfg.model='regional'; 
cfg.vol=headmodel; 
cfg.latency=[0.49 0.51]; 
cfg.nonlinear='yes'; 
cfg.numdipoles=1;
dip=ft_dipolefitting(cfg,data);

figure;ni2_topoplot(sens,dip.Vmodel(:,11));
figure;ni2_topoplot(sens,dip.Vdata(:,11));

%% Section 6

[data,time]=ni2_activation; 
[data2,time]=ni2_activation('frequency',11,'latency',0.48); 
sens=ni2_sensors('type','eeg'); 
headmodel=ni2_headmodel('type','spherical','nshell',3); 
leadfield=ni2_leadfield(sens,headmodel,[4.9 0 6.2 0 1 0]); 
leadfield2=ni2_leadfield(sens,headmodel,[-5.3 0 5.9 1 0 0]); 
sensordata=leadfield*data+leadfield2*data2+randn(91,1000)*.7e-3;


topo_observed=sensordata(:,490); 
sourcemodel=ni2_sourcemodel('type','grid','resolution',1); 
pos=sourcemodel.pos(sourcemodel.inside,:); 
load('leadfields');
sumsq=ones(size(pos,1),1); 
for k=1:size(pos,1)
  ik=(k-1)*3+(1:3); 
  dipmom=leadfield(:,ik)\topo_observed; 
  topo_modelled=leadfield(:,ik)*dipmom; 
  sumsq(k,1)=sum((topo_observed(:)-topo_modelled(:)).^2)./sum(topo_observed(:).^2); 
end

[m,ix]=min(sumsq);
pos(ix,:)

data=[]; 
data.avg=sensordata;
data.time=time; 
data.label=sens.label; 
data.elec=sens;
data.dimord='chan_time';

cfg=[]; 
cfg.gridsearch='no'; 
cfg.model='regional'; 
cfg.vol=headmodel; 
cfg.latency=[0.49 0.51]; 
cfg.nonlinear='yes'; 
cfg.numdipoles=1;
dip=ft_dipolefitting(cfg,data);

cfg.numdipoles=2;
dip=ft_dipolefitting(cfg,data);

cfg.dip.pos = [4 0 6;-4 0 6];
dip=ft_dipolefitting(cfg,data);


