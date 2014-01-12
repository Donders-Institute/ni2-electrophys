% Blind Source Separation of simulated data

%% Linear mixing of 5 orignal sources, ignoring projection
load ~/teaching/sources.mat

% The six sources are arranged in a matrix 'sources' of size
% [trials X sources X time].  We first reshape so that time is concatenated
% over trials.
sources=reshape(permute(sources,[2 3 1]),[6 100000]);
% Recall that the order of the sources is:
% 1) GammaERS
% 2) AlphaERS
% 3) AlphaERD
% 4) Heartbeat1
% 5) Heartbeat2
% 6) Pink Noise

ni2_subplot(sources(:,1:1000));

% What is the distribution type of each source?
ni2_subplot(sources(:,1:1000),2);

% Now let's create a random linear mixture of these sources.
state=randomseed(5);
mixing=randn(6);
sensors=mixing*sources;

% First try of just 2 sources.
mixing2=mixing([2 4],[2 4]);
sources2=sources([2 4],:);
sensors2=mixing2*sources2;

ni2_subplot(sensors(:,1:1000));
ni2_subplot(sensors2(:,1:1000));

% What is the distribution type of each sensor?
ni2_subplot(sensors(:,:),2);
ni2_subplot(sensors2(:,:),2);

% If we knew ground truth mixing matrix, can we reconstruct sources?
% Recall Matlab exercises from Eric Maris
estsource=mixing\sensors;

figure;plot([sources(1,:)-estsource(1,:)])
ni2_subplot(estsource(:,1:1000));

%% PCA

try
  [U,D,V]=svd(sensors);
end

[U,D,V]=svd(sensors(:,1:10000));
[U2,D2,V2]=svd(sensors2(:,1:10000));

ni2_subplot(V(1:1000,1:6)')
ni2_subplot(V(1:1000,1:6)',2)

ni2_subplot(V2(1:1000,1:2)')
ni2_subplot(V2(1:1000,1:2)',2)

% Are these components actually orthogonal?
figure;imagesc(corr(V(:,1:6)));caxis([-1 1])

% Can we solve it ourselves 'by hand'?

x(1,:)=squeeze(tlock.trial(1,1,:));
x(2,:)=squeeze(tlock.trial(2,1,:));

% Find wT such that wT*R*w is minimized
global x % This allows the variable 'x' to be seen inside the function 'project_cov'
Aunmix=reshape(fminsearch(@project_cov,randn(1,4))',[2 2]);
y=[Aunmix*x];
% figure;subplot(1,2,1);plot(x(1,:),x(2,:),'o')
% subplot(1,2,2);plot(y(1,:),y(2,:),'o')
% Dot product to assess orthogonal
% Aunmix(1,:)*Aunmix(2,:)'
% x(1,:)*x(2,:)'
y(1,:)*y(2,:)'

% Gram-Schmidt process: 
% http://en.wikipedia.org/wiki/Gram-Schmidt_process#The_Gram.E2.80.93Schmidt_process
y(1,:) = x(1,:);
e(1,:) = y(1,:)/norm(y(1,:));
y(2,:) = x(2,:) - ( y(1,:)* x(2,:)')/ ( y(1,:)* y(1,:)')*y(1,:);
e(2,:) = y(2,:)/norm(y(2,:));
% figure;subplot(1,2,1);plot(x(1,:),x(2,:),'o')
% subplot(1,2,2);plot(y(1,:),y(2,:),'o')
x(1,:)*x(2,:)'
y(1,:)*y(2,:)'
e(1,:)*e(2,:)'

    
%% ICA

ft_hastoolbox('fastica', 1);       % see http://www.cis.hut.fi/projects/ica/fastica
% ft_hastoolbox('eeglab', 1);

[fastica_mixing2, fastica_unmixing2] = fastica(sensors2);
estsources2_fastica=fastica_unmixing2*sensors2;
ni2_subplot(estsources2_fastica(:,1:1000));

[fastica_mixing, fastica_unmixing] = fastica(sensors);
estsources_fastica=fastica_unmixing*sensors;
ni2_subplot(estsources_fastica(:,1:1000));

% Has ICA recovered the original time series?
corr(sources2',estsources2_fastica')
corr(sources',estsources_fastica')
% How does this compare to how similar the original sources were to each
% other?
corr(sources',sources')
corr(sources2',sources2')

% Are mixing2 (the original mixing of sources2) and fastica_mixing2 (the estimated 
% mixing from fastica) the same?
% Try normalizing the values:
mixing2/norm(mixing2)
fastica_mixing2/norm(fastica_mixing2)
% If not, how does this work to recover the original time series?
% What does this tell you about how unique the mixing matrix is?

%% Does ICA care about order of time series?  Or only histogram?
[sources_sorted, order_sorted]=sort(sources,2);
ni2_subplot(sources_sorted(:,1:1000))

sensors_sorted = mixing*sources_sorted;
ni2_subplot(sensors_sorted(:,1:1000))

[fastica_mixing_sorted, fastica_unmixing_sorted] = fastica(sensors_sorted);
estsources_sorted=fastica_unmixing_sorted*sensors_sorted;
ni2_subplot(estsources_sorted(:,1:1000));
% Why are there only 5 estsources_sorted?

% Let's re-do with 5 sources to make it easier.
sources5=sources([1 2 3 4 6],:);
[sources5_sorted, order5_sorted]=sort(sources5,2);
sensors5_sorted = mixing([1 2 3 4 6],[1 2 3 4 6])*sources5_sorted;
[fastica_mixing5_sorted, fastica_unmixing5_sorted] = fastica(sensors5_sorted);
estsources5_sorted=fastica_unmixing5_sorted*sensors5_sorted;

% Can we reconstruct the original time series, by knowing order_sorted?
% To do so, we need to know which component of estsources_sorted matches
% with sources_sorted.
corr(estsources5_sorted',sources5_sorted')
[maxcorr,maxind]=max(abs(corr(estsources5_sorted',sources5_sorted')),[],2)
% From the correlation and taking maximums, do you think we can figure out which goes
% with which?
estsources_unsorted = estsources5_sorted(order5_sorted(maxind,:));
ni2_subplot(estsources_unsorted(:,1:1000));

% Nevertheless, are unmixing the same?

% Does type of ICA matter?
ft_hastoolbox('eeglab', 1);
[runica_mixing5_sorted, runica_unmixing5_sorted] = runica_wrapper(sensors5_sorted);
estsources5_sorted=runica_unmixing5_sorted*sensors5_sorted;
ni2_subplot(estsources_sorted(:,1:1000));

corr(estsources5_sorted',sources5_sorted')
[maxcorr,maxind]=max(abs(corr(estsources5_sorted',sources5_sorted')),[],2)
estsources_unsorted = estsources5_sorted(order5_sorted(maxind,:));
ni2_subplot(estsources_unsorted(:,1:1000));

%% What if data were perfect Gaussian-distributed to begin with?

%% sensor data
load ~/teaching/tlock.mat

cfg=[];
cfg.viewmode = 'vertical';
cfg.layout = 'elec1005';
ft_databrowser(cfg,tlock);

% What distribution does the sensor data have?
% First trial, 6 random sensors
ni2_subplot(squeeze(tlock.trial(1,[5:5:30],:)),2);

cfg=[];
cfg.method = 'fastica';
cfg.numcomponent = 20;
cfg.randomseed = 17;
comp_all = ft_componentanalysis(cfg, tlock);

% What do the components look like?
cfg=[];
cfg.viewmode = 'component';
cfg.layout = 'elec1005';
ft_databrowser(cfg,comp_all);
% Can you identify the 6 known underlying components?

% Does it help to use a greater PCA reduction?
cfg=[];
cfg.method = 'fastica';
cfg.numcomponent = 6;
cfg.randomseed = 17;
comp6 = ft_componentanalysis(cfg, tlock)

cfg=[];
cfg.viewmode = 'component';
cfg.layout = 'elec1005';
ft_databrowser(cfg,comp6);

% Which component(s) would you reject at this point?
cfg=[];
cfg.component = [ 1]; % enter here which component to reject
[data_cleanica] = ft_rejectcomponent(cfg, comp_all, tlock)

cfg=[];
cfg.viewmode = 'vertical';
cfg.layout = 'elec1005';
ft_databrowser(cfg,data_cleanica);



% Plot the difference of tlock and tlock_cleanica
cfg=[];
cfg.keeptrials='yes';
tlock_cleanica=ft_timelockanalysis(cfg,data_cleanica);
tlockdiff=tlock;
tlockdiff.trial = tlock_cleanica.trial-tlock.trial;
cfg=[];
cfg.viewmode = 'vertical';
cfg.layout = 'elec1005';
ft_databrowser(cfg,tlockdiff);

%% ICA followed by source localization
% How does removing heartbeat with ICA from sensor data affect the source
% localization results

