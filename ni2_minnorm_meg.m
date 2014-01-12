% Minimum norm estimate for source localization

%% non-uniqueness in underdetermined system
L=[1 0 1;  0 1 1];

y=[2 1]';

% 2 = x +     z
% 1 =     y + z
% For any z, we can find x and y that solve the system

% Try z = 0
% x=2, y=1;

% Try z = 3;
% x=-1, y = -2

% Try for yourself two more values of z.  Can you always find a set of x
% and y that solves the equations?

% Try z = 2;
% x=0, y = -1

snorm1=sqrt(s(1)^2+s(2)^2+s(3)^2)
snorm2=sqrt(sum(s(:).^2))
snorm3=norm(s)

%three test solutions are:

s=[2 1 0];
norm(s)
s=[-1 -2 3];
norm(s)
s=[0 -1 2];
norm(s)

% pinv

Lpinv=L'*(L*L')^-1
Lpinv=pinv(L)

Lpinv1=L\eye(2)

% Now try pinv
s_mn=Lpinv*[2 1]'
% [1 0 1]
lf*smn

%% Null space of leadfield

% 2 sources, one radial, one tangential (MEG sphere)
% see that if dot of source to LF is null, then can never recover.

L=[1 0 1;  0 1 1];

s1=[2 1 0]';
s2=[-1 -2 3]';

L*s1+L*s2
s3=[1 1 -1]';

L*(s1+s3)


%% tikhonv regularisation and sensor noise

% add sensor noise
n=[.1 -.3]'
y_n=y+n

Lpinv=L'*(L*L')^-1

Lpinv*y_n
norm(Lpinv*y_n)

Rn=diag(diag(n*n'))
% Rn=n*n'

Lpinv_reg=L'*(L*L'+Rn)^-1;
Lpinv_reg*y_n;
norm(Lpinv_reg*y_n)

Lpinv_reg=L'*(L*L'+2*Rn)^-1;
Lpinv_reg*y_n;
norm(Lpinv_reg*y_n)

Lpinv_reg=L'*(L*L'+200*Rn)^-1;
Lpinv_reg*y_n;
norm(Lpinv_reg*y_n)

Lpinv_reg=L'*(L*L'+2e15*Rn)^-1;
Lpinv_reg*y_n;
norm(Lpinv_reg*y_n)


%% effect of column norm (depth weighting)

L=[4 0 1;  0 4 1];
y=[1 1]';

Lpinv=L'*(L*L')^-1;
s_mn=Lpinv*y

% compute colnorm
Lcolnorm=sqrt(sum(L.^2,1))
% Lcolnorm=sum(L.^2,1)

% R_s=diag(Lcolnorm.^2)
R_s=diag(Lcolnorm)

Lpinv_dw=inv(R_s)*L'*(L*inv(R_s)*L')^-1

Lpinv_dw*y

L_normed=L./repmat(Lcolnorm,2,1)

pinv(L_normed)*y

[L./repmat(nut_colnorm(L),2,1)]\eye(2)*y




%% First some matrix inversion
% 
% % Invert square matrix
% a=diag([1 2 3 4])
% inv(a)
% 
% b=[1 2;3 4]
% inv(b)
% 
% 1/[b(1,1)*b(2,2)-b(1,2)*b(2,1)] * [b(2,2) -b(1,2); -b(2,1) b(1,1)]
% 
% % Can't invert non-square
% c=[1 2 3; 4 5 6]
% inv(c)
% 
% 
% % Verify that pinv is l'*inv(l*l')
% ci=pinv(c);
% ci*c
% c*ci
% 
% d=[1 2; 3 4; 5 6]
% di=pinv(d)
% 
% % is pinv still equal to l'*inv(l*l')
% 
% % what happens if lf not independent rows?
% e=[1 10;3 30]
% e=[1 10;3 30.0001]
% 
% e(1,:)=[1 4 5];
% e(2,:)=[2 -4 -3];
% e(3,:)=2*e(1,:) + 5*e(2,:);




%% apply to simulations



% 2 sources plus heartbeat from outside 
% (component of heartbeat with LF will project in)
% assumption that all activity must come from within Lc


% resolution matrix
% s_est = w * L * s
% source only at s2, what is value at s1?

%% minimization algos

% % show how source prior weights importance of that error minimization more
% % than for other sources
% 
% strue=[3; 5; 7];
% c=diag([3 1 1]);
% 
% sest=[1; 5; 9]
% a=sest-strue;
% a'*c*a
% 
% sest=[2; 6; 10]
% a=sest-strue;
% a'*c*a
% 
% % weight 2nd source more
% c=diag([1 3 1]);
% 
% sest=[1; 5; 9]
% a=sest-strue;
% a'*c*a
% 
% sest=[2; 6; 10]
% a=sest-strue;
% a'*c*a
% 


%% apply to real data

% What do we need to perform source localization with a minimum norm
% method?
% 1) Data
load ~/teaching/megdatanoisenoise.mat 

% 2) Sensor positions

load('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/dataFIC')
datanoisenoise.grad=dataFIC.grad;
clear dataFIC
save ~/teaching/megdatanoisenoise.mat datanoisenoise

%%
cd ~/teaching/ni2_jmz/
load ~/teaching/megdatanoisenoise.mat 


% 3) Head model
% For EEG, we'll be interested to test the effect of using the more
% accurate BEM versus a less accurate 3sphere and 1 sphere model.
load('/home/common/matlab/fieldtrip/data/ftp/tutorial/headmodel_meg/vol.mat')
vol=ft_convert_units(vol,'cm');

cfg=[];
cfg.method = 'singlesphere';
vol_1sph=ft_prepare_headmodel(cfg,vol.bnd);

% 4) Grid points in the head
load meg_leadfields.mat

grid=rmfield(grid2,'leadfield');

cfg = [];
cfg.grid = grid;
cfg.grad = datanoisenoise.grad;
cfg.vol = vol_1sph;
grid_1sph = ft_prepare_leadfield(cfg, datanoisenoise);

% Convert data from 'raw' to 'timelock' and to 'freq';
cfg=[];
cfg.keeptrials='yes';
cfg.covariance='yes';
tlock_tr=ft_timelockanalysis(cfg,datanoisenoise);

cfg=[];
cfg.covariance='yes';
cfg.covariancewindow = [0 0.3];
tlock_avg=ft_timelockanalysis(cfg,datanoisenoise);

% Now we can compute a source inversion, with our three different
% leadfields based on the three different head model options
cfg = [];
cfg.method = 'mne';
cfg.grid = grid2;
cfg.vol = vol;
cfg.snr = 1;
mne_grid2_avg = ft_sourceanalysis(cfg, tlock_avg);

for ll=1:size(mne_grid2_avg.pos,1)
  mne_grid2_avg.avg.normpow(ll,:) = mne_grid2_avg.avg.pow(ll,:)/trace(mne_grid2_avg.avg.noisecov{ll});
end

cfg=[];
cfg.funparameter = 'avg.normpow';
cfg.interactive = 'yes';
ft_sourceplot(cfg,mne_grid2_avg);

mri=ft_read_mri('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/Subject01.mri');
mne_grid2_avg.anatomy = mri.anatomy;
anat


cfg = [];
cfg.method = 'mne';
cfg.grid = grid_3sph;
cfg.vol = vol_3sph;
cfg.elec = elec;
cfg.lambda = '10%';
source_3sph = ft_sourceanalysis(cfg, tlock)

cfg = [];
cfg.method = 'mne';
cfg.grid = grid_1sph;
cfg.vol = vol_1sph;
cfg.elec = elec;
cfg.lambda = '10%';
source_1sph = ft_sourceanalysis(cfg, tlock)


%%
cfg=[];
cfg.keeptrials='yes';
cfg.output = 'powandcsd';
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.toi = [0:.1:2];
cfg.foi = 5:2.5:70;
cfg.t_ftimwin = .4*ones(size(cfg.foi));;
freq=ft_freqanalysis(cfg,datanoisenoise);

