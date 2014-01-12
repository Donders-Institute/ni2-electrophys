% Minimum norm estimate for source localization

% What do we need to perform source localization with a minimum norm
% method?
% 1) Data
load ~/teaching/tlock.mat

% 2) Sensor positions

elec=ft_read_sens('standard_1005.elc');
elec=ft_convert_units(elec,'cm');

% 3) Head model
% For EEG, we'll be interested to test the effect of using the more
% accurate BEM versus a less accurate 3sphere and 1 sphere model.
load standard_bem;
vol=ft_convert_units(vol,'cm');
vol_bem = vol;

cfg=[];
cfg.method = 'concentricspheres';
vol_3sph=ft_prepare_headmodel(cfg,vol);

cfg=[];
cfg.method = 'singlesphere';
vol_1sph=ft_prepare_headmodel(cfg,vol.bnd(1));

% 4) Grid points in the head
load standard_grid3d8mm;  % 'grid' with 5302 points inside MNI head
% Note, these would need to be coregistered, but here they are already

% 5) Leadfield
cfg = [];
cfg.grid = grid;
cfg.elec = elec;
cfg.vol = vol_bem;
grid_bem = ft_prepare_leadfield(cfg, tlock);
cfg = [];
cfg.grid = grid;
cfg.elec = elec;
cfg.vol = vol_3sph;
grid_3sph = ft_prepare_leadfield(cfg, tlock);
cfg = [];
cfg.grid = grid;
cfg.elec = elec;
cfg.vol = vol_1sph;
grid_1sph = ft_prepare_leadfield(cfg, tlock);

% Now we can compute a source inversion, with our three different
% leadfields based on the three different head model options
cfg = [];
cfg.method = 'mne';
cfg.grid = grid_bem;
cfg.vol = vol_bem;
cfg.elec = elec;
cfg.lambda = '10%';
source_bem = ft_sourceanalysis(cfg, tlock)

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

load standard_mri
source_3sph.anatomy = mri.anatomy;
source_3sph=rmfield(source_3sph,'anatomy');

cfg=[];
cfg.funparameter = 'avg.pow';
% cfg.method = 'slice';
ft_sourceplot(cfg,source_3sph);

