function sens = ni2_sensors(varargin)

% NI2_SENSORS create a fictive sensor-array.
%
% Use as
%  sens = ni2_sensors('type', senstype);
%
% Where senstype can be any of the following:
%  'eeg' generates a 91-channel eeg sensor array
%  'meg' generates a 301-channel meg magnetometer sensor array
%  'meg_grad_axial' generates a 301-channel meg axial gradiometer sensor
%  array.


type   = ft_getopt(varargin, 'type', 'eeg');
jitter = ft_getopt(varargin, 'jitter', 0);
n      = ft_getopt(varargin, 'n', 162); % number of vertices for the sphere, determines the number of electrodes, this is the old default

switch type
  case 'eeg'
    % create an electrode array
    [chanpos, tri] = mesh_sphere(n);
    chanpos        = chanpos*10;
    chanpos(chanpos(:,3)<0,:) = [];

    if jitter
      [th,phi,r] = cart2sph(chanpos(:,1), chanpos(:,2), chanpos(:,3));
      shift1 = 2*jitter*(rand(numel(th),1)  - 1);
      shift2 = 2*jitter*(rand(numel(phi),1) - 1);
      [chanpos(:,1),chanpos(:,2),chanpos(:,3)] = sph2cart(th+shift1(:),phi+shift2(:),r);
    end

    sens.unit     = 'cm';
    sens.coordsys = 'ctf';
    sens.chanpos  = chanpos;
    sens.elecpos  = chanpos;
    for k = 1:size(chanpos,1)
      sens.label{k,1}    = sprintf('eeg% 02d', k);
      sens.chantype{k,1} = 'eeg';
    end
    sens = ft_datatype_sens(sens);

  case {'meg'}
    % create a magnetometer array
    [chanpos, tri] = icosahedron642;
    coilori        = chanpos;
    chanpos        = chanpos*12.5;
    chanpos(:,3)   = chanpos(:,3) - 1.5;

    z = chanpos(:,3);
    coilori        = coilori(z>0,:);
    chanpos        = chanpos(z>0,:);
    nchan          = size(chanpos,1);

    sens.unit     = 'cm';
    sens.coordsys = 'ctf';
    sens.chanpos  = chanpos;
    sens.chanori  = coilori;
    sens.coilpos  = chanpos;
    sens.coilori  = coilori;
    sens.tra      = eye(nchan);
    for k = 1:nchan
      sens.label{k,1}    = sprintf('meg% 03d', k);
      sens.chantype{k,1} = 'megmag';
    end
    sens = ft_datatype_sens(sens);

  case 'ctf151'
    load('ctf151');

    % only keep the normal MEG channels
    sel = strcmp(sens.chantype, 'meggrad');
    sens.label    = sens.label(sel);
    sens.chantype = sens.chantype(sel);
    sens.chanunit = sens.chanunit(sel);
    sens.tra      = sens.tra(sel,:);
    sens.chanpos  = sens.chanpos(sel,:);
    sens.chanori  = sens.chanori(sel,:);

  case 'ctf275'
    load('ctf275');

    % only keep the normal MEG channels
    sel = strcmp(sens.chantype, 'meggrad');
    sens.label    = sens.label(sel);
    sens.chantype = sens.chantype(sel);
    sens.chanunit = sens.chanunit(sel);
    sens.tra      = sens.tra(sel,:);
    sens.chanpos  = sens.chanpos(sel,:);
    sens.chanori  = sens.chanori(sel,:);

  case 'bti248'
    load('bti248');

    % only keep the normal MEG channels
    sel = strcmp(sens.chantype, 'megmag');
    sens.label    = sens.label(sel);
    sens.chantype = sens.chantype(sel);
    sens.chanunit = sens.chanunit(sel);
    sens.tra      = sens.tra(sel,:);
    sens.chanpos  = sens.chanpos(sel,:);
    sens.chanori  = sens.chanori(sel,:);

  case 'neuromag306'
    load('neuromag306');

  case 'eeg1020'
    load('eeg1020');

  case 'eeg1010'
    load('eeg1010');

  otherwise
    error('unsupported type % s', type);
end
