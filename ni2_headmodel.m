function headmodel = ni2_headmodel(varargin)

% NI2_HEADMODEL generates a simple headmodel for forward modeling.
%
% Use as
%   headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1);
%   headmodel = ni2_headmodel('type', 'spherical', 'nshell', 3);
%
% The first call creates a singlesphere volume conduction model that can be
% used in combination with MEG sensor-arrays. The second call creates a
% three concentric spheres model that can be used in combination with EEG
% sensor-arrays.
%
% The headmodel is expressed in cm and in the CTF coordinate system.

type = ft_getopt(varargin, 'type', 'spherical');
n    = ft_getopt(varargin, 'nshell', 1);

switch type
  case 'spherical'
    switch n
      case 1
        headmodel.o    = [0 0 0];
        headmodel.r    = 10;
        headmodel.c    = 1;
        headmodel.unit = 'cm';
        headmodel.type = 'singlesphere';
        headmodel      = ft_datatype_headmodel(headmodel);
        headmodel.coordsys = 'ctf';

      case 3
        headmodel.o    = [0 0 0];
        headmodel.r    = [0.88 0.92 1.00]*10;
        headmodel.c    = [1 1/80 1];
        headmodel.unit = 'cm';
        headmodel.type = 'concentricspheres';
        headmodel      = ft_datatype_headmodel(headmodel);
        headmodel.coordsys = 'ctf';

      otherwise
        error('number of spheres other than 1 or 3 is not supported');
    end

  case 'singleshell'
    % note that this is in CTF coordinates
    load('Subject01_headmodel.mat');
    headmodel.coordsys = 'ctf';

  otherwise
    error('unknown headmodel type');
end
