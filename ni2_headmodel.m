function headmodel = ni2_headmodel(varargin)

type = ft_getopt(varargin, 'type',    'spherical');
n    = ft_getopt(varargin, 'nshells', 1);

switch type
  case 'spherical'
    switch n
      case 1
        headmodel.o    = [0 0 0];
        headmodel.r    = 10;
        headmodel.c    = 1;
        headmodel.unit = 'cm';
        headmodel      = ft_datatype_headmodel(headmodel);
      case 3
        headmodel.o    = [0 0 0];
        headmodel.r    = [0.88 0.92 1.00]*10;
        headmodel.c    = [1 1/80 1];
        headmodel.unit = 'cm';
        headmodel      = ft_datatype_headmodel(headmodel);
      otherwise
        error('number of spheres other than 1 or 3 is not supported');
    end
  otherwise
    error('unknown type requested');
end
