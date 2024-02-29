function sourcemodel = ni2_sourcemodel(varargin)

% NI2_SOURCEMODEL generates a sourcemodel for forward and inverse modeling.
%
% Use as
%   sourcemodel = ni2_sourcemodel('type', 'mesh');
%   sourcemodel = ni2_sourcemodel('type', 'grid', 'resolution', res)
%
% The first call returns a source model that consists of many dipoles that form a
% triangulated mesh that describes the cortical sheet.
%
% The second call creates a regular 3-dimensional grid of dipoles. The parameter res
% is a scalar (in cm) that specifies the resolution, i.e., the spacing between the
% dipoles.

type       = ft_getopt(varargin, 'type', 'grid');
headmodel  = ft_getopt(varargin, 'headmodel', []);
resolution = ft_getopt(varargin, 'resolution', 0.5); % in cm

switch type
  case 'grid'
    if isempty(headmodel)
      ax         = 0:resolution:9;
      ax         = [-fliplr(ax(2:end)) ax];
      [x,y,z]    = ndgrid(ax, ax, ax(ax>-1));
      inside     = false(size(x));
      inside(sqrt(x.^2+y.^2+z.^2)<=9) = true;
      
      pos        = [x(:) y(:) z(:)];
      clear x y z
      
      sourcemodel.pos     = pos;
      sourcemodel.inside  = inside(:);
      sourcemodel.dim     = [numel(ax) numel(ax) numel(ax(ax>-1))];
      sourcemodel.unit    = 'cm';

    else
      headmodel = ft_convert_units(headmodel);
      origunits = headmodel.unit;
      headmodel = ft_convert_units(headmodel, 'cm');
      
      cfg = [];
      cfg.headmodel  = headmodel;
      cfg.resolution = resolution; % in cm
      cfg.inwardshift = -1.*resolution;
      sourcemodel = ft_prepare_sourcemodel(cfg);
      
      sourcemodel = ft_convert_units(sourcemodel, origunits);
    end
    
  case 'mesh'
    load('Subject01_sourcemodel_15684');

  otherwise
    error('unsupported type of sourcemodel requested');
end