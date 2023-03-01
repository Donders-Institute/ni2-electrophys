function sourcemodel = ni2_sourcemodel(varargin)

% NI2_SOURCEMODEL generates a sourcemodel for forward and inverse modeling. 
%
% Use as:
%  sourcemodel = ni2_sourcemodel('type', 'grid', 'resolution', res)
%
% Where res is a scalar that specifies the spacing between the dipoles (in
% cm). A regular 3-dimensional grid of dipoles is created.

type       = ft_getopt(varargin, 'type',      'grid');
headmodel  = ft_getopt(varargin, 'headmodel', []);
resolution = ft_getopt(varargin, 'resolution', 0.5); % in cm

switch type
  case 'grid'
    if isempty(headmodel)
      resolution = ft_getopt(varargin, 'resolution', 0.5);
      ax         = 0:resolution:9;
      ax         = [-fliplr(ax(2:end)) ax];
      [x,y,z]    = ndgrid(ax, ax, ax(ax>-1));
      in         = false(size(x));
      in(sqrt(x.^2+y.^2+z.^2)<=9) = true;
      
      pos        = [x(:) y(:) z(:)];
      % inside     = find(in);
      % outside    = find(in==0);
      clear x y z
      
      sourcemodel.pos     = pos;
      % sourcemodel.inside  = inside;
      % sourcemodel.outside = outside;
      sourcemodel.inside  = in(:);
      sourcemodel.dim     = [numel(ax) numel(ax) numel(ax(ax>-1))];

    else
      headmodel = ft_convert_units(headmodel);
      origunits = headmodel.unit;
      headmodel = ft_convert_units(headmodel, 'cm');
      
      cfg = [];
      cfg.headmodel  = headmodel;
      cfg.resolution = resolution;
      cfg.inwardshift = -1.*resolution;
      sourcemodel = ft_prepare_sourcemodel(cfg);
      
      sourcemodel = ft_convert_units(sourcemodel, origunits);
    end
    
  case 'mesh'
    load('Subject01_sourcemodel_15684');

  otherwise
    error('unsupported type of sourcemodel requested');
end