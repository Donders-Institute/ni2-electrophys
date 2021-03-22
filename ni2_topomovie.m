function ni2_topomovie(sens, dat, time, cfg)

% ensure the sensor description to be according to what FieldTrip expects
sens = ft_datatype_sens(sens);

if nargin<4
  cfg = [];
end

if ~isfield(cfg, 'layout')
  % create layout
  if strncmp(sens.chantype{1}, 'meg', 3)
    tmpcfg.grad = sens;
  else
    tmpcfg.elec = sens;
  end
  cfg.layout = ft_prepare_layout(tmpcfg);
end
layout = cfg.layout;

figure;
[h,xi,yi,zi,x,y] = ni2_topo(layout, dat(:,1), cfg);

axis([-0.6 0.6 -0.6 0.6]);
axis off;
grid off;

th = text(0,-0.6,sprintf('time=% s',num2str(time(1))));

% determine the lower and upper limits
clim = [min(dat(:)) max(dat(:))];

for k = 1:4:size(dat,2)
  % interpolate the topographic data
  Zi = griddata(x', y', dat(:,k), xi, yi, 'v4'); 
  % replace the values outside the head with not-a-number
  Zi(isnan(zi)) = nan;
  set(h, 'CData', Zi); 
  caxis(clim);
  drawnow;
  % pause(0.05);
  set(th, 'string', sprintf('time=% s',num2str(time(k),'% 1.3f')));
end


function [h,xi,yi,zi,x,y] = ni2_topo(layout, dat, cfg)

chanX = layout.pos(1:end-2,1);
chanY = layout.pos(1:end-2,2);
hpos  = 0;
vpos  = 0;

mask    = layout.mask;
outline = layout.outline;

width         = ft_getopt(cfg, 'width',          []);
height        = ft_getopt(cfg, 'height',         []);
gridscale     = ft_getopt(cfg, 'gridscale',      67); % 67 in original
shading       = ft_getopt(cfg, 'shading',        'flat');
interplim     = ft_getopt(cfg, 'interplim',      'outline');
interpmethod  = ft_getopt(cfg, 'interpmethod',   'v4');
style         = ft_getopt(cfg, 'style',          'surfiso'); % can be 'surf', 'iso', 'isofill', 'surfiso'
tag           = ft_getopt(cfg, 'tag',            '');
isolines      = ft_getopt(cfg, 'isolines');
datmask       = ft_getopt(cfg, 'datmask');
parent        = ft_getopt(cfg, 'parent', []);

% check for nans in the data, they can be still left incase people want to
% mask non channels. 
if any(isnan(dat))
  warning('the data passed to ft_plot_topo contains NaNs, these channels will be removed from the data to prevent interpolation errors, but will remain in the mask');
  flagNaN = true;
else
  flagNaN = false;
end
NaNind = isnan(dat);

% everything is added to the current figure
holdflag = ishold;
if ~holdflag
  hold on
end

% layout units can be arbitrary (e.g. pixels for .mat files)
% so we need to compute the right scaling factor
% create a matrix with all coordinates
% from positions, mask, and outline
allCoords = [chanX chanY];
if ~isempty(mask)
  for k = 1:numel(mask)
    allCoords = [allCoords; mask{k}];
  end
end
if ~isempty(outline)
  for k = 1:numel(outline)
    allCoords = [allCoords; outline{k}];
  end
end

naturalWidth  = (max(allCoords(:,1))-min(allCoords(:,1)));
naturalHeight = (max(allCoords(:,2))-min(allCoords(:,2)));

if isempty(width) && isempty(height)
  xScaling = 1;
  yScaling = 1;
elseif isempty(width) && ~isempty(height)
  % height specified, auto-compute width while maintaining aspect ratio
  yScaling = height/naturalHeight;
  xScaling = yScaling;
elseif ~isempty(width) && isempty(height)
  % width specified, auto-compute height while maintaining aspect ratio
  xScaling = width/naturalWidth;
  yScaling = xScaling;
else
  % both width and height specified
  xScaling = width/naturalWidth;
  yScaling = height/naturalHeight;
end

% keep original channel positions
chanXorg = chanX;
chanYorg = chanY;
chanX = chanX(:) * xScaling + hpos;
chanY = chanY(:) * yScaling + vpos;

if strcmp(interplim, 'electrodes'),
  hlim = [min(chanX) max(chanX)];
  vlim = [min(chanY) max(chanY)];
elseif strcmp(interplim, 'mask') && ~isempty(mask),
  hlim = [inf -inf];
  vlim = [inf -inf];
  for i = 1:length(mask)
    hlim = [min([hlim(1); mask{i}(:,1)*xScaling+hpos]) max([hlim(2); mask{i}(:,1)*xScaling+hpos])];
    vlim = [min([vlim(1); mask{i}(:,2)*yScaling+vpos]) max([vlim(2); mask{i}(:,2)*yScaling+vpos])];
  end
else
  hlim = [min(chanX) max(chanX)];
  vlim = [min(chanY) max(chanY)];
end

% check if all mask point are inside the limits otherwise redefine mask
newpoints = [];
if length(mask)==1
  % which channels are outside
  outside = false(length(chanX),1);
  inside  = inside_contour([chanX chanY], mask{1});
  outside = ~inside;
  newpoints = [chanX(outside) chanY(outside)];
end

if ~isempty(mask)
  % convert the mask into a binary image
  maskimage = false(gridscale);
  % hlim      = [min(chanX) max(chanX)];
  % vlim      = [min(chanY) max(chanY)];
  xi        = linspace(hlim(1), hlim(2), gridscale);   % x-axis for interpolation (row vector)
  yi        = linspace(vlim(1), vlim(2), gridscale);   % y-axis for interpolation (row vector)
  [Xi,Yi]   = meshgrid(xi', yi);
  if ~isempty(newpoints) && (hpos == 0 || vpos == 0)
    warning('Some points fall outside the outline, please consider using another layout')
    % FIXME: I am not sure about it, to be tested!
    %     tmp = [mask{1};newpoints];
    %     indx = convhull(tmp(:,1),tmp(:,2));
    %     mask{1} = tmp(indx,:);
    % NOTE: if you set hpos and/or vpos, newpoints is not empty, but nothing
    % needs to be fixed (this fixme screws up things, then)
  end
  for i = 1:length(mask)
    mask{i}(:,1) = mask{i}(:,1)*xScaling+hpos;
    mask{i}(:,2) = mask{i}(:,2)*yScaling+vpos;
    mask{i}(end+1,:) = mask{i}(1,:);                   % force them to be closed
    maskimage(inside_contour([Xi(:) Yi(:)], mask{i})) = true;
  end
  
else
  maskimage = [];
end

% adjust maskimage to also mask channels as specified in maskdat
if ~isempty(datmask)
  xi           = linspace(hlim(1), hlim(2), gridscale);   % x-axis for interpolation (row vector)
  yi           = linspace(vlim(1), vlim(2), gridscale);   % y-axis for interpolation (row vector)
  maskimagetmp = griddata(chanX', chanY, datmask, xi', yi, 'nearest'); % interpolate the mask data
  if isempty(maskimage)
    maskimage = maskimagetmp;
  else
    maskimagetmp2 = maskimage + maskimagetmp;
    maskimage = maskimagetmp2 > 1.01;
  end
end

% take out NaN channels if interpmethod does not work with NaNs
if flagNaN && strcmp(interpmethod,'v4')
  dat(NaNind) = [];
  chanX(NaNind) = [];
  chanY(NaNind) = [];
end

% interpolate data
xi         = linspace(hlim(1), hlim(2), gridscale);       % x-axis for interpolation (row vector)
yi         = linspace(vlim(1), vlim(2), gridscale);       % y-axis for interpolation (row vector)
[Xi,Yi,Zi] = griddata(chanX', chanY, dat(:,1), xi', yi, interpmethod); % interpolate the topographic data

if ~isempty(maskimage)
  % apply anatomical mask to the data, i.e. that determines that the interpolated data outside the circle is not displayed
  Zi(~maskimage) = NaN;
end

if exist('maskimagetmp')
  maskimagetmp(~maskimage) = NaN;
end

% plot the outline of the head, ears and nose
for i = 1:length(outline)
  xval = outline{i}(:,1) * xScaling  + hpos;
  yval = outline{i}(:,2) * yScaling + vpos;
  ft_plot_vector(xval, yval, 'Color','k', 'LineWidth',2, 'tag', tag, 'parent', parent);
end

% Create isolines
if strcmp(style,'iso') || strcmp(style,'surfiso')
  if ~isempty(isolines)
    [cont, h] = contour(Xi,Yi,Zi,isolines,'k');
    set(h, 'tag', tag);
    if ~isempty(parent)
      set(h, 'Parent', parent);
    end
  end
end

% Plot surface
if strcmp(style,'surf') || strcmp(style,'surfiso')
  deltax = xi(2)-xi(1); % length of grid entry
  deltay = yi(2)-yi(1); % length of grid entry
  h = surf(Xi-deltax/2,Yi-deltay/2,zeros(size(Zi)), Zi, 'EdgeColor', 'none', 'FaceColor', shading);
  set(h, 'tag', tag);
  if ~isempty(parent)
    set(h, 'Parent', parent);
  end
  % if exist('maskimagetmp')
  %  set(h, 'facealpha', 'flat');
  %  set(h, 'alphadatamapping', 'scaled');
  %  set(h, 'alphadata', maskimagetmp);
  % end
end

% Plot filled contours
if strcmp(style,'isofill') && ~isempty(isolines)
  [cont,h] = contourf(Xi,Yi,Zi,isolines,'k'); 
  set(h, 'tag', tag);
  if ~isempty(parent)
    set(h, 'Parent', parent);
  end
end

if ~holdflag
  hold off
end

x = chanX;
y = chanY;
xi = Xi;
yi = Yi;
zi = Zi;
