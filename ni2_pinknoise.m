function x = c3_pinknoise(Nx,seed, methodflag, exponent, faxis, fs)
% function x = pinknoise(Nx,seed)

if nargin<3
  methodflag = 1;
  % default to the filtered white noise case
end

if exist('seed','var') && ~isempty(seed)
  ftFuncRandomseed=randomseed(seed);
end

if numel(Nx)>1,
  for k = 1:Nx(1)
    if exist('exponent', 'var')
      if size(exponent,1)==Nx(1)
        tmpexponent = exponent(k,:);
      else
        tmpexponent = exponent;
      end
    end
    
    if nargin==1
      x(k,:) = pinknoise(Nx(2));
    elseif nargin==2
      x(k,:) = pinknoise(Nx(2),seed);
    elseif nargin==3
      x(k,:) = pinknoise(Nx(2),seed, methodflag);
    elseif nargin==4
      x(k,:) = pinknoise(Nx(2),seed,methodflag,tmpexponent);
    elseif nargin==5
      x(k,:) = pinknoise(Nx(2),seed,methodflag,tmpexponent,faxis);
    elseif nargin==6
      x(k,:) = pinknoise(Nx(2),seed,methodflag,tmpexponent,faxis,fs);
    end
  end
  return;
end

switch methodflag,
  case 1
    % IIR type filter, one directional
    B = [0.049922035 -0.095993537 0.050612699 -0.004408786];
    A = [1 -2.494956002   2.017265875  -0.522189400];
    nT60 = round(log(1000)/(1-max(abs(roots(A))))); % T60 est.
    v = randn(1,Nx+nT60); % Gaussian white noise: N(0,1)
    x = filter(B,A,v);    % Apply 1/F roll-off to PSD
    x = x(nT60+1:end);    % Skip transient response
  case 2
    % fft-based
    v = randn(1,Nx);
    f = fft(v);
    m = floor(Nx./2);
    if mod(Nx,2)==0
      faxis(1:(m+1))  = 1:(m+1);
      faxis((m+1):Nx) = fliplr(faxis(2:(m+1)));
    else
      faxis(1:(m+1))  = 1:(m+1);
      faxis((m+2):Nx) = fliplr(faxis(2:(m+1)));
    end
    
    if ~exist('exponent', 'var')
      exponent = 0.5;
    end
    v = ifft(f./(faxis.^exponent));
    v = real(v);
    v = v-mean(v);
    x = v./std(v);
  case 3
    % also fft-based
    if ~exist('faxis', 'var') || ~exist('fs', 'var')
      error('''faxis'' and ''fs'' are required');
    end
    
    m      = floor(Nx./2);
    faxis2 = fs.*(0:(Nx-1))./Nx;
    
    
    if mod(Nx,2)==0
      faxis2(1:(m+1))  = 1:(m+1);
      faxis2((m+1):Nx) = fliplr(faxis2(2:(m+1)));
    else
      faxis2(1:(m+1))  = 1:(m+1);
      faxis2((m+2):Nx) = fliplr(faxis2(2:(m+1)));
    end
    
    % filtered with a frequency response
    v = randn(1,Nx);
    f = fft(v);
    
    if mod(Nx,2)==0
      faxis2(1:(m+1))  = 1:(m+1);
      faxis2((m+1):Nx) = fliplr(faxis2(2:(m+1)));
    else
      faxis2(1:(m+1))  = 1:(m+1);
      faxis2((m+2):Nx) = fliplr(faxis2(2:(m+1)));
    end
    
    % create the axis for interpolation
    faxisx   = fs.*(faxis2-1)./Nx;
    exponent = interpn(faxis,sqrt(exponent),faxisx); % assume exponent to be power, take amplitude
    exponent(~isfinite(exponent))=0;
    
    
    v = ifft(f.*exponent);
    v = real(v);
    v = v-mean(v);
    x = v./std(v);
  
  case 4
    % integrate white noise (one direction...)
    x = randn(1,Nx);
    for i = 2:Nx
      x(i) = x(i) + x(i-1);  
    end
    
  otherwise
end
