function ni2_subplot(dat, flag)
% dat matrix where first dimension indicates how many subplots
% flag = 1 for plot, flag = 2 for hist

if ~exist('flag','var')
  flag=1;
end

factors=factor(size(dat,1));

figure;

if flag==1
  for ii=1:size(dat,1)
    subplot(factors(1),prod(factors(2:end)),ii);
    plot(dat(ii,:))
  end
  
elseif flag==2
  
  for ii=1:size(dat,1)
    subplot(factors(1),prod(factors(2:end)),ii);
    hist(dat(ii,:),30)
  end
  
  
end

