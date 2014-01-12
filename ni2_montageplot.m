function ni2_montageplot(dat)

siz = size(dat);
n   = siz(3);
n1  = floor(sqrt(n));
n2  = ceil(sqrt(n));

datnew = zeros(siz(2)*n1, siz(1)*n2)+nan;
for k = 1:n
  ix2 = rem(k-1, n2) + 1;
  ix1 = floor((k-1)/n2) + 1;
  
  indx2 = (ix2-1)*siz(2) + (1:siz(2));
  indx1 = (ix1-1)*siz(1) + (1:siz(1));
 
  datnew(indx1, indx2) = dat(:,:,k);
end
figure;imagesc(datnew);axis xy
