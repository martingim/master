function [im1i,mask1i,im2i,mask2i] = iminterp2bspline(im1,mask1,xa,ya,im2,mask2,xb,yb)

  [I,J,K] = size(im1);
  [M,N] = size(xa);
  %im1 = padarray(im1,[M-I N-J],0,'post'); % 
  %im2 = padarray(im2,[M-I N-J],0,'post');
  
  im1i = zeros(M,N,K);  
  im2i = zeros(M,N,K);
  for k=1:K
    im1i(:,:,k) = interp2(im1(:,:,k),xa,ya,'spline',0);
    im2i(:,:,k) = interp2(im2(:,:,k),xb,yb,'spline',0);
  end
  
  if(isempty(mask1))
      mask1 = ones(I,J);
  end
  
  if(isempty(mask2))
      mask2 = ones(I,J);
  end
  
  mask1i = interp2(mask1,xa,ya,'nearest',0);
  mask2i = interp2(mask2,xb,yb,'nearest',0);

end