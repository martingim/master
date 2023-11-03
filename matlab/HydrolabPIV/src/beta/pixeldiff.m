function [pd1,pd2] = pixeldiff(piv,im1,mask1,im2,mask2)
 
  [MM,NN] = size(im1);  
  opt = piv.opt{end};
  W = opt.Width;
  H = opt.Height;  
  
  if(isempty(mask1))
    mask1 = ones(MM,NN);
  end
  if(isempty(mask2))
    mask2 = ones(MM,NN);
  end  
  
  % Ensure double  
  im1 = im2double(im1);  
  im2 = im2double(im2);
  mask1 = im2double(mask1);
  mask2 = im2double(mask2);
  
  [x0,y0] = meshgrid(1-W:NN+W,1-H:MM+H);     
  [xa,ya,xb,yb] = rungekutta(opt.tableau,piv,x0,y0);     
  
  im1 = interp2(im1,xa,ya,'linear',0);
  im2 = interp2(im2,xb,yb,'linear',0);
  mask1 = interp2(mask1,xa,ya,'nearest',0);
  mask2 = interp2(mask2,xb,yb,'nearest',0);
  
  diff = (im1-im2).*mask1.*mask2;
  
  [xx,yy] = meshgrid(1:NN,1:MM);
  F = scatteredInterpolant(xa(:),ya(:),diff(:));
  pd1 = reshape(F(xx,yy),MM,NN);
  F = scatteredInterpolant(xb(:),yb(:),diff(:));  
  pd2 = reshape(F(xx,yy),MM,NN);
  
  