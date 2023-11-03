function [im1i,mask1i,im2i,mask2i] = iminterp2lanczosj(im1,mask1,x1,y1,im2,mask2,x2,y2,iparam)
% IMINTERP2LANCZOSJ interpolates a PIV image pair with mask, for 
% use in distortedpass(...) to undistort images. The images are
% interpolated using lanczos resampling.
%
% [im1i,mask1i,im2i,mask2i] = IMINTERP2LANCZOSJ(im1,mask1,x1,y1,im2,mask2,x2,y2)
%
% input  - im1     - First image in PIV image pair.
%        - mask1   - Mask for first image in PIV image pair.
%        - (x1,y1) - Distortion/evaluation points of image 1.
%        - im2     - Second image in PIV image pair.
%        - mask2   - Mask for second image in PIV image pair
%        - (x2,y2) - Distortion/evaluation points of image 2.
%        - iparam  - Filter size, default value is 3.
% 
% output - im1i    - Undistorted interpolation of first image
%        - mask1i  - Undistorted interpolation of first mask
%        - im2i    - Undistorted interpolation of second image
%        - mask2i  - Undistorted interpolation of second mask
% 
% Example
%
% % Read in images
% im1 = im2double(imread('imA.png'));
% im2 = im2double(imread('imB.png'));
%
% % Normal pass
% opt1 = setpivopt();
% piv1 = normalpass([],im1,[],im2,[],opt1);
% 
% % Distorted pass
% opt2 = setpivopt('iminterp',@iminterp2lanczosj,3);
% piv2 = distortedpass(piv1,im1,[],im2,[],opt2);
%
% See also distortedpass, rungekutta, iminterp2bsplinej, iminterp2matlab 

  K = size(im1,3); 
  [M,N] = size(x1);
  
  if(nargin<9)
      iparam = 3;
  end
  
  % Interpolate image pair
  im1i = zeros(M,N,K);  
  im2i = zeros(M,N,K);
  for k=1:K
    im1l = Lanczos2(im1(:,:,k));       
    im2l = Lanczos2(im2(:,:,k));
    im1i(:,:,k) = im1l.evaluate(x1,y1,iparam);     
    im2i(:,:,k) = im2l.evaluate(x2,y2,iparam);
  end
  
  % Interpolate masks
  if(isempty(mask1))
    mask1i = im1l.mask();
  else
    mask1i = interp2(mask1,x1,y1,'nearest',0);
  end
  
  if(isempty(mask2))
    mask2i = im2l.mask();
  else    
    mask2i = interp2(mask2,x2,y2,'nearest',0);  
  end  
  
end