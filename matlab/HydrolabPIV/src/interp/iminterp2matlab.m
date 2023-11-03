function [im1i,mask1i,im2i,mask2i] = iminterp2matlab(im1,mask1,x1,y1,im2,mask2,x2,y2,iparam)
% IMINTERP2MATLAB interpolates a PIV image pair with mask, for 
% use in distortedpass(...) to undistort images. The images are
% interpolated using the inbuilt interpolation in matlab.
%
% [im1i,mask1i,im2i,mask2i] = IMINTERP2MATLAB(im1,mask1,x1,y1,im2,mask2,x2,y2,iparam)
%
% input  - im1     - First image in PIV image pair.
%        - mask1   - Mask for first image in PIV image pair.
%        - (x1,y1) - Distortion/evaluation points of image 1.
%        - im2     - Second image in PIV image pair.
%        - mask2   - Mask for second image in PIV image pair
%        - (x2,y2) - Distortion/evaluation points of image 2.
%        - iparam  - Interpolation method used by interp2, default is '*spline'.
% 
% output - im1i    - Undistorted interpolation of first image
%        - mask1i  - Undistorted interpolation of first mask
%        - im2i    - Undistorted interpolation of second image
%        - mask2i  - Undistorted interpolation of second mask
% 
% Example: Distorted pass with linear interpolation of image
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
% opt2 = setpivopt('iminterp',@iminterp2matlab,'*linear');
% piv2 = distortedpass(piv1,im1,[],im2,[],opt2);
%
% See also distortedpass, rungekutta, interp2, iminterp2bsplinej, iminterp2lanczosj 

  [I,J,K] = size(im1);
  [M,N] = size(x1);
  
  if(nargin<9)
      iparam = '*spline';
  end
  
  % Interpolate image pair
  im1i = zeros(M,N,K);  
  im2i = zeros(M,N,K);
  for k=1:K
    im1i(:,:,k) = interp2(im1(:,:,k),x1,y1,iparam,0);
    im2i(:,:,k) = interp2(im2(:,:,k),x2,y2,iparam,0);
  end
  
  % Set default masks
  if(isempty(mask1))
      mask1 = ones(I,J);
  end
  
  if(isempty(mask2))
      mask2 = ones(I,J);
  end
  
  % Interpolate masks
  mask1i = interp2(mask1,x1,y1,'nearest',0);
  mask2i = interp2(mask2,x2,y2,'nearest',0);

end