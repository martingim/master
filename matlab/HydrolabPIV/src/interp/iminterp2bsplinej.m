function [im1,mask1,im2,mask2] = iminterp2bsplinej(im1,mask1,x1,y1,im2,mask2,x2,y2,~)
% IMINTERP2BSPLINEJ interpolates a PIV image pair with mask, for 
% use in distortedpass(...) to undistort images. The images are
% interpolated using cubic bspline with not-a-knot boundary condiction.
%
% [im1i,mask1i,im2i,mask2i] = IMINTERP2BSPLINEJ(im1,mask1,x1,y1,im2,mask2,x2,y2)
%
% input  - im1     - First image in PIV image pair.
%        - mask1   - Mask for first image in PIV image pair.
%        - (x1,y1) - Distortion/evaluation points of image 1.
%        - im2     - Second image in PIV image pair.
%        - mask2   - Mask for second image in PIV image pair
%        - (x2,y2) - Distortion/evaluation points of image 2.
% 
% output - im1i    - Undistorted interpolation of first image
%        - mask1i  - Undistorted interpolation of first mask
%        - im2i    - Undistorted interpolation of second image
%        - mask2i  - Undistorted interpolation of second mask
% 
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
% opt2 = setpivopt('iminterp',@iminterp2bsplinej,[]);
% piv2 = distortedpass(piv1,im1,[],im2,[],opt2);
%
% See also distortedpass, rungekutta, iminterp2lanczosj, iminterp2matlab 

  [I,J,K] = size(im1);
  
  M = I - 3;
  N = J - 3;
  
  % Bspline nodes
  y = 1:I;
  ty = zeros(M+7,1);    
  ty(1:4)     = 1;
  ty(5:M+3)   = 3:I-2;    
  ty(M+4:end) = I;
  %ty = y(ty);
       
  x = 1:J;
  tx = zeros(N+7,1);
  tx(1:4)     = 1;
  tx(5:N+3)   = 3:J-2;    
  tx(N+4:end) = J;
   % tx = x(tx);
  
  % Assemble bspline basis matrix
  A = zeros(I);
  for i=1:I        
    A(:,i) = Bspline2.basis(i,4,y,ty); %  
  end
  [LA,UA] = lu(sparse(A)); % use sparse??
  
  B = zeros(J); 
  for j=1:J         
    B(:,j) = Bspline2.basis(j,4,x,tx);  
  end
  [LB,UB] = lu(sparse(B));  
    
  % Fit bsplines to data
  c1 = zeros(I,J,K);
  c2 = zeros(I,J,K);
      
  for k=1:K            
    for j=1:J
        c1(:,j,k) = UA\(LA\im1(:,j,k));     	
        c2(:,j,k) = UA\(LA\im2(:,j,k));     
    end
    for i=1:I
        c1(i,:,k) =  UB\(LB\c1(i,:,k)');        
        c2(i,:,k) =  UB\(LB\c2(i,:,k)');
    end
  end
  
  % Evaluate undistorted images
  [I,J] = size(x1);
  im1 = zeros(I,J,K);
  im2 = zeros(I,J,K);  
  for k=1:K
    im1bs = Bspline2(c1(:,:,k),tx,ty,4);  
    im2bs = Bspline2(c2(:,:,k),tx,ty,4);    
    im1(:,:,k) = im1bs.evaluate(x1,y1);
    im2(:,:,k) = im2bs.evaluate(x2,y2);
  end
  
  % Evaluate undistorted masks
  if(isempty(mask1))
    mask1 = im1bs.mask();
  else
    mask1 = interp2(mask1,x1,y1,'nearest',0);      
  end
  
  if(isempty(mask2))
    mask2 = im2bs.mask();
  else
    mask2 = interp2(mask2,x2,y2,'nearest',0);
  end
  
end