function [nmd,mm] = maskednmdj(f1,m1,f2,m2,idh,idw,alpha)
% MASKEDNMDJ performes masked normalized Minkowski difference on
% images f1 and f2 with corresponding masks m1 and m2 using sums in Java.
% The masks has a value between zero (masked) and one (unmasked).
%
% [nmd,mm] = MASKEDNMDJ(f1,m1,f2,m2,idh,idw,alpha)
%
% In additon to the normalized Minkowksi difference, nmd, the function 
% also returns the masked fraction, mm. 
%
% Increase the computational speed by providing the indices idh and idw  
% to partly calculate the normalized Minkowski difference such that only
% nmd(idh,idw) is calculated. 
% The number of computational threads used can be set with DiffMeasure(N).
% 
% Alpha donetes the exponent in the Minkowski difference 
% and rangees from 0 to Inf with special cases 
% - geometric difference (alpha=0)
% - absolute difference  (alpha=1) 
% - quadratic difference (alpha=2, default value)
% - maximum difference   (alpha=Inf)
%
% 
% Example
% 
% % Test with shifted subwindow
% im1 = im2double(imread('imA.png'));
% idx = (1:32);
% A = im1(idx,idx);
% B = im1(idx+4,idx+5);
%
% % Calculate normalized minimum absolute difference
% xm = 8;
% idy = (-xm:xm) + 32; 
% nmd = MASKEDNMDJ(A,[],B,[],idy,idy,1);
%
% % Find displacement dx
% [x0,delta,out] = subpixel3x3(-nmd(idy,idy));
% dx = x0-delta-xm-1
%
% % Plot normalized minimum absolute difference
% figure;
% imagesc(nmd);
%
% See Also maskedmqd, maskednmqd, maskedmdj
 
% Check input arguments and provide default values
if(nargin<4)
    error('At least four input arguments needed.')
end

[ma,na,ka] = size(f1);
[mb,nb,kb] = size(f2);
  
if(ka~=kb)
    error('The number of subwindows in f1 and f2 needs to be the same.');
end
   
if(isempty(m1))
    m1 = ones(ma,na,'like',f1);   
else
    [i,j,k1] = size(m1);
    if(i~=ma || j~=na || (k1~=1 && k1~=ka))
        error('f1 and m1 needs to be of same size.');
    end
end
  
if(isempty(m2))
    m2 = ones(mb,nb,'like',f2); 
else
    [i,j,k2] = size(m2);
    if(i~=mb || j~=nb || (k2~=1 && k2~=kb))
        error('f2 and m2 needs to be of same size.');
    end
end

if(nargin < 6 || isempty(idh) || isempty(idw))
    idh = 1:(ma+mb-1);
    idw = 1:(na+nb-1);
end

if(nargin < 7)
    alpha = 2;
end

 % Load java class
persistent dm;
if(isempty(dm))
    dm = DiffMeasure();
end

% Calculate normalized Minkowski difference
if(k1==1 && k2==1)
    nmd = zeros(ma+mb-1,na+nb-1,ka);    
    for k=1:ka    
        nmd(:,:,k) = dm.maskednmd(f1(:,:,k),m1,f2(:,:,k),m2,idh,idw,alpha);
    end
    mm = dm.getmm();
else
    nmd = zeros(ma+mb-1,na+nb-1,ka);
    mm = zeros(ma+mb-1,na+nb-1,ka);
    for k=1:ka    
        nmd(:,:,k) = dm.maskednmd(f1(:,:,k),m1(:,:,min(k,k1)),f2(:,:,k),m2(:,:,min(k,k2)),idh,idw,alpha);
        mm(:,:,k) = dm.getmm();
    end    
end
 
 
     
 
  
  
  
 
