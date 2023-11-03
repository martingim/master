function [pc,mm] = phasecorr2(f1,m1,f2,m2,~,~,alpha)
% PHASECORR computes the phase correlation between
% images f1 and f2. f1 and f2 is multiplied with the masks m1 and m2
% before calculating the phase correlation, but is not a fully masked
% generalization of phase correlation.
%
% [pc,mm] = PHASECORR(f1,m1,f2,m2,[],[],alpha)
% 
% It is possible to adjust exponent alpha in the magnitude normalization.
% The default value is alpha = 1 (normal phase correlation).
% 
% % Example
% 
% im1 = im2double(imread('imA.png'));
% idx = (1:32);
% A = im1(idx,idx);
% B = im1(idx+4,idx+5);
%
% % 
% pc = PHASECORR(A,[],B,[]);
%
% % 
% xm = 8;
% idy = (-xm:xm) + 32; 
% x0 = subpixelnone(cc(idy,idy));
% dx = x0-xm-1
%
% See also maskedcc



  [ma,na] = size(f1);
  [mb,nb] = size(f2);
    
  if (nargin < 7)
    alpha = 1;
  end
    
  % apply mask
  if(isempty(m1))
    m1 = ones(ma,na);   
  end
  if(isempty(m2))
    m2 = ones(mb,nb);
  end
  
  % reverse conjugate 
  % (change to f1 to get same as normxcorr2
  f2 = conj(f2(mb:-1:1,nb:-1:1));  
    
  % Pad images
  m  = 2^nextpow2(ma+mb);
  n  = 2^nextpow2(na+nb);  
  
  [ma,na] = size(f1);
  [mb,nb] = size(f2);
    
  FF = fft2(f1.*m1,m,n).*fft2(f2.*m2,m,n);
  MM = fft2(m1,m,n).*fft2(m2,m,n);
  %FF = conv2(FF,ones(3,3)/9,'same');
  R = FF./abs(FF).^(alpha);     
  mm = ifft2(MM);
  
  %S = MM./abs(MM).^(alpha); % allways fails if alpha~=0, the mask does not have phase information
  pc = ifft2(R)./mm.^(1-alpha);      
  %assignin('base','R',R)  
  %assignin('base','S',S)
      
  % remove areas with insufficient overlap
  %ncc(mm < 0.1*max(ma,mb)*max(na,nb)) = 0;  c
   
  pc = real(pc);
  
  % trim to standard size    
  pc(ma+mb:m,:) = [];
  pc(:,na+nb:n) = [];
  
  if(nargout>1)     
    mm = real(mm)/min(ma,mb)/min(na,nb);       
    mm(ma+mb:m,:) = [];
    mm(:,na+nb:n) = [];    
  end
  
 

  
  

  