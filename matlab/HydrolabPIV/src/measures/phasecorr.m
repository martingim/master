function [pc,mm] = phasecorr(f1,m1,f2,m2,~,~,alpha)
% PHASECORR computes the phase correlation between
% images f1 and f2. f1 and f2 is multiplied with the masks m1 and m2
% before calculating the phase correlation, but is not a fully masked
% generalization of phase correlation.
%
% [pc,mm] = PHASECORR(f1,m1,f2,m2,[],[],alpha)
% 
% It is possible to adjust exponent alpha in the magnitude normalization.
% The default value is alpha = 1 (normal phase correlation), 
% while alpha = 0 gives cross-correlation.
% 
% % Example
% 
% im1 = im2double(imread('imA.png'));
% idx = (1:32);
% A = im1(idx,idx);
% B = im1(idx+4,idx+5);
%
% % Calculate phase correlation
% pc = PHASECORR(A,[],B,[]);
%
% % Find peak
% xm = 8;
% idy = (-xm:xm) + 32; 
% x0 = subpixelnone(pc(idy,idy));
% dx = x0-xm-1
%
% figure;
% imagesc(pc)
%
% See also maskedcc

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
    [i,j,k] = size(m1);
    if(i~=ma || j~=na || (k~=1 && k~=ka))
        error('f1 and m1 needs to be of same size.');
    end
  end
  
  if(isempty(m2))
    m2 = ones(mb,nb,'like',f2); 
  else
    [i,j,k] = size(m2);
    if(i~=mb || j~=nb || (k~=1 && k~=kb))
        error('f2 and m2 needs to be of same size.');
    end
  end
  
  if (nargin < 7)
    alpha = 1;  
  end  
    
  % Apply mask
  if(isempty(m1))
    m1 = ones(ma,na,'like',f1);   
  end
  f1 = bsxfun(@times,f1,m1);
  if(isempty(m2))
    m2 = ones(mb,nb,'like',f2);
  end
  f2 = bsxfun(@times,f2,m2);
  
  % Reverse conjugate f2/m2
  % (change to f1/m1 to get simliar as normxcorr2
  f2 = conj(f2(mb:-1:1,nb:-1:1,:));  
    
  % Pad images
  m  = 2^nextpow2(ma+mb);
  n  = 2^nextpow2(na+nb);    
  
  % Calculate Phase correlation
  F1  = fft(fft(f1,m,1),n,2);
  F2  = fft(fft(f2,m,1),n,2);            
    
  FF = F1.*F2;  
  R = FF./abs(FF).^(alpha);      
  pc = ifft(ifft(R,m,1),n,2);         
  pc = real(pc);
  
  % Trim to standard size    
  pc(ma+mb:m,:,:) = [];
  pc(:,na+nb:n,:) = [];
  
  if(nargout>1)     
    % Calcluate masked area
    M1 = fft(fft(m1,m,1),n,2);
    M2 = fft(fft(m2,m,1),n,2);             
    mm = ifft(ifft(M1.*M2,m,1),n,2);
  
    % Normalize to get masked fraction
    mm = mm/min(ma,mb)/min(na,nb);    
    
    % Trim to standard size 
    mm(ma+mb:m,:,:) = [];
    mm(:,na+nb:n,:) = [];    
  end
  
 

  
  

  