function [cc,mm] = standardcc(f1,~,f2,~,~,~,pad)
% STANDARDCC performes cross-correlation on images f1 and f2.
% Padding of the images is the default option.
%
% cc = STANDARDCC(f1,[],f2,[],[],[],pad)
%
% Example
% 
% % Test with shifted subwindow
% im1 = im2double(imread('imA.png'));
% idx = (1:32);
% A = im1(idx,idx);
% B = im1(idx+4,idx+5);
%
% % Calculate cross-correlation
% cc = STANDARDCC(A,[],B,[]);
%
% % Find displacement dx
% xm = 8;
% idy = (-xm:xm) + 32; 
% [x0,delta,out] = subpixel3x3(cc(idy,idy));
% dx = x0-delta-xm-1
%
% % plot part of cross-correlation with at least 25% overlap
% figure;
% imagesc(cc)
%
% See Also maskedcc, maskedncc, maskedccj, maskednccj

  % Check input arguments and provide default values
  if(nargin<4)
    error('At least four input arguments needed.')
  end
 
  [ma,na,ka] = size(f1);
  [mb,nb,kb] = size(f2);
  
  if(ka~=kb)
    error('The number of subwindows in f1 and f2 needs to be the same.');
  end
      
  if (nargin < 7)
    pad = true;
  end   
  
  % Reverse conjugate f2/m2
  % (change to f1/m1 to get simliar as normxcorr2
  f2 = conj(f2(mb:-1:1,nb:-1:1,:));  
    
  % Set padding options
  if(pad)
    m  = 2^nextpow2(ma+mb);
    n  = 2^nextpow2(na+nb);        
  else
    m  = max(ma,mb); % = mb;
    n  = max(na,nb); % = nb;
  end  
  
  % Calculate masked cross-correlation using fft
  F1 = fft(fft(f1,m,1),n,2);
  F2 = fft(fft(f2,m,1),n,2);          
  cc = ifft(ifft(F1.*F2,m,1),n,2)/min(ma,mb)/min(na,nb);
    
  % Ensure real output for real input 
  if ~any(imag(f1(:))) && ~any(imag(f2(:)))
    cc = real(cc);
  end
  
  if(pad)
    % Trim to standard size
    cc(ma+mb:m,:,:) = [];
    cc(:,na+nb:n,:) = [];   
  else    
    % Shift zero displacement to center
    r = floor(min(ma,mb)/2);
    s = floor(min(na,nb)/2);
    tmp = circshift(cc,m-r,1);
    tmp = circshift(tmp,n-s,2);        
    %cc = tmp(1:end-1,1:end-1);         
    
    % Resize to match same size as when padded
    cc = zeros(ma+mb-1,na+nb-1,ka,'like',tmp);
    cc((1:m-1)+r,(1:n-1)+s,:) = tmp(1:end-1,1:end-1,:);      
  end
  
  if(nargout>1)              
    mm = ones(ma+mb-1,na+nb-1,ka);         
  end