function [cc,mm] = maskedccgpu(f1,m1,f2,m2,~,~,pad)
% MASKEDCC performes masked cross-correlation on
% images f1 and f2 with corresponding masks m1 and m2 using fft.
% The masks has a value between zero (masked) and one (unmasked).
% Padding of the images is the default option.
%
% [cc,mm] = MASKEDCC(f1,m1,f2,m2,[],[],pad)
%
% In additon to the cross-correlation, cc, the function also returns
% the masked fraction, mm. 
% 
% The fifth and sixth argument are for interchangability with maskedccj
% and are ignored.
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
% [cc,mm] = MASKEDCC(A,[],B,[]);
%
% % Find displacement dx
% xm = 8;
% idy = (-xm:xm) + 32; 
% [x0,delta,out] = subpixel3x3(cc(idy,idy));
% dx = x0-delta-xm-1
%
% % plot part of cross-correlation with at least 25% overlap
% figure;
% imagesc(cc.*(mm>=.25))
%
% See Also maskedncc, maskedccj, maskednccj

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
    pad = true;
  end   
  
  f1 = gpuArray(single(f1));
  f2 = gpuArray(single(f2));
  m1 = gpuArray(single(m1));
  m2 = gpuArray(single(m2));
  
  % Reverse conjugate f2/m2
  % (change to f1/m1 to get simliar as normxcorr2
  f2 = conj(f2(mb:-1:1,nb:-1:1,:));  
  m2 = conj(m2(mb:-1:1,nb:-1:1,:));
    
  % Set padding options
  if(pad)
    m  = 2^nextpow2(ma+mb);
    n  = 2^nextpow2(na+nb);        
  else
    m  = max(ma,mb); % = mb;
    n  = max(na,nb); % = nb;
  end  
  
  % Calculate masked cross-correlation using fft
  F1 = fft(fft(bsxfun(@times,f1,m1),m,1),n,2);
  F2 = fft(fft(bsxfun(@times,f2,m2),m,1),n,2);  
  M1 = fft(fft(m1,m,1),n,2);
  M2 = fft(fft(m2,m,1),n,2); 
  %clear f1 f2 m1 m2
  
  mm = ifft(ifft(M1.*M2,m,1),n,2);          
  cc = bsxfun(@rdivide,ifft(ifft(F1.*F2,m,1),n,2),mm);
  %clear F1 F2 M1 M2
  cc = gather(cc);
  mm = gather(mm);
  
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
    % Normalize to get masked fraction
    mm = mm/min(ma,mb)/min(na,nb);
    
    if(pad)
      % Trim to standard size
      mm(ma+mb:m,:,:) = [];
      mm(:,na+nb:n,:) = [];
    else
      % Shift zero displacement to center
      tmp = circshift(mm,m-r,1);
      tmp = circshift(tmp,n-s,2);      
      %mm = tmp(1:end-1,1:end-1);  
      
      % Resize to match same size as when padded
      mm = zeros(ma+mb-1,na+nb-1,ka,'like',tmp);
      mm((1:m-1)+r,(1:n-1)+s,:) = tmp(1:end-1,1:end-1,:);                 
    end
  end