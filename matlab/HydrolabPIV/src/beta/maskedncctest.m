function [ncc,mm] = maskedncctest(f1,m1,f2,m2,~,~,pad)
% MASKEDNCC performes masked normalized cross-correlation on
% images f1 and f2 with corresponding masks m1 and m2 using fft.
% The masks has a value between zero (masked) and one (unmasked).
% Padding of the images is the default option.
%
% [cc,mm] = MASKEDNCC(f1,m1,f2,m2,[],[],pad)
%
% In additon to the normalized cross-correlation, ncc, the function 
% also returns the masked fraction, mm. 
% 
% The fifth and sixth argument are for interchangability with maskednccj
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
% [ncc,mm] = MASKEDNCC(A,[],B,[]);
%
% % Find displacement dx
% xm = 8;
% idy = (-xm:xm) + 32; 
% [x0,delta,out] = subpixel3x3(ncc(idy,idy));
% dx = x0-delta-xm-1
%
% % plot part of cross-correlation with at least 25% overlap
% figure;
% imagesc(ncc.*(mm>=.25))
%
% See Also maskedcc, maskedccj, maskednccj

  % Check input arguments and provide default values
  if(nargin<4)
    error('At least four input arguments needed.')
  end
  
  [ma,na] = size(f1);
  [mb,nb] = size(f2);
  
   
  if(isempty(m1))
    m1 = ones(ma,na,'like',f1);   
  else
    [i,j] = size(m1);
    if(i~=ma || j~=na)
        error('f1 and m1 needs to be of same size.');
    end
  end
  
  if(isempty(m2))
    m2 = ones(mb,nb,'like',f2); 
  else
    [i,j] = size(m2);
    if(i~=mb || j~=nb)
        error('f2 and m2 needs to be of same size.');
    end
  end
  
  if (nargin < 7)
    pad = true;  
  end  
    
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
  
  temp1 = zeros(m,n,6);
  temp1(1:ma,1:na,1) = f1.*m1;
  temp1(1:mb,1:mb,2) = f2.*m2;
  temp1(1:ma,1:na,3) = m1;
  temp1(1:mb,1:nb,4) = m2;
  temp1(1:ma,1:na,5) = f1.*temp1(1:ma,1:na,1);
  temp1(1:mb,1:nb,6) = f2.*temp1(1:mb,1:nb,1); 
  temp1 = fft(fft(temp1));
  
  temp2 = temp1(:,:,4);
  temp1(:,:,6) = temp1(:,:,6).*temp1(:,:,3);  % MFF
  temp1(:,:,5) = temp1(:,:,5).*temp1(:,:,4);  % FFM
  temp1(:,:,4) = temp1(:,:,4).*temp1(:,:,3);  % MM
  temp1(:,:,3) = temp1(:,:,3).*temp1(:,:,2);  % MF
  temp1(:,:,1) = temp1(:,:,1).*temp1(:,:,2);  % FF
  temp1(:,:,2) = temp1(:,:,2).*temp2;         % FM   
  temp1 = ifft(ifft(temp1));
    
  % Calculate masked normalized cross-correlation using fft  
  num  = temp1(:,:,1) - temp1(:,:,2).*temp1(:,:,3)./temp1(:,:,4);
  den1 = temp1(:,:,5) - temp1(:,:,2).*temp1(:,:,2)./temp1(:,:,4);
  den2 = temp1(:,:,6) - temp1(:,:,3).*temp1(:,:,3)./temp1(:,:,4);
  
  ncc = num./sqrt(den1.*den2);  
  mm = temp1(:,:,4);
  
  % Ensure real output for real input
  if ~any(imag(f1(:))) && ~any(imag(f2(:)))
    ncc = real(ncc);
  end
  
  if(pad)
    % Trim to standard size
    ncc(ma+mb:m,:,:) = [];
    ncc(:,na+nb:n,:) = [];   
  else
    % Shift zero displacement to center 
    r = floor(min(ma,mb)/2);
    s = floor(min(na,nb)/2);
    tmp = circshift(ncc,m-r,1);
    tmp = circshift(tmp,n-s,2);        
    %ncc = tmp(1:end-1,1:end-1);  
    
    % Resize to match same size as when padded       
    ncc = zeros(ma+mb-1,na+nb-1,ka,'like',tmp);
    ncc((1:m-1)+r,(1:n-1)+s,:) = tmp(1:end-1,1:end-1,:);        
  end
  
  if(nargout>1) 
      %if ~any(imag(m1(:))) && ~any(imag(m2(:)))
      %  mm = real(mm);
      %end
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
  
  
 
