function [cc,mm] = ccorr(f1,m1,f2,m2,~,~,pad)
 % function cc = ccorr(f1,[],f2,[],[],[],pad) 
  [ma,na] = size(f1);
  [mb,nb] = size(f2);
  mn = min(ma,mb)*min(na,nb);
  
  if (nargin < 7)
    pad = 1;
  end
  
  
  % apply mask
  if(isempty(m1))
    m1 = ones(ma,na);
  else
    f1 = f1.*m1;        
  end
  if(isempty(m2))
    m2 = ones(mb,nb);
  else
    f2 = f2.*m2;
  end
  
  % reverse conjugate 
  % (change to f1 to get same as normxcorr2
  f2 = conj(f2(mb:-1:1,nb:-1:1));  
  m2 = conj(m2(mb:-1:1,nb:-1:1));
    
  % Pad images
  if(pad)
    m  = 2^nextpow2(ma+mb);
    n  = 2^nextpow2(na+nb);        
  else
    m  = max(ma,mb); % = mb;
    n  = max(na,nb); % = nb;
  end  
                         
  F1  = fft2(f1,m,n);
  F2  = fft2(f2,m,n);  
      
  cc = ifft2(F1.*F2)/mn;
  %if ~any(any(imag(f1))) && ~any(any(imag(f2)))
    cc = real(cc);
  %end
  
   % trim to standard size
  if(pad)
    cc(ma+mb:m,:) = [];
    cc(:,na+nb:n) = [];   
  else
    % only works for even min(ma,mb) & min(na,nb).    
    r = m - min(ma,mb)/2;
    s = n - min(na,nb)/2;
    cc = circshift(cc,r,1);
    cc = circshift(cc,s,2);
    cc = cc(1:end-1,1:end-1);    
      %??
  end
  
  if(nargout>1)     
    M1 = fft2(m1,m,n);
    M2 = fft2(m2,m,n);             
    mm = ifft2(M1.*M2);
  
    mm = mm/mn;
    if(pad)
      mm(ma+mb:m,:) = [];
      mm(:,na+nb:n) = [];
    else
      mm = circshift(mm,r,1);
      mm = circshift(mm,s,2);
      mm = mm(1:end-1,1:end-1);               
    end
  end
