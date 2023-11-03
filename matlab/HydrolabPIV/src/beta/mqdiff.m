function [mqd,mm] = mqdiff(f1,m1,f2,m2,~,~,pad)
 % function nmqd = maskednmqd(f1,m1,f2,m2,pad) 
  [ma,na] = size(f1);
  [mb,nb] = size(f2);
  
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
  
  if (nargin < 7)
    pad = true;
  end
  
  % reverse conjugate 
  % (change to f1 to get same as normxcorr2
  f2 = conj(f2(mb:-1:1,nb:-1:1));    
  
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
  M1  = fft2(m1,m,n);
  M2  = fft2(m2,m,n);     
  F1s = fft2(f1.^2,m,n);
  F2s = fft2(f2.^2,m,n);       
  
  ff = ifft2(F1.*F2);
  fsm = ifft2(F1s.*M2);
  mfs = ifft2(M1.*F2s);
  mm = ifft2(M1.*M2);  
    
  mqd  = (fsm-2*ff+mfs)./mm;  
  %mqd = (fsm-2*ff+mfs)/min(ma,mb)/min(na,nb);
  mqd = real(mqd);
  
  % trim to standard size
  if (pad)
    mqd(ma+mb:m,:) = [];
    mqd(:,na+nb:n) = [];
  else
    % only works for even min(ma,mb) & min(na,nb).    
    r = m - min(ma,mb)/2;
    s = n - min(na,nb)/2;
    mqd = circshift(mqd,r,1);
    mqd = circshift(mqd,s,2);
    mqd = mqd(1:end-1,1:end-1);    
    %??
  end
   
  if(nargout>1)
    
    mm = mm/min(ma,mb)/min(na,nb);
    if(pad)
      mm(ma+mb:m,:) = [];
      mm(:,na+nb:n) = [];
    else
      mm = circshift(mm,r,1);
      mm = circshift(mm,s,2);
      mm = mm(1:end-1,1:end-1);               
    end
  end
