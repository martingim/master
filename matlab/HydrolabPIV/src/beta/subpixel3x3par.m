function [x0,delta,out] = subpixel3x3par(F)
% SUBPIXEL3X3 finds highest peak position of difference measure F 
% with a 3x3 point biquadratic fit for subpixel interpolation.
% 
% [x0,delta,out] = SUBPIXEL3X3(F)
% 
% F     - Difference measure with positve peak
% x0    - Integer peak position
% delta - Subpixel correction for peak position
% out   - Returns one if the peak seems to be valid, otherwise zero.
% 
% Example
%
% im1 = im2double(imread('imA.png')); 
% im2 = im2double(imread('imB.png')); 
% A = im1(1:32,1:32);
% B = im2(1:32,1:32);
% M = ones(32);
%
% F = maskedncc(A,M,B,M);
% [x0,delta,out] = SUBPIXEL3X3(F);
% 
% figure;
% imagesc(F);
% hold on; 
% plot(x0(1)-delta(1),x0(2)-delta(2),'r+')
%
% See also SUBPIXEL3X3LS, SUBPIXEL3X3LM, SUBPIXEL5X5, SUBPIXEL3X2, SUBPIXELNONE
 
  persistent invP

  if(isempty(invP))
    [x,y]=meshgrid(-1:1); % ensure same class as F ?
    P = bivander(x,y,2,2);  
    invP = inv(P);
  end
                
  [M,N,K] = size(F);
  
  % find peaks not on edges
  G = -inf(M,N,K);
  G(2:M-1,2:N-1,:) =  F(2:M-1,2:N-1,:);    
  [~,idx] = sort(reshape(G,M*N,K),'descend');        
  kmax = M*N/2;
  
  % 
  out = zeros(1,K);
  delta = zeros(2,K);          
  %[m,n] = ind2sub([M,N],idx(1,:));   
  %x0 = [n; m]; 
  %x0 = floor([M/2; N/2]);% ??
  x0 = [ceil(M/2)*ones(1,K); ceil(N/2)*ones(1,K)];
  
  %assignin('base','idxx',idx);
  %assignin('base','out',out);
  %assignin('base','delta',delta);
  %assignin('base','x0',x0);
  
  
  for k=1:kmax
    %[m,n] = ind2sub([M,N],idx(k));   

    m = rem(idx(k,~out)-1,M)+1;
    n = (idx(k,~out)-m)/M+1;
    
    %size([)
    %size(x0(:,~out))
    x0(:,~out) = [n; m];           
    
    
    % Fit biquadratic polynomial       
    f = reshape(F(m-1:m+1,n-1:n+1,~out),9,K);       
    %c = P\f(:);     
    c = invP*f;
     
    %assignin('base','c',c);
    
    % Cholesky factorization of Hermitian
    L11 = sqrt(-2*c(7,:));
    L21 = -c(5,:)./L11;
    L22 = sqrt(-2*c(3,:)-L21.*L21);  
    
    % Check if H is positive definite    
    p = imag(L11)==0 & imag(L22)==0;
    pidx = false(1,K);
    pidx(~out) = p;
    %size(p)
    %size(L11)
    %size(L21)
    %size(L22)
    %size(c)
    % Solve with forward/backward substitution    
    tmp1 = -c(8,p)./L11(p);
    tmp2 = (-c(6,p) - L21(p).*tmp1)./L22(p);
    
    delta(2,pidx) = tmp2./L22(p);
    delta(1,pidx) = (tmp1 - L21(p).*delta(2,pidx))/L11(p);
    
    out(pidx) = max(delta(:,pidx))>1.0;
   % assignin('base','delta',delta)
   % assignin('base','x0',x0)
   %assignin('base','L11',L11)
   %assignin('base','L21',L21) 
   %assignin('base','L22',L22)
   %assignin('base','p',p)
   
  end
    