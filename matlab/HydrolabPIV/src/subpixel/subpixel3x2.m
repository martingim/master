function [x0,delta,out] = subpixel3x2(F)
% SUBPIXEL3X2 finds highest peak position of difference measure F 
% with two 3 point biquadratic fit for subpixel interpolation.
% 
% [x0,delta,out] = SUBPIXEL3X2(F)
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
% [x0,delta,out] = SUBPIXEL3X2(F);
% 
% figure;
% imagesc(F);
% hold on; 
% plot(x0(1)-delta(1),x0(2)-delta(2),'r+')
%
% See also SUBPIXEL3X3, SUBPIXEL3X3LS, SUBPIXEL3X3LM, SUBPIXEL5X5LS, SUBPIXELNONE

  [M,N] = size(F);
      
  [~,idx] = sort(F(:),'descend');        
  kmax = M*N/2;
  
  for k=1:kmax
    [m,n] = ind2sub([M,N],idx(k));       
    
    if(m==1 || m==M || n==1 || n==N)
      continue;
    end
       
    % Fit biquadratic polynomial
    delta = zeros(2,1);
    delta(1) = -(F(m,n-1) - F(m,n+1))/(F(m,n-1)-2*F(m,n)+F(m,n+1))/2;
    delta(2) = -(F(m-1,n) - F(m+1,n))/(F(m-1,n)-2*F(m,n)+F(m+1,n))/2;
             
    if (max(abs(delta))<= 1.0)
        % Peak found      
        x0 = [n; m];  
        out = 1;
        return;
    end 
    
  end
  
  % No peak found!
  out = 0;
  tmp = [(N+1)/2; (M+1)/2];
  x0 = floor(tmp);
  delta = tmp-x0;
  %delta = [0;0];          
  %[m,n] = ind2sub([M,N],idx(1));        
  %x0 = [n; m]; 
    