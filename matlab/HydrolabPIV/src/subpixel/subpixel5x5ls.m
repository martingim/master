function [x0,delta,out] = subpixel5x5ls(F)
% SUBPIXEL5X5LS finds highest peak position of difference measure F 
% with a 5x5 point biquadratic least squares fit for subpixel interpolation.
% 
% [x0,delta,out] = SUBPIXEL5X5LS(F)
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
% [x0,delta,out] = SUBPIXEL5X5LS(F);
% 
% figure;
% imagesc(F);
% hold on; 
% plot(x0(1)-delta(1),x0(2)-delta(2),'r+')
%
% See also SUBPIXEL5X5LM, SUBPIXEL3X3, SUBPIXEL3X3LS, SUBPIXEL3X3LM, SUBPIXEL3X2, SUBPIXELNONE

  persistent P

  if(isempty(P))
    [x,y]=meshgrid(-2:2);
    P = bivander(x,y,2,2);   
  end
  
  [M,N] = size(F);
      
  [~,idx] = sort(F(:),'descend');        
  kmax = M*N/2;
     
  for k=1:kmax
    [m,n] = ind2sub([M,N],idx(k));
    
    if(m<=2 || m>=M-1 || n<=2 || n>=N-1)
      continue;
    end
      
    % Fit biquadratic polynomial       
    f = (F(m-2:m+2,n-2:n+2));           
    c  = P\f(:);
    
    % Use Newtons method to find maximum
    J = -[c(8); c(6)];
    H = -[2*c(7) c(5); c(5) 2*c(3)];
    [R,p] = chol(H);
    
    % (Need to check this)
    % H is positive definite if x0 is near a local minimum    
    if (~p)      
      delta = R\(R'\J);          
      if (max(abs(delta))<= 1.0)        
        % Peak found   
        x0 = [n; m];          
        out = 1;
        return;
      end 
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
  