function [x0,delta,out] = subpixel3x3(F)
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
% See also SUBPIXEL3X3LS, SUBPIXEL3X3LM, SUBPIXEL5X5LS, SUBPIXEL3X2, SUBPIXELNONE
 
  persistent invP

  if(isempty(invP))
    [x,y]=meshgrid(-1:1);
    P = bivander(x,y,2,2);  
    invP = inv(P);
  end
                
  [M,N] = size(F);
      
  [~,idx] = sort(F(:),'descend');        
  kmax = M*N/2;
  
  for k=1:kmax
    %[m,n] = ind2sub([M,N],idx(k));       
    m = rem(idx(k)-1,M)+1;
    n = (idx(k)-m)/M+1;
    
    if(m==1 || m==M || n==1 || n==N)
      continue;
    end
    
   
    % Fit biquadratic polynomial       
    f = (F(m-1:m+1,n-1:n+1));       
    %c = P\f(:);     
    c = invP*f(:);
     
    % Use Newtons method to find minimum
    J = -[c(8); c(6)];
    H = -[2*c(7) c(5); c(5) 2*c(3)];
    
    %[L,D] = ldl(H);
    %[U,S,V] = svd(D);
    
    [R,p] = chol(H);
    
    % Check if H is positive definite
    if (~p)
      delta = R\(R'\J);    
      % Check if it seems to converge
      if (max(abs(delta))<= 1.0)
        % Peak found
        %x0 = [n - (N+1)/2; m-(M+1)/2];        
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
  %detla = [0; 0];
  %[m,n] = ind2sub([M,N],idx(1));          
  %x0 = [n; m]; 
    