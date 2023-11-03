function [x0,delta,out] = subpixel3x3lm(F)
% SUBPIXEL3X3LM finds highest peak position of difference measure F 
% with a 3x3 point non-linear least squares fit (Levenberg–Marquardt) 
% to bivariate Gaussian.
% 
% [x0,delta,out] = SUBPIXEL3X3LM(F)
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
% [x0,delta,out] = SUBPIXEL3X3LM(F);
% 
% figure;
% imagesc(F);
% hold on; 
% plot(x0(1)-delta(1),x0(2)-delta(2),'r+')
% 
% See also SUBPIXEL3X3, SUBPIXEL3X3LS, SUBPIXEL5X5LS, SUBPIXEL5X5LM, SUBPIXEL3X2, SUBPIXELNONE
  

  [x,y]=meshgrid(-1:1);
  x=x(:);
  y=y(:);
  
  [M,N] = size(F);
      
  [~,idx] = sort(F(:),'descend');        
  kmax = M*N/2;
  
  for k=1:kmax
    [m,n] = ind2sub([M,N],idx(k));       
    
    if(m==1 || m==M || n==1 || n==N)
      continue;
    end
    
    f = (F(m-1:m+1,n-1:n+1));      
    %x0 = [n; m]; 
   
    % Intial guess
    c = [0 .5 8/9 0 8/9 0 0]';
    lambda = .5;
    for i=1:10
        [fi,J] = gfit(c,x(:),y(:));    
        dc = (J*J'+lambda*eye(7))\(J*(f(:)-fi));
        c = c + dc;
    end
    
    delta = -c(6:7);
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
  
function [f,J] = gfit(c,x,y)
  xx = (x-c(6)).*(x-c(6)); 
  xy = (x-c(6)).*(y-c(7));
  yy = (y-c(7)).*(y-c(7));
  p = c(3)*xx + 2*c(4)*xy + c(5)*yy;        
  f = c(2)*exp(-p);
                
  J = ones(7,9);
  J(2,:) = f/c(2);
  J(3,:) = -xx.*f;
  J(4,:) = -2*xy.*f;
  J(5,:) = -yy.*f;
  J(6,:) = 2*(c(3)*(x-c(6)) + c(4)*(y-c(7))).*f;
  J(7,:) = 2*(c(4)*(x-c(6)) + c(5)*(y-c(7))).*f;
     
    