function [x0,delta,out] = subpixelnone(F)
% SUBPIXELNONE finds highest peak position of difference measure F 
% with no subpixel interpolation.
% 
% [x0,delta,out] = SUBPIXELNONE(F)
% 
% F     - Difference measure with positve peak
% x0    - Integer peak position
% delta - Allways returns zero, ie no subpixel interpolation
% out   - Allways returns one (valid peak)
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
% [x0,delta,out] = SUBPIXELNONE(F);
% 
% figure;
% imagesc(F);
% hold on; 
% plot(x0(1)-delta(1),x0(2)-delta(2),'r+')
%
% See also SUBPIXEL3X3, SUBPIXEL3X3LS, SUBPIXEL3X3LM, SUBPIXEL5X5LS, SUBPIXEL3X2

  % Find maximum value not on the edge               
  [f,m] = max(F(2:end-1,2:end-1));
  [~,n] = max(f);  
  
  x0 = [n+1; m(n)+1];
  delta = [0; 0];
  out = 1;  