function [M,N,x,y,tx,ty] = pivgrid(MM,NN,opt)
% PIVGRID calulate placement of subwindows and bspline nodes for use in
% normalpass and distortedpass.
%
% [M,N,x,y,tx,ty] = pivgrid(MM,NN,opt)
% 
% input  - (MM,NN) - Size of images in image pair.
%        - opt     - PIV Options set using setpivopt(), 
%                    includes information about subwindow size and overlap.
% output - (M,N)   - Size of subwindow grid.
%        - (x,y)   - Center position of subwindows.
%        - (tx,ty) - Bspline nodes for least squares fitting.
%        
% See also setpivopt, normalpass, distortedpass
     
    W = min(opt.Width);
    H = min(opt.Height);    
  
    %dx = floor(W*(1-opt.Overlap));
    %dy = floor(H*(1-opt.Overlap)); 
    dx = floor(W*(1-opt.Overlap(1)));
    dy = floor(H*(1-opt.Overlap(1))); 
    if(dx==0 || dy==0)
        error('Too high overlap, or too small subwindows');
    end
       
    M = ceil(MM/dy)+1;
    N = ceil(NN/dx)+1; 
  
    % Piv grid
    tmp = (MM-1)-(M-1)*dy;    
    y0 = (tmp+mod(tmp,2))/2 + (1+mod(H,2))/2;
    yy = y0 + (0:M-1)*dy;    
    tmp = (NN-1)-(N-1)*dx; 
    x0 = (tmp+mod(tmp,2))/2 + (1+mod(W,2))/2;
    xx = x0 + (0:N-1)*dx;    
    [x,y] = meshgrid(xx,yy);
  
    if(length(opt.Overlap)>1)                 % beta       
        H = H*(1-opt.Overlap(2));
        W = W*(1-opt.Overlap(2));       
    end
               
    % B-spline nodes  
    I = ceil(MM/H)+1;
    J = ceil(NN/W)+1;
    i0 = (MM - (I-1)*H + 1 + mod(H,2))/2;
    ii = i0 + (0:I-1)*H;
    j0 = (NN - (J-1)*W + 1 + mod(W,2))/2;
    jj = j0 + (0:J-1)*W;      
    
    ty = zeros(I+4,1);
    ty(1:4) = y(1,1);
    ty(5:I) = ii(3:I-2);
    ty(I+1:I+4) = y(end,1);
  
    tx = zeros(J+4,1);
    tx(1:4) = x(1,1);
    tx(5:J) = jj(3:J-2);
    tx(J+1:J+4) = x(1,end);
end