function [M,N,x,y,tx,ty] = pivgrid(MM,NN,opt)
% calc piv grid and bspline
     
    W = opt.Width;
    H = opt.Height;    
  
    dx = floor(W*(1-opt.Overlap));
    dy = floor(H*(1-opt.Overlap)); 
  
    M = ceil(MM/dy)+1;
    N = ceil(NN/dx)+1 ; 
  
    % Piv grid
    y0 = (MM - (M-1)*dy + 1 + mod(H,2) + mod(MM,2))/2;
    yy = y0 + (0:M-1)*dy;
    x0 = (NN - (N-1)*dx + 1 + mod(W,2) + mod(NN,2))/2;
    xx = x0 + (0:N-1)*dx;    
    [x,y] = meshgrid(xx,yy);
  
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