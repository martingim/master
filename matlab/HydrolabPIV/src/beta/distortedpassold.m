function piv = distortedpass(piv,im1,im2,opt)
% function piv = distortedpass(piv,im1pp,im2pp,opt)
  
  piv.passes = piv.passes + 1;
  piv.pass{piv.passes} = 'distorted'; 
  
  
  if (nargin < 4)
    error('Need at least four input parameters!');
  end
  
  [MM,NN] = size(im1);   
  
  %[M,N] = size(piv.x);         
  %[x,y] = meshgrid(1:MM,1:NN);    
    
  %if ~isa(class(im1bs),'Bspline2')
  %end
  
  piv.opt{piv.passes} = opt;
  
  W = opt.Width;
  H = opt.Height;    
    
  dx = floor(W*(1-opt.Overlap));
  dy = floor(H*(1-opt.Overlap));  
  
  M = ceil(MM/dy)+1;
  N = ceil(NN/dx)+1;  
  
  % Piv grid
  %y0 = (MM - (M-1)*dy + 1 + mod(H,2))/2;  
  y0 = (MM - (M-1)*dy + 1 + mod(H,2) + mod(MM,2))/2;
  yy = y0 + (0:M-1)*dy;
  %x0 = (NN - (N-1)*dx + 1 + mod(W,2))/2;  
  x0 = (NN - (N-1)*dx + 1 + mod(W,2) + mod(NN,2))/2;
  xx = x0 + (0:N-1)*dx;    
  [piv.x,piv.y] = meshgrid(xx,yy);
      
  %piv.U = zeros(M,N);
  %piv.V = zeros(M,N);
  U = zeros(M,N);
  V = zeros(M,N);
  out1 = ones(M,N);
      
  % Create window
  %[x,y] = meshgrid(-(W-1)/2:(W-1)/2,-(H-1)/2:(H-1)/2);

  %a = 2.5;
  %window = exp(-1/2*(a*x/W/2).^2).*exp(-1/2*(a*y/H/2).^2);
  
  msg = sprintf('PIV pass %d (distorted): ',piv.passes);   
  tic; 
    if(false)      
      fprintf('%s Distort images\n', msg);
      %Uspline = piv.Uspline; 
      %Vspline = piv.Vspline; 

      % Distort images using previous pass  
      h = 1/2;      
      [x0,y0] = meshgrid(1-W:NN+W,1-H:MM+H);
      kx1 = h*piv.Uspline.evaluate(x0,y0);
      ky1 = h*piv.Vspline.evaluate(x0,y0);      
      xa = x0 - kx1;
      ya = y0 - ky1;         
      xb = x0 + kx1;
      yb = y0 + ky1;      
    else
      fprintf('%s Distort images\n', msg);
       % Distort images using previous pass  
      h = 1/2;  
      [x0,y0] = meshgrid(1-W:NN+W,1-H:MM+H);       
      kx1 = h*piv.Uspline.evaluate(x0,y0);
      ky1 = h*piv.Vspline.evaluate(x0,y0);
      kx2 = h*piv.Uspline.evaluate(x0-kx1/2,y0-ky1/2);
      ky2 = h*piv.Vspline.evaluate(x0-kx1/2,y0-ky1/2);
      kx3 = h*piv.Uspline.evaluate(x0-kx2/2,y0-ky2/2);
      ky3 = h*piv.Vspline.evaluate(x0-kx2/2,y0-ky2/2);
      kx4 = h*piv.Uspline.evaluate(x0-kx3,y0-ky3);
      ky4 = h*piv.Vspline.evaluate(x0-kx3,y0-ky3);
      xa = x0 - (kx1+2*kx2+2*kx3+kx4)/6;
      ya = y0 - (ky1+2*ky2+2*ky3+ky4)/6;   
                    
      %kx1 = h*piv.Uspline.evaluate(x0,y0);
      %ky1 = h*piv.Vspline.evaluate(x0,y0);
      kx2 = h*piv.Uspline.evaluate(x0+kx1/2,y0+ky1/2);
      ky2 = h*piv.Vspline.evaluate(x0+kx1/2,y0+ky1/2);
      kx3 = h*piv.Uspline.evaluate(x0+kx2/2,y0+ky2/2);
      ky3 = h*piv.Vspline.evaluate(x0+kx2/2,y0+ky2/2);
      kx4 = h*piv.Uspline.evaluate(x0+kx3,y0+ky3);
      ky4 = h*piv.Vspline.evaluate(x0+kx3,y0+ky3);
      xb = x0 + (kx1+2*kx2+2*kx3+kx4)/6;
      yb = y0 + (ky1+2*ky2+2*ky3+ky4)/6;      
      
    end
    
    switch(opt.iminterp)
      case 'lanczos'
        %[im1,mask1] = lanczos2(im1,xa,ya,opt.iparam);
        %[im2,mask2] = lanczos2(im2,xb,yb,opt.iparam);
        im1l = Lanczos2(im1);
        im1 = im1l.evaluate(xa,ya,opt.iparam);
        mask1 = im1l.mask();
        im2l = Lanczos2(im2);
        im2 = im2l.evaluate(xb,yb,opt.iparam);
        mask2 = im2l.mask();            
      otherwise                    
        [im1bs,im2bs] = im2bspline(im1,im2);
        im1 = im1bs.evaluate(xa,ya);
        mask1 = im1bs.mask();    
        im2 = im2bs.evaluate(xb,yb);
        mask2 = im2bs.mask();     
    end
    
 
     
  toc
  %assignin('base','im1b',im1); 
  %assignin('base','im2b',im2);
  
  
  idw = (opt.Umin:opt.Umax) + W;
  idh = (opt.Vmin:opt.Vmax) + H;
    
  tic;    
    nl = size(sprintf('%d',M*N),2);        
    fprintf(sprintf('%%s %%%dd of %%%dd',nl,nl), msg, 1, M*N);
    lastupdate=toc;
  
    %figure;
    %imagesc(im1)
    %figure;
    %imagesc(im2)
    
    
    for m=1:M
      for n=1:N    
                   
        idx = (1:H) + piv.x(m,n) - .5 + H/2; % Need to fix ?
        idy = (1:W) + piv.y(m,n) - .5 + W/2;
        
        f1 = im1(idy,idx); 
        m1 = mask1(idy,idx);
        f2 = im2(idy,idx);
        m2 = mask2(idy,idx);
        
        
        % Difference measure                        
        [F,mm] = opt.measure(f1,m1,f2,m2,idh,idw,opt.mparam);          
    
                                            
        % Find peak and do sub-interpolation
        if(max(mm(:))<1/16)
            U(m,n) = 0; 
            V(m,n) = 0;
                        
            piv.peak(m,n) = 0;
            piv.snr(m,n) = 0;                    
            piv.masked(m,n) = 0;  
            out1(m,n) = 0;       
        else
        subF = opt.prefun(F(idh,idw));
        [x0,delta,out1(m,n)] = opt.subpixel(subF);              
        tmp = median(subF(:));        
        piv.peak(m,n) =  abs(interp2(subF,x0(1)-delta(1),x0(2)-delta(2)) -tmp);
        piv.snr(m,n) = piv.peak(m,n)/median(abs(subF(:)-tmp));
        %[x0,delta,out1(m,n)] = opt.subpixel(opt.prefun(F(idh,idw)));    
        
        piv.masked(m,n) = interp2(mm(idh,idw),x0(1)-delta(1),x0(2)-delta(2));
        
        if(false)
          figure(3)
          [i,j]= size(F);
          imagesc(F,'Xdata',[-1 1]*(j-1)/2,'Ydata',[-1 1]*(i-1)/2);                
          colorbar
          assignin('base','F',F);
          pause;
        end        
        
        U(m,n) = -(x0(1) - delta(1)+opt.Umin-1);
        V(m,n) = -(x0(2) - delta(2)+opt.Vmin-1);
        end      
        
        % Update every half second
        if (toc - lastupdate > .5 || (m==M && n==N))
          fprintf(repmat('\b',1,2*nl+4));
          fprintf(sprintf('%%%dd of %%%dd',nl,nl), N*(m-1)+n, M*N);
          lastupdate = toc;
        end
      end
    end
    
    U = piv.Uspline.evaluate(piv.x,piv.y) + opt.alpha*U;
    V = piv.Vspline.evaluate(piv.x,piv.y) + opt.alpha*V;        
    
    
    % Outlier detection            
    %[piv.globalres,out3] = globalfilt(U,V,Inf,5);    
    [piv.extres,out3] = localfilt(piv.x,piv.y,U,V,0.1,2,9);  
    [piv.localres,out2] = localfilt(piv.x,piv.y,U,V,0.1,2,3);    
    %[~,out2] = outlier([piv.x(:) piv.y(:)],[U(:) V(:)],0.1,9,2);
    %out2 = reshape(out2,M,N);                
    out = out1 & out2 & out3;
    piv.out = out;
    
    % LS B-spline fit
    %if(opt.ignoreoutliers)
    %  U = U(out);
    %  V = V(out);
    %  x = piv.x(out);
    %  y = piv.y(out);    
    %else
    %  U = U(:);
    %  V = V(:);
    %  x = piv.x(:);
    %  y = piv.y(:);    
    %end
    
    % B-spline nodes  
    I = ceil(MM/H)+1;
    J = ceil(MM/W)+1;
    i0 = (MM - (I-1)*H + 1 + mod(H,2))/2;
    ii = i0 + (0:I-1)*H;
    j0 = (NN - (J-1)*W + 1 + mod(W,2))/2;
    jj = j0 + (0:J-1)*W; 
    
  
    ty = zeros(I+4,1);
    ty(1:4) = piv.y(1,1);
    ty(5:I) = ii(3:I-2);
    ty(I+1:I+4) = piv.y(end,1);
    %piv.bspline.ty = ty;
  
    tx = zeros(J+4,1);
    tx(1:4) = piv.x(1,1);
    tx(5:J) = jj(3:J-2);
    tx(J+1:J+4) = piv.x(1,end); 
    %tx = piv.bspline.tx;
    %ty = piv.bspline.ty;
    
    %assignin('base','out1',out1)
    %[piv.Uspline,piv.Vspline,piv.localres,piv.globalres] = LSBSfit2m2(U,V,out1,piv.x,piv.y,tx,ty);
    
    %r1 = max(piv.localres,[],3);
    %r2 = max(piv.globalres,[],3);
    % change next line?
    %piv.out = (r1<2) & abs(r2) < 4.685/.6745*median(r2(:));
    
    [piv.Uspline,piv.Vspline,piv.LSres] = LSBSfit2(U,V,out,piv.x,piv.y,tx,ty);              
    
    piv.out = piv.out & piv.LSres <= 9;
    
    piv.U = U;
    piv.V = V;   
    
    %I = length(ty)-4;
    %J = length(tx)-4;
   % K = length(U);
   % 
    %A = sparse([],[],[],K,I*J,K*ceil(4/(1-opt.Overlap)));
    %for i=1:I   
    %  for j=1:J
    %    k = (j-1)*I + i;      
    %    A(:,k) = basis2(i,4,y,ty).*basis2(j,4,x,tx);      
    %  end
    %end
    
    %piv.bspline.Ucoeff = reshape(A\U,I,J);
    %piv.bspline.Vcoeff = reshape(A\V,I,J);
    %
    %piv.U = bseval2(piv.bspline.Ucoeff,tx,ty,piv.x,piv.y,4);
    %piv.V = bseval2(piv.bspline.Vcoeff,tx,ty,piv.x,piv.y,4);
    
    fprintf('\n');
  toc
 
    
  