function piv = distortedpass(piv,im1,mask1,im2,mask2,opt)
% DISTORTEDPASS Add two values together.
%   piv = DISTORTEDPASS(piv,im1,mask1,im2,mask2,opt) adds A to itself.
%   piv = DISTORTEDPASS(piv,im1,mask1,im2,mask2,pive) Ensemble.
%
%   input
%   output
%   piv
% -------------------------------------------------------------------------
% Example: Shallow water wave
% -------------------------------------------------------------------------
% % Read in wave images
% im1 = imread('wave1.png');
% mask1 = imread('wave1mask.png');
% im2 = imread('wave2.png');
% mask2 = imread('wave2mask.png');
%   
% % First pass with masks
% opt1 = setpivopt('range',[-8 8 -8 8],'subwindow',64,64,.75);
% piv1 = normalpass([],im1,mask1,im2,mask2,opt1); 
%  
% %Shifted pass using smaller subwindows
% opt2 = setpivopt('range',[-4 4 -4 4],'subwindow',32,32,.75);
% piv2 = normalpass(piv1,im1,mask1,im2,mask2,opt2);
%
% % Distorted passes
% opt3 = setpivopt('range',[-4 4 -4 4],'subwindow',32,32,.75);
% piv3 = distortedpass(piv2,im1,mask1,im2,mask2,opt3);
%  
% opt4 = setpivopt('range',[-4 4 -4 4],'subwindow',32,32,.75);
% piv4 = distortedpass(piv3,im1,mask1,im2,mask2,opt4);
%
% % Coordinate transformation
% load wavetform
% dt = 1/375;
% [Up,Vp,xp,yp] = replaceoutliers(piv4,mask1&mask2);     
% [Uw,Vw,xw,yw] = pixel2world(tform,Up,Vp,x,y,dt);
% 
% % Plot result
% figure;
% idx = 1:4:129;
% scale = 1/10;
% h= 0.2;
% quiver(xw(idx,idx)/h,yw(idx,idx)/h,scale*Uw(idx,idx),scale*Vw(idx,idx),0);
% xlabel(' x/h ');
% ylabel(' y/h ');
% 
% See also NORMALPASS, INITPASS, SETPIVOPT.

% function piv = distortedpass(piv,im1pp,im2pp,opt)
% see SETPIVOPT
  if (nargin < 6)
    error('Need at least six input parameters!');
  end
  
  [MM,NN,K] = size(im1);  
  piv.passes = piv.passes + 1;
  if(isfield(opt,'pass'))
      pive = opt;
      opt = pive.opt{end};
      piv.ensemblepasses = pive.ensemblepasses + 1;             
      if(~isfield(pive,'F'))
          error(['The difference measure peak from the previous pass needs to be stored when using ensemble PIV,' ...
              ' use setpivopt(...,''savepeaks'',true);']);
      end                
  else
      pive = [];            
      piv.ensemblepasses = 1;
  end    
  if(K>1 || ~isempty(pive))
    piv.pass{piv.passes} = 'distorted (ensemble)';
  else      
    piv.pass{piv.passes} = 'distorted';
  end      
  piv.opt{piv.passes} = opt;
  if(isempty(pive))  
    msg = sprintf('PIV %s pass %d: ',piv.pass{piv.passes},piv.passes);   
  else    
    msg = sprintf('PIV %s pass %d-%d: ',piv.pass{piv.passes},piv.passes,piv.ensemblepasses);   
  end
  
  [M,N,piv.x,piv.y,tx,ty] = pivgrid(MM,NN,opt);
   if(isempty(mask1))
    mask1 = ones(MM,NN);
  end
  if(isempty(mask2))
    mask2 = ones(MM,NN);
  end  
  %assignin('base','masktest',mask1);
  
  % Ensure double  
  im1 = im2double(im1);  
  im2 = im2double(im2);
  mask1 = im2double(mask1);
  mask2 = im2double(mask2);
  
  W = opt.Width;
  H = opt.Height;          
  window = bsxfun(@times,opt.window(H),opt.window(W)');
  window = window/mean(window(:));
  
  U = zeros(M,N);
  V = zeros(M,N);
  out1 = true(M,N);
  if(~strcmp(opt.snr,'none'))
    piv.peak = zeros(M,N);
    piv.snr = zeros(M,N);
  else
     piv.peak = [];
     piv.snr = [];
  end
  piv.masked = zeros(M,N); 
  if(opt.savepeaks)
    piv.F =  zeros(2*H-1,2*W-1,M,N);   
  end
            
  tic;     
    fprintf('%s Distort images\n', msg);
    
    if(isempty(pive) || ~opt.savedistortions)
        [x0,y0] = meshgrid(1-W:NN+W,1-H:MM+H);     
        [xa,ya,xb,yb] = opt.rkfun(opt.tableau,piv,x0,y0);                    
    else         
        xa = pive.distortions.xa;
        ya = pive.distortions.ya;
        xb = pive.distortions.xb;
        yb = pive.distortions.yb;
    end
    if(opt.savedistortions)        
        piv.distortions.xa = xa;
        piv.distortions.ya = ya;
        piv.distortions.xb = xb;
        piv.distortions.yb = yb;
    end
        
    [im1,mask1,im2,mask2] = opt.iminterp(im1,mask1,xa,ya,im2,mask2,xb,yb,opt.iparam);    
    
    %piv.pixeldiff = im1-im2;
    
    
    %switch(opt.iminterp)
    %  case 'lanczos'
    %    [im1,mask1,im2,mask2] = iminterp2lanczosj(im1,mask1,xa,ya,im2,mask2,xb,yb);           
    %  otherwise    
    %    %[im1,mask1,im2,mask2] = iminterp2bsplinej(im1,mask1,xa,ya,im2,mask2,xb,yb);   
    %    [im1,mask1,im2,mask2] = iminterp2bspline(im1,mask1,xa,ya,im2,mask2,xb,yb);                      
    %end     
    
    
  toc
  %assignin('base','im1c',im1); 
  %assignin('base','im2c',im2);
  %assignin('base','xa2',xa); 
  %assignin('base','ya2',ya);
  %assignin('base','xb2',xb); 
  %assignin('base','yb2',yb);
        
  
  %return
  idw = (opt.Umin:opt.Umax) + W;
  idh = (opt.Vmin:opt.Vmax) + H;
    
  tic;    
    nl = size(sprintf('%d',M*N),2);        
    fprintf(sprintf('%%s %%%dd of %%%dd',nl,nl), msg, 1, M*N);
    lastupdate=toc;      
    
    for m=1:M
      for n=1:N    
                   
        idx = (1:W) + piv.x(m,n) - .5 + W/2; 
        idy = (1:H) + piv.y(m,n) - .5 + H/2;
        
        %f1 = im1(idy,idx); 
        %m1 = mask1(idy,idx);
        %f2 = im2(idy,idx);
        %m2 = mask2(idy,idx);
                      
                                                        
        % Extract subwindow
        m1 = mask1(idy,idx).*window;
        m2 = mask2(idy,idx).*window;
        f1 = im1(idy,idx,:);                
        f2 = im2(idy,idx,:);   
            
        % Difference measure  
        [F,mm] = opt.measure(f1,m1,f2,m2,idh,idw,opt.mparam);  
        
        F = opt.prefun(opt.ensfun(F,3));         
        %[F,mm] = opt.measure(f1,m1,f2,m2,idh,idw,opt.mparam);          
        if(~isempty(pive))                       
           F = (F + pive.ensemblepasses*pive.F(:,:,m,n))/piv.ensemblepasses;
        end
        if(opt.savepeaks)
           piv.F(:,:,m,n) = F;
        end
                                            
        % Find peak and do sub-interpolation
        if(max(mm(:))<1/16)
            U(m,n) = 0; 
            V(m,n) = 0;
                        
            piv.peak(m,n) = 0;
            piv.snr(m,n) = 0;                    
            piv.masked(m,n) = 0;  
            out1(m,n) = false;       
        else
            subF = F(idh,idw); %subF = opt.prefun(opt.ensfun(F(idh,idw,:),3));     
            %subF = opt.prefun(F(idh,idw));
            [x0,delta,out1(m,n)] = opt.subpixel(subF);    
            switch (opt.snr)
              case 'sub'
                tmp = nanmedian(subF(:));        
                %piv.peak(m,n) =  abs(interp2(subF,x0(1)-delta(1),x0(2)-delta(2)) -tmp);
                piv.peak(m,n) =  abs(myinterp2(subF,x0,delta) - tmp);
                piv.snr(m,n) = piv.peak(m,n)/nanmedian(abs(subF(:)-tmp));
              case 'full'
                tmp = nanmedian(F(:));        
                %piv.peak(m,n) =  abs(interp2(subF,x0(1)-delta(1),x0(2)-delta(2)) -tmp);
                piv.peak(m,n) =  abs(myinterp2(subF,x0,delta) - tmp);
                piv.snr(m,n) = piv.peak(m,n)/nanmedian(abs(F(:)-tmp)); 
              otherwise                        
                    % ...
            end     
           
            %tmp = median(subF(:));        
            
            %piv.peak(m,n) =  abs(interp2(subF,x0(1)-delta(1),x0(2)-delta(2)) -tmp);
            %piv.snr(m,n) = piv.peak(m,n)/median(abs(subF(:)-tmp));
               
            %piv.masked(m,n) = interp2(mm(idh,idw),x0(1)-delta(1),x0(2)-delta(2));               
            piv.masked(m,n) = myinterp2(mm(idh,idw),x0,delta);
        
             
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
    
    out = piv.masked>.2 & out1;
    
    piv.U = opt.alpha*piv.Uspline.evaluate(piv.x,piv.y) + U; 
    piv.V = opt.alpha*piv.Vspline.evaluate(piv.x,piv.y) + V;        
           
    [piv.Uspline,piv.Vspline,piv.localres,piv.globalres,piv.out] = lsbsfit2(piv.U,piv.V,out,piv.x,piv.y,tx,ty);
                     
    fprintf('\n');
  toc
 
    
function fi = myinterp2(f,x0,delta)
% MYINTERP2 faster linear interpolation for single value
    [I,J] = size(f);
    if(x0(1)==1 || x0(2)==1 || x0(1)==J || x0(2)==I)
        fi = 0;
    else 
        idx = (-1:1);    
        f = f(idx+x0(2),idx+x0(1));     
        dx = idx'*ones(1,3) + delta(1);
        dy = ones(3,1)*idx  + delta(2);    
        w = max(1-abs(dx),0).*max(1-abs(dy),0);    
        fi = sum(f(:).*w(:));
    end
  