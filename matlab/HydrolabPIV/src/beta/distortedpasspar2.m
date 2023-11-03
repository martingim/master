function piv = distortedpasspar2(piv,im1,mask1,im2,mask2,opt)
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
  
  W1 = opt.Width(1);
  W2 = opt.Width(end);
  H1 = opt.Height(1);
  H2 = opt.Height(end);  
  W = max(W1,W2);
  H = max(H1,H2); 
  
  window1 = bsxfun(@times,opt.window(H1),opt.window(W1)');
  window1 = window1/mean(window1(:));
  window2 = bsxfun(@times,opt.window(H2),opt.window(W2)');
  window2 = window2/mean(window2(:));
  
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
    piv.F =  zeros(H1+H2-1,W1+W2-1,M,N);   
  end
            
  tic;     
    fprintf('%s Distort images %3d%%', msg, 0);
       
    MH = MM+2*H;
    NH = NN+2*W;
    if(isempty(pive) || ~opt.savedistortions)
        [x0,y0] = meshgrid(1-W:NN+W,1-H:MM+H);   
        xa = zeros(MH,NH);    
        ya = zeros(MH,NH);    
        xb = zeros(MH,NH);    
        yb = zeros(MH,NH);
        
        L = ceil(MH*NH/65536);        
        for l=1:L
            idx = l:L:MH*NH;
            [xa(idx),ya(idx),xb(idx),yb(idx)] = opt.rkfun(opt.tableau,piv,x0(idx),y0(idx));                
            fprintf('\b\b\b\b%3d%%',floor(100*l/L));
        end
        
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
    fprintf('\n');
  toc
  %assignin('base','imi1',im1);
  %assignin('base','imi2',im2);
  %assignin('base','imask1',mask1);
  %assignin('base','imask2',mask2);
  
  %return
  idw = (opt.Umin:opt.Umax) + .5*(W1+W2);
  idh = (opt.Vmin:opt.Vmax) + .5*(H1+H2);
    
  tic;    
    nl = size(sprintf('%d',M*N),2);        
    fprintf(sprintf('%%s %%%dd of %%%dd',nl,nl), msg, 1, M*N);
        
    n_blocks = ceil(M*N/2048);
    processed_subwindows = 0;
    
    for b=1:n_blocks
        % select block of subwindows
        block = b:n_blocks:M*N;
        block_subwindows = length(block);
        processed_subwindows = processed_subwindows + block_subwindows;
        
        % Define subwindow      
        idx1 = bsxfun(@plus,(1:W1)' - .5 + W - W1/2, piv.x(block)); 
        idy1 = bsxfun(@plus,(1:H1)' - .5 + H - H1/2, piv.y(block));      
        idx2 = bsxfun(@plus,(1:W2)' - .5 + W - W2/2, piv.x(block)); 
        idy2 = bsxfun(@plus,(1:H2)' - .5 + H - H2/2, piv.y(block));                  
                       
        idx1 = reshape(idx1,1,W1,block_subwindows);
        idy1 = reshape(idy1,H1,1,block_subwindows);
        idx2 = reshape(idx2,1,W2,block_subwindows);
        idy2 = reshape(idy2,H2,1,block_subwindows);
        
%         assignin('base','idx1',idx1);
%         assignin('base','idx2',idx2);
%         assignin('base','idy1',idy1);
%         assignin('base','idy2',idy2); % seems correct
        
        % Extract subwindows
        idx1 = bsxfun(@plus,(idx1-1)*MH,idy1); %?
        idx2 = bsxfun(@plus,(idx2-1)*MH,idy2);                    
        m1 = bsxfun(@times,mask1(idx1),window1);
        m2 = bsxfun(@times,mask2(idx2),window2); 
        
        f1 = im1(idx1);                
        f2 = im2(idx2);  
        
    
%         
%         assignin('base','f1',f1); % not correct
%         assignin('base','f2',f2);
%         assignin('base','m1',m1);
%         assignin('base','m2',m2);
%         return
        [F,mm] = opt.measure(f1,m1,f2,m2,idh,idw,opt.mparam);  
              
        %F = opt.prefun(opt.ensfun(F,4));                        
        if(~isempty(pive))           
            % TODO: find a way of dealing with nans/infs
            % replace nan->0, [-Inf,-1)->-1, (1,Inf]->1 ?? (only ncc)               
            % weighted average?
            F = (F + pive.ensemblepasses*pive.F(:,:,block))/piv.ensemblepasses;
        end
        
        if(opt.savepeaks)
            piv.F(:,:,block) = F;
        end
        
                        
        maxmm = max(max(mm));
        for a=1:block_subwindows
            
            if(maxmm(a)<1/16)
                U(block(a)) = 0; 
                V(block(a)) = 0;
                        
                piv.peak(block(a)) = 0;
                piv.snr(block(a)) = 0;                    
                piv.masked(block(a)) = 0;  
                out1(block(a)) = false;       
            else       
                % Find peak and do sub-interpolation     
                subF = F(idh,idw,a); 
                [x0,delta,out1(block(a))] = opt.subpixel(subF);                
                switch (opt.snr)
                  case 'sub'
                    tmp = nanmedian(subF(:));        
                    piv.peak(block(a)) =  abs(myinterp2(subF,x0,delta) - tmp);
                    piv.snr(block(a)) = piv.peak(block(a))/nanmedian(abs(subF(:)-tmp));
                  case 'full'          
                    fullF = F(:,:,a); 
                    tmp = nanmedian(fullF(:));                            
                    piv.peak(block(a)) =  abs(myinterp2(subF,x0,delta) - tmp);                    
                    piv.snr(block(a)) = piv.peak(block(a))/nanmedian(abs(fullF(:)-tmp)); 
                  otherwise                 
                      
                    % ...
                end                
                
                piv.masked(block(a)) =  myinterp2(mm(idh,idw,a),x0,delta);
                U(block(a)) =  -(x0(1) - delta(1) +opt.Umin-1); 
                V(block(a)) =  -(x0(2) - delta(2) +opt.Vmin-1); 
            end
        end                       
        
        fprintf(repmat('\b',1,2*nl+4));
        fprintf(sprintf('%%%dd of %%%dd',nl,nl), processed_subwindows, M*N);                
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
  