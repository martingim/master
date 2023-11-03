function piv = normalpass(piv,im1,mask1,im2,mask2,opt)
% NORMALPASS performs PIV on a 
%
% piv = normalpass(piv,im1,mask1,im2,mask2,opt)
% piv = normalpass(piv,im1,mask1,im2,mask2,pive)
%
% input  - piv   - piv result from previos pass, 
%                  first pass use empty [].         
%        - im1   - First image
%        - mask1 - First mask, if empty mask i
%        - im2   - Second image
%        - mask2 - second mask
%        - opt   - see <a href="matlab:help setpivopt">setpivopt</a> 
%        - pive  - ensemble piv 
% output - piv   - 
%
% 
%
%
% -------------------------------------------------------------------------
% Example: Single pass PIV
% -------------------------------------------------------------------------
% piv = NORMALPASS([],im1,[],im2,[],opt) Single pass piv
%
% -------------------------------------------------------------------------
% Example: Multipass PIV (shifted pass)
% -------------------------------------------------------------------------
%
% piv = NORMALPASS(piv,im1,[],im2,[],opt) Shifted pass.
% -------------------------------------------------------------------------
% Example: Ensemble PIV using single call
% -------------------------------------------------------------------------
%   piv = NORMALPASS([],im1,mask1,im2,mask2,opt) with masking.
%
% -------------------------------------------------------------------------
% Example: Ensemble PIV using single call
% -------------------------------------------------------------------------
% I=512; J=512; K=25;
% im1 = zeros(I,J,K);
% im2 = zeros(I,J,K);
% for k=1:K
%   im1(:,:,k) = imread(sprintf('sparse%02dA.png',k));
%   im2(:,:,k) = imread(sprintf('sparse%02dB.png',k));  
% end
% 
% opt = setpivopt('ensemble',@nanmedian);
% piv = normalpass([],im1,[],im2,[],opt);
% [U,V,x,y] = replaceoutliers(piv);
% 
% figure;
% quiver(x,y,U,V);
%
% -------------------------------------------------------------------------
% Example: Ensemble piv using multiple calls 
% -------------------------------------------------------------------------
% opt1 = setpivopt('savepeaks',true);
% for k=1:25
%   im1 = imread(sprintf('sparse%02dA.png',k));
%   im2 = imread(sprintf('sparse%02dB.png',k));
%   
%   if(k==1)    
%     piv1 = normalpass([],im1,[],im2,[],opt1);
%   else
%     piv1 = normalpass([],im1,[],im2,[],piv1);
%   end
% end
% [U1,V1,x1,y1] = replaceoutliers(piv1);
% 
% figure;
% quiver(x1,y1,U1,V1);
%
%   
% See also NORMALPASS, INITPASS, SETPIVOPT.

%
%function piv = normalpass(piv,im1,mask1,im2,mask2,opt)
%function piv = normalpass(piv,im1,mask1,im2,mask2,pive)
    
  
% Notes
% allow indepent masks on ensembles??
  if (nargin < 6)
    error('Need at least six input parameters!');
  end
  
  
  [MM,NN,K] = size(im1);
  
  if (isempty(piv))
    piv.passes = 1;
  else
    piv.passes = piv.passes + 1;    
  end
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
    piv.pass{piv.passes} = 'normal (ensemble)';
  else      
    piv.pass{piv.passes} = 'normal';
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
  
  % Ensure double
  im1 = im2double(im1);  
  im2 = im2double(im2);
  mask1 = im2double(mask1);
  mask2 = im2double(mask2);
  
  %im1 = im2single(im1);  
  %im2 = im2single(im2);
  %mask1 = im2single(mask1);
  %mask2 = im2single(mask2);
  % check if im2 has same size?
  
  W = opt.Width;
  H = opt.Height;
      
  im1   = padarray(im1,[H W]);     
  mask1 = padarray(mask1,[H W]);
  im2   = padarray(im2,[H W]);     
  mask2 = padarray(mask2,[H W]);       
    
  window = bsxfun(@times,opt.window(H),opt.window(W)');
  window = window/mean(window(:));
  
  out1 = true(M,N);
  if(~strcmp(opt.snr,'none'))
    piv.peak = zeros(M,N);
    piv.snr = zeros(M,N);
  else
     piv.peak = [];
     piv.snr = [];
  end
  piv.masked = zeros(M,N);
      
  idw = (opt.Umin:opt.Umax) + W; % Adjust for subpixel ? 3x3 (+1)or 5x5 (+2)
  idh = (opt.Vmin:opt.Vmax) + H;
         
  piv.U = zeros(M,N);
  piv.V = zeros(M,N);  
  if(opt.savepeaks)
    piv.F =  zeros(2*H-1,2*W-1,M,N);   
  end
  
  if(piv.passes==1)
    U = zeros(M,N);
    V = zeros(M,N);
  else
    U = round(piv.Uspline.evaluate(piv.x,piv.y)/2);
    V = round(piv.Vspline.evaluate(piv.x,piv.y)/2);            
  end
  
  
  tic;
    nl = size(sprintf('%d',M*N),2);        
    fprintf(sprintf('%%s %%%dd of %%%dd',nl,nl), msg, 1, M*N);
    lastupdate=toc;
    for m=1:M
      for n=1:N        
                           
        % Define subwindow         
        idx1 = (1:W) + piv.x(m,n) - .5 + W/2; 
        idy1 = (1:H) + piv.y(m,n) - .5 + H/2;              
        if(piv.passes > 1)            
          idx2 = idx1 + U(m,n);
          idy2 = idy1 + V(m,n);
          idx1 = idx1 - U(m,n);
          idy1 = idy1 - V(m,n);                            
        else
          idx2 = idx1;
          idy2 = idy1;
        end                                           
        
        % check if indices are out of bound (Need to test this) 
        if(idx1(1) < 1 || idy1(1) < 1 || idx1(W) > NN+2*W || idy1(H) > MM+2*H || ...
           idx2(1) < 1 || idy2(1) < 1 || idx2(W) > NN+2*W || idy2(H) > MM+2*H)       
                  
            U(m,n) = 2*U(m,n); 
            V(m,n) = 2*V(m,n);
            
            piv.peak(m,n) = 0;
            piv.snr(m,n) = 0;                    
            piv.masked(m,n) = 0;  
            out1(m,n) = false;            
        else                                                                             
            % Extract subwindow
            m1 = mask1(idy1,idx1).*window;
            m2 = mask2(idy2,idx2);%.*window;    % only multiply m1??          
            f1 = im1(idy1,idx1,:);                
            f2 = im2(idy2,idx2,:);   
            
            % Difference measure  
            [F,mm] = opt.measure(f1,m1,f2,m2,idh,idw,opt.mparam);  
            
            F = opt.prefun(opt.ensfun(F,3));                        
            if(~isempty(pive))           
                % TODO: find a way of dealing with nans/infs
                % replace nan->0, [-Inf,-1)->-1, (1,Inf]->1 ?? (only ncc)               
                % weighted average?
                F = (F + pive.ensemblepasses*pive.F(:,:,m,n))/piv.ensemblepasses;
            end
            if(opt.savepeaks)
                piv.F(:,:,m,n) = F;
            end
    	    %[F,mm] = opt.measure(f1,m1,f2,m2,idh,idw,opt.mparam);  
            %assignin('base','F',F);
            %assignin('base',sprintf('f%d',k),f1);
            %assignin('base',sprintf('g%d',k),f2);
            %return
            
            % Find peak and do sub-interpolation
            if(max(mm(:))<1/16) % add value to options?
                U(m,n) = 2*U(m,n); 
                V(m,n) = 2*V(m,n);
                        
                piv.peak(m,n) = 0;
                piv.snr(m,n) = 0;                    
                piv.masked(m,n) = 0;  
                out1(m,n) = false;       
            else                                
                subF = F(idh,idw); 
                [x0,delta,out1(m,n)] = opt.subpixel(subF);        % does this work if Umax ~= -Umin?               
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
                  %case 'masked'
                  %  tmp = sum(F(:).*mm(:))/sum(mm(:));
                  %  piv.peak(m,n) =  abs(interp2(subF,x0(1)-delta(1),x0(2)-delta(2)) -tmp);
                  %  piv.snr(m,n) = piv.peak(m,n)/sum(abs(F(:)-tmp).*mm(:))*sum(mm(:)); 
                  otherwise                 
                      
                    % ...
                end                
                                
                %piv.masked(m,n) = interp2(mm(idh,idw),x0(1)-delta(1),x0(2)-delta(2));                
                piv.masked(m,n) = myinterp2(mm(idh,idw),x0,delta);
        
                U(m,n) = 2*U(m,n) - (x0(1) - delta(1) +opt.Umin-1); % Fix this in subpixel
                V(m,n) = 2*V(m,n) - (x0(2) - delta(2) +opt.Vmin-1); % instead ?
            end
        end
        
        % Update every half second
        if (toc - lastupdate > .5 || (m==M && n==N))
          fprintf(repmat('\b',1,2*nl+4));
          fprintf(sprintf('%%%dd of %%%dd',nl,nl), N*(m-1)+n, M*N);
          lastupdate = toc;
        end
      end
    end                  
     
    out = piv.masked>.2 & out1; % add option in setpivopt?
    
    [piv.Uspline,piv.Vspline,piv.localres,piv.globalres,piv.out] = lsbsfit2(U,V,out,piv.x,piv.y,tx,ty);
    piv.U = U;
    piv.V = V;
    
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
    
    
  
