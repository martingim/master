function piv = normalpasspar(piv,im1,mask1,im2,mask2,opt)
% NORMALPASS Add two values together.
%   piv = NORMALPASS([],im1,[],im2,[],opt) Single pass piv
%   piv = NORMALPASS([],im1,mask1,im2,mask2,opt) with masking.
%   piv = NORMALPASS(piv,im1,[],im2,[],opt) Shifted pass.
%   piv = NORMALPASS([],im1,[],im2,[],pive) Ensemble piv.
%
%   examples
%   
%   
%   See also NORMALPASS, INITPASS, SETPIVOPT.

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
  
  % Ensure single
  im1 = im2single(im1);  
  im2 = im2single(im2);
  mask1 = im2single(mask1);
  mask2 = im2single(mask2);
  %im1 = im2double(im1);  
  %im2 = im2double(im2);
  %mask1 = im2double(mask1);
  %mask2 = im2double(mask2);
  
  % check if im2 has same size?
  
  W = opt.Width;
  H = opt.Height;
      
  im1   = padarray(im1,[H W]);     
  mask1 = padarray(mask1,[H W]);
  im2   = padarray(im2,[H W]);     
  mask2 = padarray(mask2,[H W]);       
  MM = MM + 2*H;
  NN = NN + 2*W;
  
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
    %lastupdate=toc;
    
    n_blocks = ceil(M*N/2048/2);
    processed_subwindows = 0;
    for b=1:n_blocks
        % select block of subwindows
        block = b:n_blocks:M*N;
        block_subwindows = length(block);
        processed_subwindows = processed_subwindows + block_subwindows;
        
        % Define subwindow         
        idx1 = bsxfun(@plus,(1:W)' - .5 + W/2, piv.x(block)); 
        idy1 = bsxfun(@plus,(1:H)' - .5 + H/2, piv.y(block));            
        
        %size(idx1)
        %size(U(block)')
        if(piv.passes > 1)
          idx2 = bsxfun(@plus, idx1, U(block));
          idy2 = bsxfun(@plus, idy1, V(block));
          idx1 = bsxfun(@minus, idx1, U(block));
          idy1 = bsxfun(@minus, idy1, V(block));                            
        else
          idx2 = idx1;
          idy2 = idy1;
        end                                           
        
        
        % check if indices are out of bound (Need to test this) 
        ido = idx1(1,:) < 1 | idy1(1,:) < 1 | idx1(W,:) > NN | idy1(H,:) > MM | ...
              idx2(1,:) < 1 | idy2(1,:) < 1 | idx2(W,:) > NN | idy2(H,:) > MM;
                  
        U(block(ido)) = 2*U(block(ido)); 
        V(block(ido)) = 2*V(block(ido));
            
        piv.peak(block(ido)) = 0;
        piv.snr(block(ido)) = 0;                    
        piv.masked(block(ido)) = 0;  
        out1(block(ido)) = false;            
        
        idx1 = idx1(:,~ido);
        idy1 = idy1(:,~ido);
        idx2 = idx2(:,~ido);
        idy2 = idy2(:,~ido);
               
        idx1 = reshape(idx1,1,W,block_subwindows);
        idy1 = reshape(idy1,H,1,block_subwindows);
        idx2 = reshape(idx2,1,W,block_subwindows);
        idy2 = reshape(idy2,H,1,block_subwindows);
                        
        % Extract subwindows
        idx1 = bsxfun(@plus,(idx1-1)*MM,idy1);
        idx2 = bsxfun(@plus,(idx2-1)*MM,idy2);                    
        m1 = bsxfun(@times,mask1(idx1),window);
        m2 = bsxfun(@times,mask2(idx2),window);                              
        
        
        f1 = im1(idx1);                
        f2 = im2(idx2);  
                    
        % check masking before calling measure?
        %assignin('base','f1',f1)
        %assignin('base','f2',f2)
        %assignin('base','m1',m1)
        %assignin('base','m2',m2)
        
        % if(K>1)
        % reshape((0:K-1)*MM*NN,1,1,1,K)
        % end
  
        % Difference measure    
        %F = zero(1)
        %for k=1:K
        %end
        [F,mm] = opt.measure(f1,m1,f2,m2,idh,idw,opt.mparam);  
                    
        %assignin('base','F',F)
        %assignin('base','mm',mm)
              
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
                
                U(block(a)) = 2*U(block(a)); 
                V(block(a)) = 2*V(block(a));
                        
                piv.peak(block(a)) = 0;
                piv.snr(block(a)) = 0;                    
                piv.masked(block(a)) = 0;  
                out1(block(a)) = false;       
            else       
                
                % Find peak and do sub-interpolation     
                subF = F(idh,idw,a).*(mm(idh,idw,a)>=1/16); 
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
                U(block(a)) = 2*U(block(a)) - (x0(1) - delta(1) +opt.Umin-1); 
                V(block(a)) = 2*V(block(a)) - (x0(2) - delta(2) +opt.Vmin-1); 
            end
                
        end
        
        
        
        
          fprintf(repmat('\b',1,2*nl+4));
          fprintf(sprintf('%%%dd of %%%dd',nl,nl), processed_subwindows, M*N);
        
        
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
    
