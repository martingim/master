function [Ubspline,Vbspline,localres,globalres,out] = lsbsfit2(u,v,mask,x,y,tx,ty)
% LSBSFIT2 performs robust bsplines fit using iterative weighted least
% squares, the weights are determined with Tukey biweight using the residuals 
% normalized by the local median residual. The inital weights are set using
% a local median outlier filter.
%
% [Ubspline,Vbspline,localres,globalres,out] = LSBSFIT2(u,v,mask,x,y,tx,ty)
%
% input  - (u,v)             - Velocity vectors
%        - mask              - Masked vectors (must be logical)
%        - (x,y)             - Positions
%        - (tx,ty)           - Bspline nodes
% output - (Uspline,Vspline) - Fitted bspline
%        - localres          - Residuals for local median filter
%        - globalres         - Normalized residuals for bspline fit
%        - out               - Outliers based on localres and globalres
%
% See also normalpass, distortedpass, replaceoutliers, localfilt

    [M,N] = size(u);
    K = M*N;
    I = numel(ty)-4;
    J = numel(tx)-4;
    
    % Build Least squares matrix from Bspline basis functions
    Y = zeros(K,I);
    for i=1:I        
       Y(:,i) = Bspline2.basis(i,4,y(:),ty);
    end  
    Y = sparse(Y);
  
    X = zeros(K,J);
    for j=1:J       
       X(:,j) =  Bspline2.basis(j,4,x(:),tx);
    end
    X = sparse(X);
  
    idy = (1:I)'*ones(1,J);
    idx = ones(I,1)*(1:J);

    A  = Y(:,idy(:)).*X(:,idx(:)); 
    
    % Calculate residuals using local median filter
    tol = 0.1;
    rmax = 2;  
    
   % Local median outlier detection
   % Use 5x5 local filter if partially masked
   [r3,mi3] = localfilt(x,y,u,v,mask,tol,rmax,3);   
   [r5,~] = localfilt(x,y,u,v,mask,tol,rmax,5);
   localres = r5;
   idm = repmat(mi3,1,1,2)==8; 
   localres(idm) = r3(idm); 
   idm = repmat(mask,1,1,2);
   localres(~idm) = rmax;
   localres = max(localres,[],3)/rmax;
          
    % Use Local residuals as intial weights
    W = max(1-localres.^2/(1.171.^2),eps).^2;  
      
    % Weighted Least squares fit
    AW = bsxfun(@times,A,W(:))';   
    [L,isposdef] = chol(AW*A,'lower');
    if(~isposdef)
        cu = L'\(L\(AW*u(:)));
        cv = L'\(L\(AW*v(:)));
    else        
        cu = A\u(:);
        cv = A\v(:);    
    end
    
    % Correction constant for making the local median, 
    % mean-unbiased for normal distributed errors.
    phi = 2.1981; %=1/icdf('chi2',1/2,1);
    
    for k=2:5   
        % Calculate residual to fit
        res = reshape((u(:)-A*cu).^2 + (v(:)-A*cv).^2,M,N);        
        res(~mask) = NaN;         
        
        % 9x9 local median of residuals
        mres = phi*localmedian(res,9);         
        
        % Set weight of masked vectors to a small number 0.03
        % (This ensures a result even if all is masked),        
        % for other vectors use Tukey biweight using residuals
        % scaled with local median of the residuals.        
        W = 0.03*ones(M,N); 
        W(mask) = max(1-res(mask)./(4.685^2*mres(mask)),0).^2; % Biweight
        
        % Weighted Least squares fit
        AW = bsxfun(@times,A,W(:))';   
        [L,isposdef] = chol(AW*A,'lower');
        if(~isposdef)
            cu = L'\(L\(AW*u(:)));
            cv = L'\(L\(AW*v(:)));
        else      
            break;
        end
    end 
    
    % Return fitted least squares spline
    Ubspline = Bspline2(reshape(cu,I,J),tx,ty,4);
    Vbspline = Bspline2(reshape(cv,I,J),tx,ty,4);   
 
    % Normalize residual from last iteration for outlier detection
    res = reshape((u(:)-A*cu).^2 + (v(:)-A*cv).^2,M,N);
    res(~mask) = NaN;   
    mres = phi*localmedian(res,9);     
    globalres = sqrt(res./mres)/4.685;

    % Mark outliers
    out = (localres<1) & (globalres<1);
            
    
function Xm = localmedian(X,M)
% LOCALMEDIAN find the median of local M*M points 
% in the regular grid X, (excluding the middle point).

    N = (M-1)/2;
    [I,J] = size(X);
    idy = [N+2:M 1:I I-M+1:I-N-1];
    idx = [N+2:M 1:J J-M+1:J-N-1];
    Xext = X(idy,idx);
    
    Xi = zeros(I,J,M*M-1);
       
    k=0;
    for i=0:M-1
        for j=0:M-1
            if(i~=1 || j~=1)
                k=k+1;
                Xi(:,:,k)= Xext((1:I)+i,(1:J)+j);
            end
        end
    end
    
    Xm = nanmedian(Xi,3);