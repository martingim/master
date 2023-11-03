function [r,masked,out] = localfilt(X,Y,U,V,mask,tol,rmax,M)
% LOCALFILT performs normalized local median outlier detection
% on regular grid with mask. For more information, see 
% <a href="matlab: 
% web('http://link.springer.com/article/10.1007/s00348-005-0016-6')"
% >Universal outlier detection for PIV data (Westerweel & Scarano, 2005)</a>.
%
% [r,masked,out] = localfilt(X,Y,U,V,mask,tol,rmax,M)
% 
% input  - (X,Y)  - Positions
%        - (U,V)  - Velocities
%        - mask   - Masked vectors
%        - tol    - Tolerance, default value is 0.1
%        - rmax   - Residual threshhold value, default is 2
%        - M      - Filter size, usually 3 (3x3 neighbour grid)
% output - r      - Residuals 
%        - masked - Number of masked neighbours
%        - out    - Outliers (|r|<rmax)
% 
% See also replaceoutliers, lsbsfit2

% Check and set default input arguments
if(nargin<8)
  M = 3;
  if(nargin<7)
    rmax = 2;
    if(nargin<6)
      tol = .1;
    end
  end
else
  if(mod(M,2)==0) 
    M=M+1;
  end
end
N = (M-1)/2;

% Apply mask
U(~mask) = NaN;
V(~mask) = NaN;

% Pad velocity field
[I,J] = size(X);
idy = [N+2:M 1:I I-M+1:I-N-1];
idx = [N+2:M 1:J J-M+1:J-N-1];
Xext = X(idy,idx);
Yext = Y(idy,idx);
Uext = U(idy,idx);
Vext = V(idy,idx);
Mext = mask(idy,idx);

% Find the M*M-1 nearest neighbors (including masked ones)
k=0;
Xi = zeros(I,J,M*M-1);
Yi = zeros(I,J,M*M-1);
Ui = zeros(I,J,M*M-1);
Vi = zeros(I,J,M*M-1);
Mi = zeros(I,J,M*M-1);
for i=0:M-1
  for j=0:M-1
    if(i~=1 || j~=1)
      k=k+1;
      Xi(:,:,k)= Xext((1:I)+i,(1:J)+j);
      Yi(:,:,k)= Yext((1:I)+i,(1:J)+j);
      Ui(:,:,k)= Uext((1:I)+i,(1:J)+j);
      Vi(:,:,k)= Vext((1:I)+i,(1:J)+j);
      Mi(:,:,k)= Mext((1:I)+i,(1:J)+j);
    end
  end
end

% Calculate Euclidian distance
dist = sqrt(bsxfun(@minus,Xi,X).^2 + bsxfun(@minus,Yi,Y).^2);

% Normalized tolerance with median distance
dm = median(dist,3); % mean??
q = -(dm + sign(dm).*sqrt(dm.^2 + 4*tol))/2;
ntol = -tol./q; % max(q,-tol./q); % non-negative root  

% Normalize velocity vectors
U0 = U./(dm+ntol);
V0 = V./(dm+ntol);
Ui = Ui./bsxfun(@plus,dist,ntol);
Vi = Vi./bsxfun(@plus,dist,ntol);

% Calculate residual from normalized velocity vectors
r = zeros(I,J,2);
Um = nanmedian(Ui,3);
Vm = nanmedian(Vi,3);
r(:,:,1) = abs(U0-Um)./(nanmedian(abs(bsxfun(@minus,Ui,Um)),3)+ntol);
r(:,:,2) = abs(V0-Vm)./(nanmedian(abs(bsxfun(@minus,Vi,Vm)),3)+ntol);

% Number of masked vectors
masked = sum(Mi,3);

% Threshold the residuals to find outliers
out = max(r,[],3)< rmax;

