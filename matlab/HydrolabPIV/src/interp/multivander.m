function [P,E,C] = multivander(X,N)
% MULTIVANDER creates multivariate Vandermonde matrix
%
% For a given set of data point X of size MxD, 
% where M is the number of points with dimension D,
% [P,E,diffconst] = MULTIVANDER(X,N)
% creates a multivariate Vandermonde matrix of polynomial order N.
%
% E gives the exponents for polynomial for each variable.
% C gives a constant, useful for evaluateing the derivatives of 
% the polynomial in (0,0,...)
% 
% Example
% 
% % Test function
% [x,y,z] = peaks();
% 
% % Ninth order bivariate least squares fit
% P = bivander(x(:),y(:),9,9);
% c = P\z(:);
% zi = reshape(P*c,size(z));
% 
% figure
% surf(x,y,zi)
% 
% See also VANDER, BIVANDER

  if(nargin<2)
    error('Needs two input arguments: P = multivander(X,N)');    
  end

  [M,D] =size(X);
  
  % Find polygonal numbers
  K = max(D,N)+1;
  T = ones(K);
  for i=2:K
    for j=i:K
      T(i,j) = T(i-1,j)+T(i,j-1);
      T(j,i) = T(i,j);
    end
  end
  L=T(N+1,D+1);
  
  % Assemble Vandermonde matrix
  %fprintf('M=%d, L=%d\n',M,T(N+1,D+1));
  P = ones(M,L);  
  if (nargout>1)
    E = zeros(D,L);
    I = eye(D);
  end
  
  
  i=1;
  j=1;
  for n=1:N
    i=i+T(n,D);
    for d=1:D       
      J = T(n,D-d+1);
             
      idi = L+1-(i-(J:-1:1));
      idj = L+1-(j+(1:J));
            
      P(:,idj) = bsxfun(@times,X(:,d),P(:,idi));
      if (nargout>1)
        E(:,idj) = bsxfun(@plus,E(:,idi),I(:,d));
      end
      
      j=j+J;
    end    
  end

  if(nargout>2) 
    fact = ones(N+1,1);
    
    for n=3:N+1;
      fact(n)=fact(n-1)*(n-1);
    end
    
    C = ones(1,L);
    for d=1:D      
      C = C./fact(E(d,:)+1)';
    end
  end