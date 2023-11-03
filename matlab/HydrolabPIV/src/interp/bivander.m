function P = bivander(X,Y,M,N)
% BIVANDER creates a bivariate Vandermonde matrix
%
% For a given set of data points (X,Y) use P = BIVANDER(X,Y,M,N)
% to create a bivariate Vandermonde matrix where (M,N) is the polynomial order.
% 
% Example
% 
% % Test function
% [x,y,z] = peaks();
% 
% % Ninth order multivariate least squares fit
% X = [x(:) y(:)];
% P = multivander(X,9);
% c = P\z(:);
% zi = reshape(P*c,size(z));
% 
% figure
% surf(x,y,zi)
%
% See also VANDER, MULTIVANDER

if(nargin<4)
    error('Needs four input arguments: P = bivander(X,Y,M,N)');    
end

M=M+1;
N=N+1;

X=X(:); I = size(X,1);
Y=Y(:); J = size(Y,1);

if (I~=J)
  error('X,Y must be of same size!');
end

% Single variate vandermonde matrix for X
A = ones(I,M,'like',X);
for m=M-1:-1:1
   A(:,m) = A(:,m+1).*X;
end

% Single variate vandermonde matrix for Y
B = ones(J,N,'like',Y);
for n=N-1:-1:1
   B(:,n) = B(:,n+1).*Y;
end
 
% Bivariate vandermonde matrix
P = zeros(I,M*N);
for n=1:N
  for m=1:M    
    i = m+M*(n-1);      
    P(:,i) = A(:,m).*B(:,n);  
  end
end
