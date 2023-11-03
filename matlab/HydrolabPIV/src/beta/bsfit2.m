function [c,tx,ty,A,B] = bsfit2(x,y,f,tx,ty)
  
  [I,J] = size(f);
   
  %I = length(ty);
  %J = length(tx);
  M = I - 3;
  N = J - 3;
     
  %if(nargin<5)
    ty = zeros(M+7,1);    
    ty(1:4)     = y(1,1);
    ty(5:M+3)   = y(3:I-2,1);    
    ty(M+4:end) = y(I,1);
    
    tx = zeros(N+7,1);
    tx(1:4)     = x(1,1);
    tx(5:N+3)   = x(1,3:J-2);    
    tx(N+4:end) = x(1,J);
  %end

  
  A = sparse([],[],[],I,I,I*4);
  for i=1:I        
    A(:,i) = basis(i,4,y(:,1),ty);    %Bspline2.basis(i,4,y(:,1),ty); %  
  end
  [L,U] = lu(A);
  for j=1:J
    c(:,j) = U\(L\f(:,j)); 
  end
  
  B = sparse([],[],[],J,J,J*4);
  for j=1:J         
    B(:,j) = basis(j,4,x(1,:)',tx); %Bspline2.basis(j,4,x(1,:)',tx);  
  end
  [L,U] = lu(B);  
  for i=1:I
    c(i,:) =  U\(L\c(i,:)');
  end
   
  %if(nargout==1)
  %  c = Bspline2(c,tx,ty,4);
  %end
  
function N = basis(i,k,x,t)
  M = size(t,1);         
  
  if (k==1)   
              
    N = (t(i) <= x & x < t(i+1));
    
  else
    a = (t(i+k-1)-t(i));
    if (a~=0)
      N = (x-t(i))/a.*basis(i,k-1,x,t);
    else
      N = 0;
    end    
    b = (t(i+k)-t(i+1));
    if (b~=0)
      N = N + (t(i+k)-x)/b.*basis(i+1,k-1,x,t);  
    end
  end
  
    
  % Fix endpoint ?
  if (x(end)==t(end) & i==M-k)    
    N(end) = 1;    
  end