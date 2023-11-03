function [im1bs,im2bs] = im2bspline(im1,im2,x,y)
  
  
  [I,J] = size(im1);
  
  if(nargin<3)
    x = 1:J;
    y = 1:I;
  end
   
  M = I - 3;
  N = J - 3;
     
  %if(nargin<5)
    ty = zeros(M+7,1);    
    ty(1:4)     = 1;
    ty(5:M+3)   = 3:I-2;    
    ty(M+4:end) = I;
    ty = y(ty);
  
    tx = zeros(N+7,1);
    tx(1:4)     = 1;
    tx(5:N+3)   = 3:J-2;    
    tx(N+4:end) = J;
    tx = x(tx);
    %end
  

%assignin('base','tx',tx)
%assignin('base','ty',ty)
  
  A = zeros(I); %sparse([],[],[],I,I,I*4);
  for i=1:I        
    A(:,i) = Bspline2.basis(i,4,y,ty); %  
  end
  %A = sparse(A);  
  [L,U] = lu(sparse(A));
  
  c1 = zeros(I,J);
  if(nargin>1)
    c2 = zeros(I,J);
  end
  
  for j=1:J
    c1(:,j) = U\(L\im1(:,j)); 
    if(nargin>1)
        c2(:,j) = U\(L\im2(:,j)); 
    end
  end
  
  B = zeros(J); %sparse([],[],[],J,J,J*4);
  for j=1:J         
    B(:,j) = Bspline2.basis(j,4,x,tx);  
  end
  [L,U] = lu(sparse(B));  
  for i=1:I
    c1(i,:) =  U\(L\c1(i,:)');    
    if(nargin>1)
        c2(i,:) =  U\(L\c2(i,:)');
    end
  end
       
  im1bs = Bspline2(c1,tx,ty,4);
  
  if(nargin>1)
    im2bs = Bspline2(c2,tx,ty,4);
  end
  %end
  
