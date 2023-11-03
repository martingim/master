function [x0,delta,out] = subpixelLM(F)
  
 [I,J] = size(F); 
 
 
  [x,y]=meshgrid(1:J,1:I);
  %P = bivander(x,y,2,2);
  x=x(:);
  y=y(:);
  
  [M,N] = size(F);
      
  [~,idx] = sort(F(:),'descend');        
  
  [m,n] = ind2sub([M,N],idx(1));       
        
    f = F; %(F(m-1:m+1,n-1:n+1));      
    x0 = [n; m]; 
   
    % Intial guess
    c = [0 .5 8/9 0 8/9 n m]';
    lambda = .5;
    for i=1:10
        [fi,J] = gfit(c,x(:),y(:));    
        %size(fi)
        %size(J)
        dc = (J*J'+lambda*eye(7))\(J*(f(:)-fi));
        c = c + dc;
    end
    assignin('base','c',c);
    assignin('base','J',J);
    assignin('base','fi',fi);
    assignin('base','f',f);
    
    delta = -c(6:7);
    if (max(abs(delta))<= 1.0)
        % Peak found
        %x0 = [n - (N+1)/2; m-(M+1)/2];        
        x0 = [n; m];  
        out = 1;
        return;
    end 
    
    
    
    %return    
  
  
  % No peak found!
  out = 0;
  delta = [0;0];          
  [m,n] = ind2sub([M,N],idx(1));          
  x0 = [n; m]; 
  
function [f,J] = gfit(c,x,y)
  xx = (x-c(6)).*(x-c(6)); 
  xy = (x-c(6)).*(y-c(7));
  yy = (y-c(7)).*(y-c(7));
  p = c(3)*xx + 2*c(4)*xy + c(5)*yy;        
  f = c(2)*exp(-p);
                
  J = ones(7,numel(x));
  J(2,:) = f/c(2);
  J(3,:) = -xx.*f;
  J(4,:) = -2*xy.*f;
  J(5,:) = -yy.*f;
  J(6,:) = 2*(c(3)*(x-c(6)) + c(4)*(y-c(7))).*f;
  J(7,:) = 2*(c(4)*(x-c(6)) + c(5)*(y-c(7))).*f;
     
    