function [xa,ya,xb,yb] = rungekutta_less_memory(tableau,piv,x0,y0)
%[xa,ya,xb,yb] = rungekutta(tableau,piv,x0,y0)
%
%
%
% Implicit Runge-Kutta given by the Butcher-Tableau
% ____________________
%  0 | 
% c2 | a21
% c3 | a31 a32 
% ...| ... ... ... 
% ci | ai1 ai2 ... aij
% ___|________________
%    | b1  b2  ... bj
%
% Note 1: Any explicit coeffisients are ignored.
%
% Note 2: The c-coeffisients are not used as the velocity field does not
%         depend on t.


    [M,N] = size(x0);
    I = size(tableau,1)-1;

    kx = zeros(M,I);
    ky = zeros(M,I);

    % Half time-step
    h = piv.opt{end}.alpha/2;
    tableau = h*tableau;

    % Integrate backward
    xa = x0; % zeros(M,N);
    ya = y0; %zeros(M,N);
    xb = x0; % zeros(M,N);
    yb = y0; %zeros(M,N);
    for n=1:N
        for i=1:I
            %dt = t0 + h*tableau(i,1);
            xi = x0(:,n);
            yi = y0(:,n);   
            for j=1:i-1
                xi = xi - tableau(i,j+1)*kx(:,j);
                yi = yi - tableau(i,j+1)*ky(:,j);
            end
    
            kx(:,i) = piv.Uspline.evaluate(xi,yi);
            ky(:,i) = piv.Vspline.evaluate(xi,yi);    
    
            xa(:,n) = xa(:,n) - tableau(I+1,i+1)*kx(:,i);
            ya(:,n) = ya(:,n) - tableau(I+1,i+1)*ky(:,i);
        end

        % Integrate forward
        xb(:,n) = xb(:,n) + tableau(I+1,2)*kx(:,1);
        yb(:,n) = yb(:,n) + tableau(I+1,2)*ky(:,1);
        for i=2:I
            %dt = t0 + h*tableau(i,1);
            xi = x0(:,n); 
            yi = y0(:,n);   
            for j=1:i-1
                xi = xi + tableau(i,j+1)*kx(:,j);
                yi = yi + tableau(i,j+1)*ky(:,j);
            end
    
            kx(:,i) = piv.Uspline.evaluate(xi,yi);
            ky(:,i) = piv.Vspline.evaluate(xi,yi);    
    
            xb(:,n) = xb(:,n) + tableau(I+1,i+1)*kx(:,i);
            yb(:,n) = yb(:,n) + tableau(I+1,i+1)*ky(:,i);    
        end
    end
    %assignin('base','kx',kx)
    %assignin('base','ky',ky)
end

