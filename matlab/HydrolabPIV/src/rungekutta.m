function [x1,y1,x2,y2] = rungekutta(tableau,piv,x0,y0)
% RUNGEKUTTA calculates the distortion of image 1 (x1,y1) and image 2 (x2,y2)
% from the initial position (x0,y0).
%
% [x1,y1,x2,y2] = rungekutta(tableau,piv,x0,y0)
%
% The distortion is found by integrating the velocity field given by 
% the piv pass using implicit Runge-Kutta given by the Butcher-Tableau
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
% 
% See also DISTORTEDPASS


    [M,N] = size(x0);
    I = size(tableau,1)-1;

    kx = zeros(M,N,I);
    ky = zeros(M,N,I);

    % Half time-step
    h = piv.opt{end}.alpha/2;
    tableau = h*tableau;

    % Integrate backward
    x1 = x0; % zeros(M,N);
    y1 = y0; %zeros(M,N);
    for i=1:I
        %dt = t0 + h*tableau(i,1);
        xi = x0; 
        yi = y0;   
        for j=1:i-1
            xi = xi - tableau(i,j+1)*kx(:,:,j);
            yi = yi - tableau(i,j+1)*ky(:,:,j);
        end
    
        kx(:,:,i) = piv.Uspline.evaluate(xi,yi);
        ky(:,:,i) = piv.Vspline.evaluate(xi,yi);    
    
        x1 = x1 - tableau(I+1,i+1)*kx(:,:,i);
        y1 = y1 - tableau(I+1,i+1)*ky(:,:,i);
    end

    % Integrate forward
    x2 = x0 + tableau(I+1,2)*kx(:,:,1);
    y2 = y0 + tableau(I+1,2)*ky(:,:,1);
    for i=2:I
        %dt = t0 + h*tableau(i,1);
        xi = x0; 
        yi = y0;   
        for j=1:i-1
            xi = xi + tableau(i,j+1)*kx(:,:,j);
            yi = yi + tableau(i,j+1)*ky(:,:,j);
        end
    
        kx(:,:,i) = piv.Uspline.evaluate(xi,yi);
        ky(:,:,i) = piv.Vspline.evaluate(xi,yi);    
    
        x2 = x2 + tableau(I+1,i+1)*kx(:,:,i);
        y2 = y2 + tableau(I+1,i+1)*ky(:,:,i);    
    end
end

