function [X,dXdx,dXdy] = mapping(x,y,transformtype)
% MAPPING defines a mapping function from on coordinate system to another,
% usually pixel to world, for use in the tform created by createcoordsystem.
%
% [X,dXdx,dXdy] = mapping(x,y,transformtype)
% 
% Returns X = polynomial(x,y) where the transformtype can be 'linear',
% 'quadratic','cubic','bilinear', 'biquadratic', 'bicubic', 'linear-quadratic',
% 'quadratic-linear','linear-cubic' or 'cubic-linear'.
%
% dXdx and dXdy (Jacobian) can be used to transform vectors.
% 
% See also createcoordsystem, pixel2world

  M = size(x,1); 
  I = ones(M,1);
  O = zeros(M,1);
    
  switch (transformtype)
   case 'linear' 
    X =    [x y I];
    dXdx = [I O O];
    dXdy = [O I O];        
   case 'quadratic' 
    X =    [x.^2 x.*y y.^2 x y I];
    dXdx = [2*x  y    O    I O O];
    dXdy = [O    x    2*y  O I O];  
   case 'cubic' 
    X =    [x.^3   x.^2.*y x.*y.^2 y.^3   x.^2 x.*y y.^2 x y I];
    dXdx = [3*x.^2 2*x.*y  y.^2    O      2*x  y    O    I O O];
    dXdy = [O      x.^2    2.*x.*y 3*y.^2 O    x    2*y  O I O];  
   case 'bilinear' 
    X =    [x.*y x y I];   
    dXdx = [y    I O O];
    dXdy = [x    O I O];    
   case 'linear-quadratic' 
    X =    [x.*y.^2 x.*y y.^2 x y I];
    dXdx = [y.^2    y    O    I O O];
    dXdy = [2*x.*y  x    2*y  O I O];       
   case 'linear-cubic'           
    X =    [x.*y.^3   x.*y.^2 y.^3   x.*y y.^2 x y I];
    dXdx = [y.^3      y.^2    O      y    O    I O O];
    dXdy = [3*x.*y.^2 2.*x.*y 3*y.^2 x    2*y  O I O];  
   case 'quadratic-linear' 
    X =    [x.^2.*y x.^2 x.*y x y I];
    dXdx = [2*x.*y  2*x  y    I O O];
    dXdy = [x.^2    O    x    O I O];  
   case 'biquadratic' 
    X =    [x.^2.*y.^2 x.^2.*y x.*y.^2 x.^2 x.*y y.^2 x y I];
    dXdx = [2*x.*y.^2  2*x.*y  y.^2    2*x  y    O    I O O];
    dXdy = [2*x.^2.*y  x.^2    2*x.*y  O    x    2*y  O I O];  
   case 'cubic-linear'     
    X =    [x.^3.*y   x.^2.*y x.^3   x.^2 x.*y x y I];    
    dXdx = [3*x.^2.*y 2*x.*y  3*x.^2 2*x  y    I O O];    
    dXdy = [x.^3      x.^2    O      O    x    O I O];  
   case 'bicubic' 
    X = [x.^3.*y.^3 x.^3.*y.^2 x.^2.*y.^3  ...
         x.^3.*y x.^2.*y.^2 x.*y.^3 ...     
         x.^3 x.^2.*y x.*y.^2 y.^3 ... 
         x.^2 x.*y y.^2 x y I];
    dXdx = [3.*x.^2.*y.^3 3*x.^2.*y.^2 2*x.*y.^3  ...
            3*x.^2.*y 2*x.*y.^2 y.^3 ...
            3*x.^2 2*x.*y y.^2 O ...
            2*x y O I O O];
    dXdy = [3.*x.^3.*y.^2 2*x.^3.*y 3*x.^2.*y.^2  ...
            x.^3 2*x.^2.*y 3*x.*y.^2 ...
            O x.^2 2*x.*y 3*y.^2 ...
            O x 2*y O I O];  
   otherwise
    error(['Transformtype "%s" is not supported.\n' ...
           'Must be linear, quadratic, cubic, ' ...
           'bilinear, biquadratic, bicubic,'
           'linear-quadratic, quadratic-linear, ' ...
           'liner-cubic or cubic-linear.'],transformtype);   
  end
  
  
  