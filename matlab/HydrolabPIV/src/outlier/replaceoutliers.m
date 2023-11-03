function [U,V,x,y] = replaceoutliers(piv,mask)
% REPLACEOUTLIERS replaces the outliers marked in piv.out by 
% evaluating the weighted least squares bsplines (piv.Uspline,piv.Vspline)
%
% [U,V,x,y] = REPLACEOUTLIERS(piv,mask)
% 
% input  - piv     - PIV structure returned by normalpass, distortedpass.
%        - mask    - Use image mask, or specify a minimum unmasked fraction 
%                    (scalar value) to set masked vectors to nan.
% output - (U,V)   - Velocity field with replaced outliers.
%        - (x,y)   - Positions of vectors.
%
% See also normalpass, distortedpass, lsbsfit2, localfilt

U = piv.U;
V = piv.V;

idx =~piv.out;
if(sum(idx(:))>0)
    x = piv.x(idx);
    y = piv.y(idx);
    U(idx) = piv.Uspline.evaluate(x,y);
    V(idx) = piv.Vspline.evaluate(x,y);
end

x = piv.x;
y = piv.y;

if(nargin>1)    
    if(isscalar(mask))
        % Set vector with minimum unmasked fraction less than mask to NaN
        idm = (piv.masked<mask);
        U(idm) = NaN;
        V(idm) = NaN;    
    else
        % Set masked vector to NaN using image mask
        % Pad array to get correct behavior at edges.
        I=32;
        mask = padarray(double(mask),[I I],'symmetric'); 
        idm = interp2(mask,x+I,y+I,'nearest')==0;  
        U(idm) = NaN;
        V(idm) = NaN;
    end
end
