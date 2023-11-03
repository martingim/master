function [uw,vw,xw,yw] = pixel2world(tform,up,vp,xp,yp,dt)
% PIXEL2WORLD transform vectror from pixel coordinate to world coordinates, 
% using tform created by createcoordsystem or other tforms.
%
% [u,v,x0,y0] = pixel2world(tform,dx,dy,x,y,dt)
% 
% input  - tform   - coordinate transformation
%        - (up,vp) - pixel displacement
%        - (xp,yp)   - pixel position
%        - dt      - time step
% output - (uw,vw) - velocity
%        - (xw,yw) - position
%
% See also createcoordsystem, mapping, maketform

[I,J] = size(xp);

if(isfield(tform,'tdata') && isfield(tform.tdata,'transformtype'))   
    % Use exact Jacobian
    xp = double(xp(:));
    yp = double(yp(:));
    up = double(up(:));
    vp = double(vp(:));

    [X,dXdx,dXdy] = mapping(xp,yp,tform.tdata.transformtype);

    Y = X*tform.tdata.T;
    xw = reshape(Y(:,1),I,J);
    yw = reshape(Y(:,2),I,J);
    
    dYdx = dXdx*tform.tdata.T;
    dYdy = dXdy*tform.tdata.T;
    
    uw=reshape(up.*dYdx(:,1) + vp.*dYdy(:,1),I,J)/dt;
    vw=reshape(up.*dYdx(:,2) + vp.*dYdy(:,2),I,J)/dt;    
else
  % Use numerical approximation 
  [xw,yw]=tformfwd(tform,double(xp),double(yp));  
  [x1,y1]=tformfwd(tform,double(xp-.5*up),double(yp-.5*vp));  
  [x2,y2]=tformfwd(tform,double(xp+.5*up),double(yp+.5*vp));
     
  uw = (x2-x1)/dt;
  vw = (y2-y1)/dt;
end


