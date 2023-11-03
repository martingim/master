function [tform,err,errinv] = createcoordsystem(pixel,world,transformtype,name)
%CREATECOORDSYSTEM creates a coordinate system transfromation.
%  
% Usage:
%
% [tform,err,errinv] = createcoordsystem(pixel,world,transformtype,name)  
%  
% input  - pixel         - Points in image coordinates [pixel]
%        - world         - Points in world coordinates [units]
%        - transformtype - linear, quadratic, cubic,
%                          bilinear, biquadratic, bicubic,
%                          linear-quadratic, quadratic linear,
%                          linear-cubic or cubic-linear.
%        - name          - filename without extension, 
%                          output is save to name.mat.
% output - tform         - Coordinate transformation for use with
%                          pixel2world,  tformfwd or tforminv.
%        - err           - rms err forward transformation (pixel->world)
%        - errinv        - rms err inverse transformation (world->pixel)
%
% % Example
%
% % Load coordinate image
% coord = imread('wavecoord.png');
% imagesc(coord)
%
% % Select reference points in pixel coordinate.
% % It is often convinient to select the pixel control point from 
% % right to left, and then upwards matching the control points made 
% % with ndgrid.
% h=impoly;
% pixel = h.getPosition;
%
% % Refine pixel positions (optional)
% c = graythresh(coord);
% im2bw = im2bw(coord);
% cc = bwconncomp(bw);
% stats = regionprops(cc,'Centroid');
% xc = vertcat(stats.Centroid);
% idx = knnsearch(xc,pixel);
% pixel = xc(idx,:);
%  
% % Define matching reference points in world coordinate
% [wx,wy] = ndgrid((9:-1:0)*0.02 + 4.33,(-10:0)*0.02 + 0.103);
% world = [wx(:) wy(:)];
%   
% % Create coordinate transformation
% [tform,err,errinv] = createcoordsystem(pixel,world,'cubic')   
%
% See also pixel2world, mapping, tformfwdm, tforminv, maketform

  if (nargin<3)
    error('Need three input arguments.')
  end
    
  tdata.transformtype = transformtype;
    
  % Forward transformation
  X = mapping(pixel(:,1),pixel(:,2),transformtype);
  [Q,R] = qr(X);
  tdata.T = R\(Q'*world);
  err = sqrt(mean(sum((mapping(pixel(:,1),pixel(:,2),transformtype)*tdata.T - world).^2,2)));
  
  % Backward transformation
  U = mapping(world(:,1),world(:,2),transformtype);
  [Q,R] = qr(U);
  tdata.Tinv = R\(Q'*pixel); 
  errinv =  sqrt(mean(sum((mapping(world(:,1),world(:,2),transformtype)*tdata.Tinv - pixel).^2,2)));   
  
  % Create tform
  tform = maketform('custom',2,2,@fwd_custom,@inv_custom,tdata);

  % Write to mat-file  
  if (nargin>=4)
    save(name,'-v7.3','pixel','world','tform','err','errinv');
  end 
     
function world = fwd_custom(pixel,T)
  world = mapping(pixel(:,1),pixel(:,2),T.tdata.transformtype)*T.tdata.T;
  
function pixel = inv_custom(world,T)
  pixel = mapping(world(:,1),world(:,2),T.tdata.transformtype)*T.tdata.Tinv;  
  
