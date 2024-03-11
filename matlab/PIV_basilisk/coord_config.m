function [tform] = coord_config()
%% script for making the coordinate transforms from pixel to world
close all
coord_name1 = "MEK4600_G3_21_10_02/2021-02-10_coordinationphoto_C001H001S0001000001.bmp";
coord_name2 = "MEK4600_G3_21_17_02/2021-02-17_coordinationphoto_C001H001S0001000001.bmp";
coord_name3 = "MEK4600_G3_21_24_02/2021-02-24_coordphoto_C001H001S0001000001.bmp";

height_from_surface = [0.005 0.00029 0.0031]; %the height from the free surface to the top coordinate dots

coord = imread(coord_name1);

x0 = 73;
x1 = 1988;
y0 = 472;
y1 = 2020;
[x_pos, y_pos] = ndgrid(round(x1:(x0-x1)/26:x0), round(y1:(y0-y1)/21:y0));
pixel = [reshape(x_pos, [], 1) reshape(y_pos, [], 1)];

%refine pixel positions
c = graythresh(coord);
bw = im2bw(coord, 0.9);
cc = bwconncomp(bw);
stats = regionprops(cc,'Centroid');
xc = vertcat(stats.Centroid);
idx = knnsearch(xc,pixel);
pixel = xc(idx,:);

% Define matching reference points in world coordinate
[wx,wy] = ndgrid((13:-1:-13)*0.01,(-21:1:0)*0.01-height_from_surface(1));
%[wx,wy] = ndgrid((1:-1:-1)*0.01,(-1:1:1)*0.01-height_from_surface);
world = [wx(:) wy(:)];
[tform1, err, env] = createcoordsystem(pixel, world, 'linear');
% %show coord image and chosen poinbts
% imshow(coord)
% hold on 
% plot(pixel(:,1),pixel(:,2), 'x')

%% coordinate system from second lab day 
coord2 = imread(coord_name2);

%approximate dot positions in pixel location
x0 = 21;
x1 = 1836;
y0 = 470;
y1 = 2041;
x_points = 25;
y_points = 22;
[x_pos, y_pos] = ndgrid(round(x1:(x0-x1)/(x_points-1):x0), round(y1:(y0-y1)/(y_points-1):y0));
pixel = [reshape(x_pos, [], 1) reshape(y_pos, [], 1)];
pixel
%refine pixel positions
c = graythresh(coord2);
bw = im2bw(coord2, 0.9);
cc = bwconncomp(bw);
stats = regionprops(cc,'Centroid');
xc = vertcat(stats.Centroid);
idx = knnsearch(xc,pixel);
pixel2 = xc(idx,:);
 
% %show the initial dot guess and refined positions
% imshow(bw)
% hold on
% plot(pixel(:,1), pixel(:,2), 'rx')
% plot(pixel2(:,1), pixel2(:,2), 'yx')

% Define matching reference points in world coordinate
[wx,wy] = ndgrid((12:-1:-12)*0.01,(-21:1:0)*0.01-height_from_surface(2));
%[wx,wy] = ndgrid((1:-1:-1)*0.01,(-1:1:1)*0.01-height_from_surface);
world = [wx(:) wy(:)];
[tform2, err, env] = createcoordsystem(pixel2, world, 'linear');


% %% coordinate system from third lab day
coord3 = imread(coord_name3);
% 
%approximate dot positions in pixel location
x0 = 77;
x1 = 1996;
y0 = 466;
y1 = 2016;
x_points = 27;
y_points = 22;
[x_pos, y_pos] = ndgrid(round(x1:(x0-x1)/(x_points-1):x0), round(y1:(y0-y1)/(y_points-1):y0));
pixel = [reshape(x_pos, [], 1) reshape(y_pos, [], 1)];
% 
% %refine pixel positions
c = graythresh(coord3);
bw = im2bw(coord3, 0.9);
cc = bwconncomp(bw);
stats = regionprops(cc,'Centroid');
xc = vertcat(stats.Centroid);
idx = knnsearch(xc,pixel);
pixel3 = xc(idx,:);
%  
%show the initial dot guess and refined positions
imshow(bw)
hold on
plot(pixel(:,1), pixel(:,2), 'rx')
plot(pixel3(:,1), pixel3(:,2), 'yx')

%Define matching reference points in world coordinate
distance_between_points = 0.01;
x_positions = (floor(x_points/2):-1:-ceil(x_points/2)+1)*distance_between_points;
y_positions = (-y_points+1:1:0)*distance_between_points-height_from_surface(3);

[wx,wy] = ndgrid(x_positions, y_positions);
world = [wx(:) wy(:)];
[tform3, err, env] = createcoordsystem(pixel3, world, 'linear');


%%
tform = [tform1 tform2 tform3];
end


% %make pixel to world transform by selecting points
% if exist('tform', 'var') ==1
% 
% else
%     coord = imread(coord_name);
%     imagesc(coord)
%     h = impoly;
%     title('select from top right towards left')
%     pixel = h.getPosition;
%     
%     %refine pixel positions
%     c = graythresh(coord);
%     bw = im2bw(coord);
%     cc = bwconncomp(bw);
%     stats = regionprops(cc,'Centroid');
%     xc = vertcat(stats.Centroid);
%     idx = knnsearch(xc,pixel);
%     pixel = xc(idx,:);
%     
%     % Define matching reference points in world coordinate
%     [wx,wy] = ndgrid((10:-5:-10)*0.01,(-10:5:-5)*0.01);
%     world = [wx(:) wy(:)];
%     [tform, err, env] = createcoordsystem(pixel, world, 'linear');
% end
