%Setup
current_path = cd;
cd '/home/martin/Documents/master/matlab/HydrolabPIV';
setup_hlpiv
cd(current_path)

%parameters from the article
a = 0.0205;     %measured amplitude of the wave
omega = 8.95;   %frequency of the waves
h = 0.6;        %the water level in the tank
water_depth = 0.6; %the water level in the tank at rest
frequency = [omega/(2*pi)];
time_between_frames = [1/60]; %%%%%%change just a guess

%the paths to the images on the same format as my previous lab excerise

image_name = string(zeros(1,2));
image_name(1,1) = 'Jensen_images_2001/mpim1b.jpg';
image_name(1,2) = 'Jensen_images_2001/mpim1c.jpg';
coord_name = 'Jensen_images_2001/';% add this image later
%% create coord config for the lab

coord = imread(coord_name);
top_point_distance_from_surface = 0;

%max and min pixel indexes of the points on the coord photo
x0 = 73;
x1 = 1988;
y0 = 472;
y1 = 2020;
%number of points in each direction on the coord photo
n_points_x = 27;
n_points_y = 22;
%used to guess the dot positions, later refined with knnsearch
[x_pos, y_pos] = ndgrid(round(x1:(x0-x1)/(n_points_x-1):x0), round(y1:(y0-y1)/(n_points_y-1):y0));
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
distance_between_points = 0.01;
x_positions = (floor(n_points_x/2):-1:-ceil(n_points_x/2)+1)*distance_between_points;
y_positions = (-n_points_y+1:1:0)*distance_between_points-top_point_distance_from_surface;
[wx,wy] = ndgrid(x_positions,y_positions);
world = [wx(:) wy(:)];
[tform1, err, env] = createcoordsystem(pixel, world, 'cubic');

%% perform PIV
perform_PIV(1,1)