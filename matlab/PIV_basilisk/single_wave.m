%Setup
current_path = cd;
cd '/home/martin/Documents/master/matlab/HydrolabPIV';
setup_hlpiv
cd(current_path)

run_number = 16;
%parameters from the article
a = 0.0205;     %measured amplitude of the wave
omega = 8.95;   %frequency of the waves
k = 7.95;       %wavenumber
h = 0.6;        %the water level in the tank
water_depth = zeros(16,1);
frequency = zeros(16,1);
time_between_frames = zeros(16,1);
tform = [];
water_depth(run_number) = 0.6; %the water level in the tank at rest
frequency(run_number) = omega/(2*pi);
time_between_frames(run_number) = 0.012;

%the paths to the images on the same format as my previous lab excerise

image_name = string(zeros(16,2));
image_name(run_number,1) = 'Jensen_images_2001/mpim1b.bmp';
image_name(run_number,2) = 'Jensen_images_2001/mpim1c.bmp';
coord_name = 'Jensen_images_2001/mpwoco.bmp';
%% create coord config for the lab

coord = imread(coord_name);
top_point_distance_from_surface = 0.0495;

%max and min pixel indexes of the points on the coord photo
x0 = 88;
x1 = 1264;
y0 = 524;
y1 = 820;
%s = 235
%number of points in each direction on the coord photo
n_points_x = 5;
n_points_y = 2;
%used to guess the dot positions, later refined with knnsearch
[x_pos, y_pos] = ndgrid(round(x1:(x0-x1)/(n_points_x-1):x0), round(y1:(y0-y1)/(n_points_y-1):y0));
pixel_guess = [reshape(x_pos, [], 1) reshape(y_pos, [], 1)];
%refine pixel positions
c = graythresh(coord);
bw = im2bw(coord, 0.95);
cc = bwconncomp(bw);
stats = regionprops(cc,'Centroid');
xc = vertcat(stats.Centroid);
idx = knnsearch(xc,pixel_guess);
pixel = xc(idx,:);

% Define matching reference points in world coordinate
distance_between_points = 0.05;
x_positions = (floor(n_points_x/2):-1:-ceil(n_points_x/2)+1)*distance_between_points;
y_positions = (-n_points_y+1:1:0)*distance_between_points-top_point_distance_from_surface;
[wx,wy] = ndgrid(x_positions,y_positions);
world = [wx(:) wy(:)]; 
[tform1, err, env] = createcoordsystem(pixel, world, 'linear');
tform = [tform1 tform1 tform1 tform1];

% %show the initial dot guess and refined positions
% imshow(coord)
% hold on
% plot(pixel_guess(:,1), pixel_guess(:,2), 'rx')
% plot(pixel(:,1), pixel(:,2), 'yx')


%save all the parameters to image_params container
image_params = containers.Map;
image_params('image_name') = image_name;
image_params('frequency') = frequency;
image_params('time_between_frames') = time_between_frames;
image_params('water_depth') = water_depth;
image_params('tform') = tform; 

%% perform PIV
perform_PIV(run_number,1,image_params);

%% save some constants to params
load params.mat params
p = params(run_number);
p('a') = a;
p('std_a') = 0.0005;
p('water_depth') = water_depth;
params(run_number) = p;
save('params.mat', 'params')

%% plot velocity under crest
% timestep = 1;
close all
plot_velocity_under_crest(run_number,1,image_params);
print('~/Documents/master/movies_and_figures/PIV_horizontal_velocity_under_crest', '-dpng')
%%
%compare with basilisk results
% MOVING PISTON
timestep=0; 
basilisk_moving_piston_velocity_profile(timestep, 0.0244, omega, h)
timestep=1;
basilisk_moving_piston_velocity_profile(timestep, 0.02, omega, h)
% plot_basilisk_velocity_profile(timestep, a, omega, h)
% timestep = 19;
% plot_basilisk_velocity_profile(timestep, a, omega, h)

% plot_basilisk_multilayer_velocity_profile(256,20,a,omega,h)

%% plot alpha

plot_alpha(run_number, 1, image_params, true);
timestep=0;
basilisk_moving_piston_alpha(timestep, 0.0244, k, omega, h)

timestep=1;
basilisk_moving_piston_alpha(timestep, 0.025, k, omega, h)

%compare with basilisk results
% timestep = 0;
% plot_basilisk_alpha(timestep, a, k, omega, h);
% timestep = 10;
% plot_basilisk_alpha(timestep, a, k, omega, h);
% timestep = 20;
% plot_basilisk_alpha(timestep, a, k, omega, h);
% basilisk_multilayer_alpha(a, k, omega, h, 512, 20);
% basilisk_multilayer_alpha(a, k, omega, h, 512, 21);
% basilisk_multilayer_alpha(a, k, omega, h, 512, 40);

%% plot energy
figure;
nx = [128 256 512];
nl = [10 20 40];
for i=nx
    for j=nl
        plot_energy_multilayer(i,j)
    end
end
% plot_energy_multilayer(512, 40)
%% 
figure;
plot_energy_multilayer(512, 20) 
plot_energy_ns(32, 256)
plot_energy_ns(32, 512)
plot_energy_ns(64, 256)
