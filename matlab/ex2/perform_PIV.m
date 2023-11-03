function [UVw, params_this_run] = perform_PIV(run_number, pair_number)
%%perform the piv passes 

%read images and constants

image_names
coord_config

im1 = imread(image_name(run_number, pair_number*2-1));
im2 = imread(image_name(run_number, pair_number*2));
mask1_name = image_name(run_number, pair_number*2-1) + ".mask.mat";
mask2_name = image_name(run_number, pair_number*2) + ".mask.mat";


height = 0.33;
f = frequency(run_number);
[omega,T,k,LAMBDA,CP,CG]=wparam(f,height);

%a = surface_height(run_number); % run the surface height script to find the amplitude
dt = 1/60;     %based on the framerate of the videos
if run_number ==1
    dt = 1/120;
end

g = 9.81;
t = 0;

window_sizes = [64 32];
window_overlap = [0.5 0.5 0.5 0.5 0.5 0.5];

%try to load the image masks 
try
    load(mask1_name)
    load(mask2_name)
catch %if it doesn't work mask and save the mask.
    mask1 = mask_image(im1);
    mask2 = mask_image(im2);
    save(mask1_name, 'mask1')
    save(mask2_name, 'mask2')
end

%% Perform PIV passes
piv1 = [];
for i= 1:size(window_sizes, 2)
    search_range = [-round(window_sizes(i)/3) round(window_sizes(i)/3) -round(window_sizes(i)/3) round(window_sizes(i)/3)];

    opt = setpivopt('range',search_range,'subwindow',window_sizes(i),window_sizes(i),window_overlap(i));
    piv1 = normalpass(piv1,im1,mask1,im2,mask2,opt);
    piv2 = piv1;
end
%%
%Replace outliers
[U,V, x, y] = replaceoutliers(piv1);

%Convert to world coordinates
[Uw, Vw, xw, yw] = pixel2world(tform(ceil(run_number/5)), U, V, x, y, dt);

idx = interp2(double(mask1&mask2),piv1.x,piv1.y)==1;





%% Save the results
try
    load('velocities.mat');
catch
    velocities = containers.Map('KeyType', 'double', 'ValueType', 'any');
end

try
    run_velocities = velocities(run_number);
catch
    run_velocities = containers.Map('KeyType', 'double', 'ValueType', 'any');
end

UVw(1,:,:) = Uw;
UVw(2,:,:) = Vw;
UVw(3,:,:) = xw;
UVw(4,:,:) = yw;
UVw(5,:,:) = idx;
run_velocities(pair_number) = UVw;
velocities(run_number) = run_velocities;
save('velocities.mat', 'velocities');

try 
    load('params.mat')
catch
    params = containers.Map('KeyType', 'double', 'ValueType', 'any');
end
params_this_run = containers.Map;
params_this_run('a') = surface_height(run_number);
params_this_run('k') = k;
params_this_run('omega') = omega;

params(run_number) = params_this_run;

save('params.mat', 'params')






end