%Setup
current_path = cd;
cd '/home/martin/Documents/master/matlab/HydrolabPIV';
setup_hlpiv
cd(current_path)
image_names
coord_config
h = 0.33; %heightof the water surface at rest

close all
%% Perform PIV
n_runs = 15;
n_waves = 1;
max_arrows = 50; %max number of arrows in one dim in the quiver plots
%height = 0.33;

%%Perform PIV
for run_number=11:n_runs
    for wave=1:n_waves
        perform_PIV(run_number, wave);
        disp('PIV performed');
    end
end

%% Plot velocities
run_number =8;
pair_number = 1;

im1 = imread(image_name(run_number, pair_number*2-1));
im2 = imread(image_name(run_number, pair_number*2));
mask1_name = image_name(run_number, pair_number*2-1) + ".mask.mat";
mask2_name = image_name(run_number, pair_number*2) + ".mask.mat";

load(mask1_name);
load(mask2_name);
%surface 
[idx1,eta1] = max(mask1);
[idx2,eta2] = max(mask2);

[etax,etay] = tformfwd(tform(ceil(run_number/5)),1:size(im1,2),(eta1+eta2)/2);

load('velocities.mat')
load('params.mat')
run_number = 8;
alpha = containers.Map('KeyType', 'double', 'ValueType', 'any');
crest_mask = containers.Map('KeyType', 'double', 'ValueType', 'any');
y_scaled = containers.Map('KeyType', 'double', 'ValueType', 'any');
for wave = 1:3
    [y_scaled(wave), u_crest_scaled, alpha(wave), crest_mask(wave)] = plot_velocities(run_number, wave, 50);
    %close all
end

mask = crest_mask(1)&crest_mask(2)&crest_mask(3);

y1 = y_scaled(1);
y2 = y_scaled(2);
y3 = y_scaled(3);
a1 = alpha(1);
a2 = alpha(2);
a3 = alpha(3);

[tmp, i] = min([size(y1, 1) size(y2, 1) size(y3, 1)]);
y = y_scaled(i);

alpha_mean = (a1 + a2 + a3)/3;

plot(alpha_mean(mask), y(:,1), 'x')
hold on


plot(a1(crest_mask(1)), y1(:,1), 'x')
plot(a2(crest_mask(2)), y2(:,1), 'x')
plot(a3(crest_mask(3)), y3(:,1), 'x')
p = params(run_number);
k = p('k');
a = p('a');
plot(exp(k*y1(:,1)*h), y1(:,1))
legend('mean', '1', '2', '3', 'analytical')


