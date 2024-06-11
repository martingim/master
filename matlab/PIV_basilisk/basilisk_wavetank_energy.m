
%basilisk_folder = "~/piston-moving";
basilisk_folder = "~/Documents/master/basilisk/piston-moving";
%basilisk_folder = "~/2piston";
timestep = 25;

omega = 8.95;
h = 0.6;

mat_file = append(basilisk_folder, sprintf("/vtu/matlab/moving_piston_timestep_%d.mat", timestep));
vtu_file = append(basilisk_folder, sprintf("/vtu/ascii/%06d_0.vtu", timestep));
%% convert the files from binary vtu to mat files
if isfile(vtu_file)
    disp("already made ascii file")
else
    disp("did not find ascii vtu file converting vtu file")
    command = append("/home/martin/ParaView-5.12.0-RC3-MPI-Linux-Python3.10-x86_64/bin/pvpython ~/Documents/master/basilisk/paraview_script.py ", basilisk_folder);
    re = system(command);
end

if isfile(mat_file)
    disp("already converted to mat file")
else
    %convert the vtu files to ascii
    %convert the vtu for the timestep to mat file
    command = sprintf("python3 /home/martin/Documents/master/basilisk/read_vtu_python.py %s %d", basilisk_folder, timestep);
    system(command);
end
load(mat_file)
disp("file loaded")
%% find amplitude using surface_probes.csv
surface_probes = load(append(basilisk_folder, "/surface_probes.csv"));
surface_probes  = surface_probes(2000:end,:);
%plot(surface_probes(:,2))
min_separation = 10;
LMax = islocalmax(surface_probes(:,2), 'MinSeparation', min_separation);
LMin = islocalmin(surface_probes(:,2), 'MinSeparation', min_separation);
a = (mean(surface_probes(LMax,2))-mean(surface_probes(LMin,2)))/2;


%% plot energy
%plot 5th order stokes alpha
figure;
[k, ~,~,~] = Stokes5th_alpha(a, omega,h,true);
hold on
% plot the energy from the basilisk results
x_start = 10; %where to plot the mean kinetic energy from and to
x_end = 12;


mask(X(:,:,1)<x_start)= 0;
mask(X(:,:,1)>x_end)= 0;
y = X(:,:,2);
u = U(:,:,1);
v = U(:,:,2);

g = 9.81;
alpha = omega/(a*g*k)*(u.^2 + v.^2).^.5;
mean_alpha = zeros(size(alpha, 1), 1)*NaN;
for i=1:size(alpha, 1)
    horizontal_slice_alpha = alpha(i,:);
    mean_alpha(i) = mean(horizontal_slice_alpha(mask(i,:)), 'omitnan');
end
crest_idx = find_crest_index_from_mask(y, mask);

mean_alpha = mean_alpha(mask(:,crest_idx));
crest_y = y(:,crest_idx);
crest_y = crest_y(mask(:,crest_idx));
plot(mean_alpha, crest_y/h, 'DisplayName',sprintf('ns t:%d', timestep))
legend()
title("mean energy")
xlabel('$\alpha=\frac{\omega}{agk}(u^2+v^2)^{\frac{1}{2}}$', 'interpreter', 'latex', 'FontSize', 20)
ylabel('$\frac{y}{h}$', 'interpreter', 'latex', 'FontSize', 20, 'rotation', 0)

%%
hold on
x_start = 8;
x_end = 10;
omega = 8.95;
h = 0.6;


basilisk_moving_piston_velocity_profile(mat_file,a, omega, h, x_start, x_end)