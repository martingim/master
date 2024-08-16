
figure;
hold on;
basilisk_folders = [];
titles = [];
%% run 5 
%test refineparam
% basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL14_refineparam1/"];     titles = [titles "LEVEL 14 refineparam 1"];
% basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL14_refineparam01/"];    titles = [titles "LEVEL 14 refineparam 0.1"];
% basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL14_refineparam001/"];   titles = [titles "LEVEL 14 refineparam 0.01"];
% basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL14_refineparam0001/"];  titles = [titles "run 5 LEVEL 14 refineparam 0.001"];

%test levels
% basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL15/"];  titles = [titles "run 5 LEVEL 15"];
% basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL14_refineparam0001/"];  titles = [titles "run 5 LEVEL 14 refineparam 0.001"];
% basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL13_refineparam0001/"];  titles = [titles "run 5 LEVEL 13 refineparam 0.001"];
% basilisk_folders = [basilisk_folders "~/Documents/results/run5/LEVEL12/"];  titles = [titles "run 5 LEVEL 12"];
%% run 1
%test levels
basilisk_folders = [basilisk_folders "~/Documents/results/run1/LEVEL12/"];  titles = [titles "run 1 LEVEL 12"];
basilisk_folders = [basilisk_folders "~/Documents/results/run1/LEVEL13/"];  titles = [titles "run 1 LEVEL 13"];
basilisk_folders = [basilisk_folders "~/Documents/results/run1/LEVEL13_2pad/"];  titles = [titles "run 1 LEVEL 13 2pad"];

%% run 4
%test levels
basilisk_folders = [basilisk_folders "~/Documents/results/run4/LEVEL12/"];  titles = [titles "run 4 LEVEL 12"];
basilisk_folders = [basilisk_folders "~/Documents/results/run4/LEVEL13/"];  titles = [titles "run 4 LEVEL 13"];
basilisk_folders = [basilisk_folders "~/Documents/results/run4/LEVEL14/"];  titles = [titles "run 4 LEVEL 14"];

%%

timestep = 25;
omega = 8.95;
h = 0.6;

mat_files = [];
vtu_files = [];
for i=1:size(basilisk_folders, 2)
    basilisk_folder = basilisk_folders(i);
    mat_files = [mat_files append(basilisk_folder, sprintf("/vtu/matlab/moving_piston_timestep_%d.mat", timestep))];
    vtu_files = [vtu_files append(basilisk_folder, sprintf("/vtu/ascii/%d_0.vtu", timestep))];
    vtu_file = vtu_files(i);
    mat_file = mat_files(i);
    %% convert the files from binary vtu to mat files
    if isfile(vtu_file)
        disp("already made ascii file")
    else
        disp("did not find ascii vtu file converting vtu file")
        command = sprintf("/home/martin/ParaView-5.12.0-RC3-MPI-Linux-Python3.10-x86_64/bin/pvpython ~/Documents/master/basilisk/paraview_script.py %s %d", basilisk_folder, timestep);
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

    % find amplitude using surface_probes.csv
    surface_probes = load(append(basilisk_folder, "/surface_probes.csv"));
    surface_probes = surface_probes(surface_probes(:,1)>timestep-1,:);
    surface_probes = surface_probes(surface_probes(:,1)<timestep+1,:);
    %plot(surface_probes(:,2))
    min_separation = 10;
    LMax = islocalmax(surface_probes(:,2), 'MinSeparation', min_separation);
    LMin = islocalmin(surface_probes(:,2), 'MinSeparation', min_separation);
    a = (mean(surface_probes(LMax,2))-mean(surface_probes(LMin,2)))/2;

    % plot energy
    %plot 5th order stokes alpha
    if i==1
        [k, ~,~,~] = Stokes5th_alpha(a, omega,h,true);
    else
        [k, ~,~,~] = Stokes5th_alpha(a, omega,h,false);
    end
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
    for j=1:size(alpha, 1)
        horizontal_slice_alpha = alpha(j,:);
        mean_alpha(j) = mean(horizontal_slice_alpha(mask(j,:)), 'omitnan');
    end
    crest_idx = find_crest_index_from_mask(y, mask);
    
    mean_alpha = mean_alpha(mask(:,crest_idx));
    crest_y = y(:,crest_idx);
    crest_y = crest_y(mask(:,crest_idx));
    plot(mean_alpha, crest_y/h, 'DisplayName',titles(i))
    legend()
    title("mean energy")
    xlabel('$\alpha=\frac{\omega}{agk}(u^2+v^2)^{\frac{1}{2}}$', 'interpreter', 'latex', 'FontSize', 20)
    ylabel('$\frac{y}{h}$', 'interpreter', 'latex', 'FontSize', 20, 'rotation', 0)
end
